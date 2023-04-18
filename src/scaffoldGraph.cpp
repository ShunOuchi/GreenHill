/*
Copyright (C) 2022 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GreenHill.

GreenHill is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GreenHill is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GreenHill; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "scaffoldGraph.h"
#include <omp.h>
#include <climits>
#include <cfloat>
#include <string>
#include <sstream>

using std::vector;
using std::string;
using std::cerr;
using std::endl;

//////////////////////////////////////////////////////////////////////////////////////
// const parameter define
//////////////////////////////////////////////////////////////////////////////////////
// this value must be 2^N (o <= N <= 64)
const unsigned ScaffoldGraph::TABLE_DIVID = 1024;
const double ScaffoldGraph::MAX_DIFF_RATE = 0.75;
const double ScaffoldGraph::EDGE_EXPECTED_RATE_TH = 0.5;
const double ScaffoldGraph::EDGE_EXPECTED_RATE_UPPER_TH = 4.0;
const double ScaffoldGraph::CHECK_USING_LONGER_LIB_TH = 0.1;
const unsigned ScaffoldGraph::SC_REP = 0x1;
const unsigned ScaffoldGraph::SC_INC = 0x2;
const unsigned ScaffoldGraph::SC_DEL = 0x4;
const unsigned ScaffoldGraph::SC_RR = 0x8;
const unsigned ScaffoldGraph::SC_RL = 0x10;
const double ScaffoldGraph::MAX_HOMO_RATE = 1.5;
const double ScaffoldGraph::MAX_HETERO_RATE = 0.75;
const double ScaffoldGraph::MAX_OVERLAP_IDENTITY_DIFF = 0.05;
const long ScaffoldGraph::MIN_NUM_MAPPED_TAG = 3;


//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
ScaffoldGraph::ScaffoldGraph()
: seedLength(0), minOverlap(0), hashOverlap(0), indexLength(0), minLink(0)
, tolerence(0) , minTolerenceFactor(0), numContig(0), numNode(0)
, averageCoverage(0), bubbleThreshold(0)
, contigFP(NULL), bubbleFP(NULL), overlapFP(NULL), graphLinkFP(0)
, contig(), contigUnlinkSet(), allLibraryMT(NULL), longReadLibraryMT(NULL), tagLibraryMT(NULL), HiCLibraryMT(NULL), targetLibraryIndex(0), coverage(), seqPool(), node() //edited by ouchi (HiCLibraryMT(NULL))
, numBubble() , contigPositionInScaffold(), overlapTable(TABLE_DIVID)
{
    bubbleFP = platanus::makeTemporaryFile();
    bubbleOpositeFP = platanus::makeTemporaryFile();
    deletedHeteroFP = platanus::makeTemporaryFile();
    remainingHeteroFP = platanus::makeTemporaryFile();
}

ScaffoldGraph::~ScaffoldGraph()
{
    fclose(bubbleFP);
    fclose(bubbleOpositeFP);
    fclose(deletedHeteroFP);
    fclose(remainingHeteroFP);
}


//////////////////////////////////////////////////////////////////////////////////////
// destroy graph
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::destroyGraph(void)
{
    node.clear();
}


//////////////////////////////////////////////////////////////////////////////////////
// find overlap between contigs
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::saveOverlap(const Mapper &map, const long hashOverlapValue, const long cutoffLength, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    long bufferSize;

    cerr << "saving overlaps... (LEN_CUTOFF=" << cutoffLength << ")" << endl;

    hashOverlap = hashOverlapValue;

    if (overlapFP != NULL)
        fclose(overlapFP);
    overlapFP = platanus::makeTemporaryFile();

    omp_lock_t lock[TABLE_DIVID];

	keyLength = map.getKeyLength();

    for (unsigned i = 0; i < TABLE_DIVID; ++i)
        omp_init_lock(&lock[i]);
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1) private(positionBuffer, bufferSize)
    for (long threadID = 0; threadID < numThread; ++threadID) {
        for (long i = threadID; i < numContig; i += numThread) {
            Kmer31 leftKmer(keyLength);
            Kmer31 rightKmer(keyLength);
            long k;
            bool isLeftUnknown = false;
            bool isRightUnknown = false;
            Overlap overlap;
            Overlap preOverlap;
            if (contig[i].length < cutoffLength) continue;

            // check contig terminal contains N
			if (leftKmer.setKmer(contig[i], 0, keyLength, keyLength))
				isLeftUnknown = true;
			if (rightKmer.setKmer(contig[i], contig[i].length - keyLength, keyLength, keyLength))
				isRightUnknown = true;

            if (!isLeftUnknown) {
                bufferSize = map.mapKey(leftKmer, positionBuffer);
                preOverlap.id1 = 0;
                for (long j = 0; j < bufferSize; ++j) {
                    // avoid double count from i to j and form j to i
                    // avoid ownself mapping too
                    if (abs(positionBuffer[j].id) - 1 <= i
                        || (positionBuffer[j].id - 1 == i && positionBuffer[j].offset == 0))
                        continue;

                    // kmer extension
                    if (positionBuffer[j].id > 0) {
                        overlap.length = contig[positionBuffer[j].id - 1].length - positionBuffer[j].offset;
                        if (overlap.length < minOverlap
                           || overlap.length > contig[i].length
                           || contig[positionBuffer[j].id - 1].length < cutoffLength)
                            continue;
                        for (k = keyLength; k < overlap.length; ++k) {
                            if (contig[i].base[k] != contig[positionBuffer[j].id - 1].base[positionBuffer[j].offset + k])
                                break;
                        }
                    }
                    else {
                        overlap.length = positionBuffer[j].offset + keyLength;
                        if (overlap.length < minOverlap
                            || overlap.length > contig[i].length
                            || contig[-(positionBuffer[j].id) - 1].length < cutoffLength)
                            continue;
                        for (k = keyLength; k < overlap.length; ++k) {
                            if (contig[i].base[k] != (0x3 ^ contig[-(positionBuffer[j].id) - 1].base[positionBuffer[j].offset + keyLength - k - 1]))
                                break;
                        }
                    }

                    if (k < overlap.length)
                        continue;

                    overlap.id1 = -(i + 1);
                    overlap.id2 = -(positionBuffer[j].id);
                    if (overlap.id1 == preOverlap.id1 && overlap.id2 == preOverlap.id2) {
                        if (overlap.length > preOverlap.length)
                            preOverlap.length = overlap.length;
                    }
                    else {
                        if (preOverlap.id1 != 0) {
                        unsigned mod = decideTableID(std::abs(preOverlap.id1));
                        omp_set_lock(&lock[mod]);
                        overlapTable[mod][std::pair<int, int>(preOverlap.id1, preOverlap.id2)] = preOverlap;
                        omp_unset_lock(&lock[mod]);
                        }
                        preOverlap = overlap;
                    }
                }
                if (preOverlap.id1 != 0) {
                        unsigned mod = decideTableID(std::abs(preOverlap.id1));
                        omp_set_lock(&lock[mod]);
                        overlapTable[mod][std::pair<int, int>(preOverlap.id1, preOverlap.id2)] = preOverlap;
                        omp_unset_lock(&lock[mod]);
                }
            }

            if (!isRightUnknown) {
                bufferSize = map.mapKey(rightKmer, positionBuffer);
                preOverlap.id1 = 0;
                for (int j = 0; j < bufferSize; ++j) {
                    if (abs(positionBuffer[j].id) - 1 <= i
                        || (positionBuffer[j].id - 1 == i && positionBuffer[j].offset == contig[i].length - keyLength))
                        continue;

                    if (positionBuffer[j].id > 0) {
                        overlap.length = positionBuffer[j].offset + keyLength;
                        if (overlap.length < minOverlap
                            || overlap.length > contig[i].length
                            || contig[positionBuffer[j].id - 1].length < cutoffLength)
                            continue;
                        for (k = keyLength; k < overlap.length; ++k) {
                            if (contig[i].base[contig[i].length - k - 1] != contig[positionBuffer[j].id - 1].base[positionBuffer[j].offset + keyLength -  k - 1])
                                break;
                        }
                    }
                    else {
                        overlap.length = contig[-(positionBuffer[j].id) - 1].length - positionBuffer[j].offset;
                        if (overlap.length < minOverlap
                            || overlap.length > contig[i].length
                            || contig[-(positionBuffer[j].id) - 1].length < cutoffLength)
                            continue;
                        for (k = keyLength; k < overlap.length; ++k) {
                            if (contig[i].base[contig[i].length - k - 1] != (0x3 ^ contig[-(positionBuffer[j].id) - 1].base[positionBuffer[j].offset + k]))
                                break;
                        }
                    }
                    if (k < overlap.length)
                        continue;

                    overlap.id1 = i + 1;
                    overlap.id2 = positionBuffer[j].id;

                    if (overlap.id1 == preOverlap.id1 && overlap.id2 == preOverlap.id2) {
                        if (overlap.length > preOverlap.length)
                            preOverlap.length = overlap.length;
                    }
                    else {
                        if (preOverlap.id1 != 0) {
                        unsigned mod = decideTableID(std::abs(preOverlap.id1));
                        omp_set_lock(&lock[mod]);
                        overlapTable[mod][std::pair<int, int>(preOverlap.id1, preOverlap.id2)] = preOverlap;
                        omp_unset_lock(&lock[mod]);
                        }
                        preOverlap = overlap;
                    }
                }
                if (preOverlap.id1 != 0) {
                        unsigned mod = decideTableID(std::abs(preOverlap.id1));
                        omp_set_lock(&lock[mod]);
                        overlapTable[mod][std::pair<int, int>(preOverlap.id1, preOverlap.id2)] = preOverlap;
                        omp_unset_lock(&lock[mod]);
                }
            }

        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// getter overlap length
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::getOverlap(long id1, long id2)
{
    long index;
    if (abs(id1) > abs(id2)) {
        index = id1;
        id1 = -id2;
        id2 = -index;
    }
    unsigned mod = decideTableID(abs(id1));

    std::pair<int, int> key(id1, id2);
    auto it = overlapTable[mod].find(key);
/*
    if (it != overlapTable[mod].end() && it->second.id1 == id1 && it->second.id2 == id2)
        return it->second.length;
	else if (minOverlap >= keyLength)
		return 0;
	else
        return getShortOverlap(id1, id2);
*/
    if (it != overlapTable[mod].end() && it->second.id1 == id1 && it->second.id2 == id2)
        return it->second.length;
	else
        return getShortOverlap(id1, id2);
}


//////////////////////////////////////////////////////////////////////////////////////
// getter overlap length
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::getShortOverlap(long id1, long id2) const
{
    long i, j;
    long contig1Add, contig2Add;
    long contig1Code, contig2Code;
    long contig1EXOR, contig2EXOR;

    if (id1 > 0) {
        contig1Add = contig[id1 - 1].length - 1;
        contig1Code = 1;
        contig1EXOR = 0x0;
    } else {
        contig1Add = 0;
        contig1Code = -1;
        contig1EXOR = 0x3;
    }

    if (id2 > 0) {
        contig2Add = -1;
        contig2Code = 1;
        contig2EXOR = 0x0;
    } else {
        contig2Add = contig[-id2 - 1].length;
        contig2Code = -1;
        contig2EXOR = 0x3;
    }

    for (i = hashOverlap; i >= minOverlap; --i) {
        for (j = 0; j < i; ++j) {
            if (((contig1EXOR) ^ (contig[id1 * contig1Code - 1].base[(contig1Code * -1 * j) + contig1Add])) != ((contig2EXOR) ^ (contig[id2 * contig2Code - 1].base[(contig2Code) * (i - j) + contig2Add]))) {
                break;
            }
        }
        if (i == j) {
            return i;
        }
    }

    return 0;

}


//////////////////////////////////////////////////////////////////////////////////////
// estimate the number of link (share kmer)
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::estimateLink(void)
{
    unsigned long long sum = 0;
    for (unsigned i = 0; i < contig.size(); ++i) {
        sum += contig[i].length;
    }
    return std::max(1.0, calcExpectedLink(averageCoverage, (double)sum / numContig, (double)sum / numContig, (double)((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize()) - 2.0 * (*allLibraryMT)[targetLibraryIndex][0].getAverageLength()));
}


//////////////////////////////////////////////////////////////////////////////////////
// getter overlap length
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::getScaffoldOverlap(long id1, long id2)
{
    if (id1 > 0)
        id1 = node[id1 - 1].contig[node[id1 - 1].numContig - 1].id;
    else
        id1 = -(node[-id1 - 1].contig[0].id);

    if (id2 > 0)
        id2 = node[id2 - 1].contig[0].id;
    else
        id2 = -(node[-id2 - 1].contig[node[-id2 - 1].numContig - 1].id);
    return getOverlap(id1, id2);
}


//////////////////////////////////////////////////////////////////////////////////////
// set init parameter
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::initScaffolding(std::vector<unsigned short> &cov, Mapper &mapper, const double ave, const double bubble)
{
    numContig = mapper.getNumSeq();
    mapper.moveSeq(this->contig);
    coverage = cov;
    numNode = numContig;
    averageCoverage = ave;
    bubbleThreshold = bubble;
    node.clear();
    node.resize(numNode);
    contigPositionInScaffold.clear();
    contigPositionInScaffold.resize(numContig);
    numBubble.clear();
    numBubble.resize(numContig, 0);
    contigState.clear();
	contigState.resize(numContig, 0);

    // # pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numNode; ++i) {
        node[i].numContig = 1;
        node[i].contig.resize(node[i].numContig);
        node[i].length = node[i].contig[0].end =this->contig[i].length;
        node[i].contig[0].id = i + 1;
        contigPositionInScaffold[i].resize(1);   //added by ouchi (contigPositionInScaffold)
        contigPositionInScaffold[i][0].id = i + 1; //edited by ouchi (contigPositionInScaffold)
        contigPositionInScaffold[i][0].offset = 0; //edited by ouchi (contigPositionInScaffold)
    }

}


//////////////////////////////////////////////////////////////////////////////////////
// count the number of bubble for each contig
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::countBubble(const platanus::Contig &bubble, const vector<platanus::Position> &bubblePosition)
{
    for (long i = 0; i < bubble.numSeq; ++i) {
        if (bubblePosition[i].id == 0) continue;
        long j = std::abs(bubblePosition[i].id) - 1;
        if ((static_cast<double>(std::abs(coverage[j] - bubble.coverage[i])) / coverage[j]) <= MAX_DIFF_RATE)
            ++(numBubble[j]);
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// classify node whether homo or hetero
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::classifyNode(void)
{
    const long MIN_HOMO_BUBBLE_NUM = 1;
    const long MIN_HOMO_COV = static_cast<long>(averageCoverage * 1.0 + 0.5);

    for (long i = 0; i < numNode; ++i) {
        long bubble = 0;
        for (long j = 0; j < node[i].numContig; ++j) {
            if (contigPositionInScaffold[std::abs(node[i].contig[j].id) - 1].size() == 1) //edited by ouchi (contigPositionInScaffold)
                bubble += numBubble[std::abs(node[i].contig[j].id) - 1];
        }

        if (bubble >= MIN_HOMO_BUBBLE_NUM || this->calcNodeCoverage(node[i]) >= MIN_HOMO_COV)
            node[i].isHomo = true;
        else
            node[i].isHomo = false;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
double ScaffoldGraph::calcNodeCoverage(const GraphNode &node)
{
    unsigned long long sum = 0;
    unsigned long long num = 0;
    for (long i = 0; i < node.numContig; ++i) {
        num += contig[abs(node.contig[i].id) - 1].length;
        sum += coverage[abs(node.contig[i].id) - 1] * contig[abs(node.contig[i].id) - 1].length;
    }

    if (num == 0)
        return 0.0;
    else
        return static_cast<double>(sum) / num;

}

void ScaffoldGraph::storeGraphLinkFromMappedPair(vector<GraphLink> &graphLinkPool, long numThread)
{
	if (allLibraryMT == NULL)
		return;

    long cutoffLength = std::max(tolerence, seedLength) * 2;

	for (long threadID = 0; threadID < numThread; ++threadID)
		rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

    # pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
        platanus::Position forwardResult;
        platanus::Position reverseResult;

        while (fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP)) {
            fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

            long i = abs(forwardResult.id) - 1;
            if (contigPositionInScaffold[i].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
            forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i][0].id : -(contigPositionInScaffold[i][0].id); //edited by ouchi (contigPositionInScaffold)
            forwardResult.offset = contigPositionInScaffold[i][0].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
            forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[i][0].offset].start; //edited by ouchi (contigPositionInScaffold)

            long j = abs(reverseResult.id) - 1;
            if (contigPositionInScaffold[j].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
            reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j][0].id : -(contigPositionInScaffold[j][0].id); //edited by ouchi (contigPositionInScaffold)
            reverseResult.offset = contigPositionInScaffold[j][0].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
            reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[j][0].offset].start; //edited by ouchi (contigPositionInScaffold)

            GraphLink graphLink;
            if (abs(forwardResult.id) < abs(reverseResult.id)) {
                graphLink.id1 = forwardResult.id;
                graphLink.offset1 = contigPositionInScaffold[i][0].offset; //edited by ouchi (contigPositionInScaffold)
                graphLink.id2 = -(reverseResult.id);
                graphLink.offset2 = contigPositionInScaffold[j][0].offset; //edited by ouchi (contigPositionInScaffold)
            }
            else {
                graphLink.id1 = reverseResult.id;
                graphLink.offset1 = contigPositionInScaffold[j][0].offset; //edited by ouchi (contigPositionInScaffold)
                graphLink.id2 = -(forwardResult.id);
                graphLink.offset2 = contigPositionInScaffold[i][0].offset; //edited by ouchi (contigPositionInScaffold)
            }

            // calc gap length (average length minus own node length plus offset)
            // image
            //                    ---------average gap----------
            // -----------------|node length |NNNNNNNNNNNNNNNNNN----------------
            //                  -- <- offset
            graphLink.gap = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
            long gapMinus;
            if (forwardResult.id > 0) {
                gapMinus = node[forwardResult.id - 1].length - forwardResult.offset;
            } else {
                gapMinus = forwardResult.offset + 1;
            }
            if (node[std::abs(forwardResult.id) - 1].length < cutoffLength) continue;
            graphLink.gap -= gapMinus;


            if (reverseResult.id > 0) {
                if (node[reverseResult.id - 1].length < cutoffLength || abs(forwardResult.id) == abs(reverseResult.id)) continue;
                graphLink.gap -= node[reverseResult.id-1].length - reverseResult.offset;
            } else {
                if (node[(-1 * reverseResult.id) - 1].length < cutoffLength || abs(forwardResult.id) == abs(reverseResult.id)) continue;
                graphLink.gap -= reverseResult.offset + 1;
            }

            // check contigs aren't overlapeed too large area
            if (-graphLink.gap > tolerence + this->getScaffoldOverlap(graphLink.id1, graphLink.id2)) continue;
            # pragma omp critical (push)
            {
                graphLinkPool.push_back(graphLink);
            }
        }
    }
}

void ScaffoldGraph::countMappedTagForEachContig(long numThread)
{
	if (tagLibraryMT == NULL)
		return;

    omp_set_num_threads(numThread);

	int forwardTagID;
	int reverseTagID;
	int forwardContigIndex;
	int reverseContigIndex;

	contigTagCounterUchar.resize(numThread);
	for (auto itr = contigTagCounterUchar.begin(); itr != contigTagCounterUchar.end(); ++itr)
		itr->clear();

	for (long threadID = 0; threadID < numThread; ++threadID)
		rewind((*tagLibraryMT)[threadID].mappedFP);
	for (long threadID = 0; threadID < numThread; ++threadID) {
		while (fread(&forwardTagID, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP)) {
			fread(&forwardContigIndex, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP);
			fread(&reverseTagID, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP);
			fread(&reverseContigIndex, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP);

			auto itr = contigTagCounterUchar[forwardTagID % numThread].find(std::make_pair(forwardTagID, forwardContigIndex));
			if (itr != contigTagCounterUchar[forwardTagID % numThread].end())
				itr->second = std::min(itr->second + 1, UCHAR_MAX - 1);
			else
				contigTagCounterUchar[forwardTagID % numThread][std::make_pair(forwardTagID, forwardContigIndex)] = 1;

			itr = contigTagCounterUchar[reverseTagID % numThread].find(std::make_pair(reverseTagID, reverseContigIndex));
			if (itr != contigTagCounterUchar[reverseTagID % numThread].end())
				itr->second = std::min(itr->second + 1, UCHAR_MAX - 1);
			else
				contigTagCounterUchar[reverseTagID % numThread][std::make_pair(reverseTagID, reverseContigIndex)] = 1;
		}
	}

	contigTagCounterUint.resize(numThread);
	for (auto itr = contigTagCounterUint.begin(); itr != contigTagCounterUint.end(); ++itr)
		itr->clear();

    # pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto itr = contigTagCounterUchar[threadID].begin(); itr != contigTagCounterUchar[threadID].end(); ++itr) {
			if (itr->second >= UCHAR_MAX - 1) {
				contigTagCounterUint[threadID][itr->first] = 0;
				contigTagCounterUchar[threadID].erase(itr);
			}
		}
	}

	for (long threadID = 0; threadID < numThread; ++threadID)
		rewind((*tagLibraryMT)[threadID].mappedFP);
    for (long threadID = 0; threadID < numThread; ++threadID) {
        while (fread(&forwardTagID, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP)) {
			fread(&forwardContigIndex, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP);
            fread(&reverseTagID, sizeof(platanus::Position), 1, (*tagLibraryMT)[threadID].mappedFP);
			fread(&reverseContigIndex, sizeof(int), 1, (*tagLibraryMT)[threadID].mappedFP);

			auto itr = contigTagCounterUint[forwardTagID % numThread].find(std::make_pair(forwardTagID, forwardContigIndex));
			if (itr != contigTagCounterUint[forwardTagID % numThread].end())
				itr->second = std::min(itr->second + 1, UINT_MAX - 1);

			itr = contigTagCounterUint[reverseTagID % numThread].find(std::make_pair(reverseTagID, reverseContigIndex));
			if (itr != contigTagCounterUint[reverseTagID % numThread].end())
				itr->second = std::min(itr->second + 1, UINT_MAX - 1);
        }
    }
}

void ScaffoldGraph::countMappedTagForEachScaffold(long numThread)
{
	if (tagLibraryMT == NULL)
		return;

    omp_set_num_threads(numThread);

	for (auto itr = this->node.begin(); itr != this->node.end(); ++itr)
		itr->numMappedTag.clear();

//	const unsigned minNumMapped = this->minLink;
	const unsigned minNumMapped = MIN_NUM_MAPPED_TAG;

	vector<std::unordered_map<std::pair<int, int>, unsigned char, platanus::PairHash, platanus::PairEqual> > nodeTagCounterUchar(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto contigItr = contigTagCounterUchar[threadID].begin(); contigItr != contigTagCounterUchar[threadID].end(); ++contigItr) {
			if (contigPositionInScaffold[contigItr->first.second].size() != 1) //edited by ouchi (contigPositionInScaffold)
				continue;

			int nodeIndex = abs(contigPositionInScaffold[contigItr->first.second][0].id) - 1; //edited by ouchi (contigPositionInScaffold)
			auto nodeItr = nodeTagCounterUchar[threadID].find(std::make_pair(contigItr->first.first, nodeIndex));
			if (nodeItr != nodeTagCounterUchar[threadID].end())
				nodeItr->second = std::min(static_cast<unsigned>(nodeItr->second + contigItr->second), static_cast<unsigned>(UCHAR_MAX - 1));
			else
				nodeTagCounterUchar[threadID][std::make_pair(contigItr->first.first, nodeIndex)] = contigItr->second;
		}
	}


    # pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto contigItr = contigTagCounterUint[threadID].begin(); contigItr != contigTagCounterUint[threadID].end(); ++contigItr) {
			if (contigPositionInScaffold[contigItr->first.second].size() != 1) //edited by ouchi (contigPositionInScaffold)
				continue;

			int nodeIndex = abs(contigPositionInScaffold[contigItr->first.second][0].id) - 1; //edited by ouchi (contigPositionInScaffold)
			auto nodeItr = nodeTagCounterUchar[threadID].find(std::make_pair(contigItr->first.first, nodeIndex));
			if (nodeItr != nodeTagCounterUchar[threadID].end())
				nodeItr->second = UCHAR_MAX - 1;
			else
				nodeTagCounterUchar[threadID][std::make_pair(contigItr->first.first, nodeIndex)] = UCHAR_MAX - 1;
		}
	}

	vector<std::unordered_map<std::pair<int, int>, unsigned int, platanus::PairHash, platanus::PairEqual> > nodeTagCounterUint(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto itr = nodeTagCounterUchar[threadID].begin(); itr != nodeTagCounterUchar[threadID].end(); ++itr) {
			if (itr->second >= UCHAR_MAX - 1) {
				nodeTagCounterUint[threadID][itr->first] = 0;
				nodeTagCounterUchar[threadID].erase(itr);
			}
			else if (itr->second < minNumMapped)
				nodeTagCounterUchar[threadID].erase(itr);
		}
	}


    # pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto contigItr = contigTagCounterUchar[threadID].begin(); contigItr != contigTagCounterUchar[threadID].end(); ++contigItr) {
			if (contigPositionInScaffold[contigItr->first.second].size() != 1) //edited by ouchi (contigPositionInScaffold)
				continue;

			int nodeIndex = abs(contigPositionInScaffold[contigItr->first.second][0].id) - 1; //edited by ouchi (contigPositionInScaffold)
			auto nodeItr = nodeTagCounterUint[threadID].find(std::make_pair(contigItr->first.first, nodeIndex));
			if (nodeItr != nodeTagCounterUint[threadID].end())
				nodeItr->second = std::min(static_cast<unsigned>(nodeItr->second + contigItr->second), static_cast<unsigned>(UINT_MAX - 1));
		}
	}

    # pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto contigItr = contigTagCounterUint[threadID].begin(); contigItr != contigTagCounterUint[threadID].end(); ++contigItr) {
			if (contigPositionInScaffold[contigItr->first.second].size() != 1) //edited by ouchi (contigPositionInScaffold)
				continue;

			int nodeIndex = abs(contigPositionInScaffold[contigItr->first.second][0].id) - 1; //edited by ouchi (contigPositionInScaffold)
			auto nodeItr = nodeTagCounterUint[threadID].find(std::make_pair(contigItr->first.first, nodeIndex));
			if (nodeItr != nodeTagCounterUint[threadID].end())
				nodeItr->second = std::min(static_cast<unsigned>(nodeItr->second + contigItr->second), static_cast<unsigned>(UINT_MAX - 1));
		}
	}

	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto itr = nodeTagCounterUchar[threadID].begin(); itr != nodeTagCounterUchar[threadID].end(); ++itr) {
			this->node[itr->first.second].numMappedTag[itr->first.first] = itr->second;
		}
	}
	nodeTagCounterUchar.clear();

	for (long threadID = 0; threadID < numThread; ++threadID) {
		for (auto itr = nodeTagCounterUint[threadID].begin(); itr != nodeTagCounterUint[threadID].end(); ++itr) {
			this->node[itr->first.second].numMappedTag[itr->first.first] = itr->second;
		}
	}
	nodeTagCounterUint.clear();
}

long ScaffoldGraph::getCommonTagBetweenNodePair(const long nodeID1, const long nodeID2)
{
	long sum = 0;

	GraphNode *node1 = &(this->node[abs(nodeID1) - 1]);
	GraphNode *node2 = &(this->node[abs(nodeID2) - 1]);

	if (node1->numMappedTag.size() < node2->numMappedTag.size())
		std::swap(node1, node2);

	for (auto itr1 = node1->numMappedTag.begin(); itr1 != node1->numMappedTag.end(); ++itr1) {
		auto itr2 = node2->numMappedTag.find(itr1->first);
		if (itr2 != node2->numMappedTag.end())
			sum += std::min(itr1->second, itr2->second);
	}

	return sum;
}

//////////////////////////////////////////////////////////////////////////////////////
// calculate the number of link
// link means that the number of read which forward and reverse mappend each contig
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::calcLink(const long libraryIndex, const long linkThreshold, const long numThread)
{
	graphLinkFP.resize((*allLibraryMT).size(), NULL);
    if (graphLinkFP[libraryIndex] != NULL)
        fclose(graphLinkFP[libraryIndex]);
    graphLinkFP[libraryIndex] = platanus::makeTemporaryFile();

//    cerr << "linking scaffolds (MIN_LINK = " << minLink << ")" << endl;
    cerr << "linking scaffolds..." << endl;

    vector<GraphLink> graphLinkPool;

	storeGraphLinkFromMappedPair(graphLinkPool, numThread);

    long totalLink = graphLinkPool.size();
    cerr << "sorting links in contigID order..." << endl;
    sort(graphLinkPool.begin(), graphLinkPool.end());

    graphLinkPool.emplace_back();
    graphLinkPool[totalLink].id1 = 0;

    cerr << "estimating contig distances..." << endl;
    vector<GraphLinkPoolIndex> indexes(1, 0);
    for (long idx = 1; idx < totalLink; ++idx) {
        if (graphLinkPool[idx - 1].id1 != graphLinkPool[idx].id1 || graphLinkPool[idx - 1].id2 != graphLinkPool[idx].id2) {
            if (indexes.back().numLink < linkThreshold) {
                indexes.pop_back();
            }
            indexes.emplace_back(idx);
        }
        ++indexes.back().numLink;
    }
    if (indexes.back().numLink < linkThreshold) {
        indexes.pop_back();
    }
    sort(indexes.begin(), indexes.end(), GraphLinkPoolIndexGreater());

	if (numThread == 1) {
		for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
			calcLinkAndWriteGraphLinkFile(graphLinkPool, indexes[idxIndex], this->targetLibraryIndex, false);
		}
	}
	else {
		# pragma omp parallel for schedule(dynamic)
		for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
			calcLinkAndWriteGraphLinkFile(graphLinkPool, indexes[idxIndex], this->targetLibraryIndex, true);
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::calcLinkAndWriteGraphLinkFile(const std::vector<GraphLink>& graphLinkPool, const GraphLinkPoolIndex& index, const long libraryIndex, const bool multiThreadFlag)
{
    unsigned long indexBegin = index.index;
    long numLink = index.numLink;

    GraphLink graphLink = graphLinkPool[indexBegin];
    vector<long> breakdown1(node[abs(graphLink.id1) - 1].numContig, 0);
    vector<long> breakdown2(node[abs(graphLink.id2) - 1].numContig, 0);
    for (long idxLink = indexBegin, endLink = idxLink + numLink; idxLink < endLink; ++idxLink) {
        ++(breakdown1[graphLinkPool[idxLink].offset1]);
        ++(breakdown2[graphLinkPool[idxLink].offset2]);
    }

//    graphLink.gap = estimateGapSize(graphLinkPool, indexBegin, numLink);
    graphLink.gap = estimateGapSizeAverage(graphLinkPool, indexBegin, numLink);
    if (graphLink.gap == LONG_MIN) return;

    long breakdownSize = breakdown1.size() + breakdown2.size();
    std::unique_ptr<long[]> tmp(new long[breakdownSize]());
    std::copy(breakdown1.begin(), breakdown1.end(), tmp.get());
    std::copy(breakdown2.begin(), breakdown2.end(), &(tmp[breakdown1.size()]));

	if (multiThreadFlag) {
		#pragma omp critical (write_link)
		{
			fwrite(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex]);
			fwrite(&graphLink, sizeof(GraphLink), 1, graphLinkFP[libraryIndex]);
			fwrite(tmp.get(), sizeof(long), breakdownSize, graphLinkFP[libraryIndex]);
		}
	}
	else {
		fwrite(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex]);
		fwrite(&graphLink, sizeof(GraphLink), 1, graphLinkFP[libraryIndex]);
		fwrite(tmp.get(), sizeof(long), breakdownSize, graphLinkFP[libraryIndex]);
	}
}


//////////////////////////////////////////////////////////////////////////////////////
// estimate gap size
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::estimateGapSize(const std::vector<GraphLink>& links, const unsigned index, const unsigned size) const
{
    platanus::NormalDistribution<unsigned> distribution((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize(), (*allLibraryMT)[targetLibraryIndex][0].getSDInsSize());
    long distance = 0;
//    const long distanceMedian = links[index + size / 2].gap % 2 && links[index + size / 2].gap ?
//            (links[index + size / 2 - 1].gap + links[index + size / 2].gap) / 2 : links[index + size / 2].gap;
    const long distanceMedian = size % 2 ? links[index + size / 2].gap : (links[index + size / 2 - 1].gap + links[index + size / 2].gap) / 2;
    const long distanceWidth = std::min(tolerence, (*allLibraryMT)[targetLibraryIndex][0].getSDInsSize() * 3);
    const long distanceShort = std::max(links[index].gap, distanceMedian - distanceWidth);
    const long distanceLong  = std::min(links[index + size - 1].gap, distanceMedian + distanceWidth) + 1;

    double bestLikelihood = -DBL_MAX;
    for (long tmpDistance = distanceShort; tmpDistance < distanceLong; ++tmpDistance) {
        double likelihood = 0;
        for (unsigned idx = 0; idx < size; ++idx) {
            const long tmpInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize() + tmpDistance - links[index + idx].gap;
            likelihood += log(distribution(tmpInsSize));
        }
        if (likelihood > bestLikelihood) {
            distance = tmpDistance;
            bestLikelihood = likelihood;
        }
    }
    if (bestLikelihood <= -DBL_MAX) {
        return LONG_MIN;
    }
    return distance;
}


//////////////////////////////////////////////////////////////////////////////////////
// estimate gap size
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::estimateGapSizeAverage(const std::vector<GraphLink>& links, const unsigned index, const unsigned size) const
{
    double gap = 0;
//    long tmpNumLink = 0.9 * size;
//    long startIndex = 0.05 * tmpNumLink + index;
//    for (long idx = 0; idx < tmpNumLink; ++idx) {
//        gap += links[startIndex + idx].gap;
//    }
//    return gap / tmpNumLink + 0.5;

    for (long i = 0; i < size; ++i)
        gap += links[i].gap;

    return gap / size + 0.5;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::makeGraph(const long numThread)
{
    calcLink(this->targetLibraryIndex, this->minLink, numThread);

    cerr << "constructing scaffold graph" << endl;

    for (int i = 0; i < numNode; ++i) {
		node[i].edge.clear();
		node[i].numEdge = 0;
	}

    rewind(graphLinkFP[this->targetLibraryIndex]);
    long numLink;
    while (fread(&numLink, sizeof(long), 1, graphLinkFP[this->targetLibraryIndex])) {
        GraphLink link;
        fread(&link, sizeof(GraphLink), 1, graphLinkFP[this->targetLibraryIndex]);

        long i = abs(link.id1) - 1;
        long j = abs(link.id2) - 1;
        long m = node[i].numEdge;
        long n = node[j].numEdge;

		node[i].edge.resize(m + 1);
		node[j].edge.resize(n + 1);

        node[i].edge[m].numLink = numLink;
        node[j].edge[n].numLink = numLink;

        for (unsigned int k = 0; k < node[i].numContig; ++k) {
            fread(&(node[i].edge[m].breakdown[k]), sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        }
        for (unsigned int k = 0; k < node[j].numContig; ++k) {
            fread(&(node[j].edge[n].breakdown[k]), sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        }

        node[i].edge[m].length = link.gap;
        node[j].edge[n].length = link.gap;
        node[i].edge[m].direction = static_cast<char>(link.id1 / (i + 1));
        node[j].edge[n].direction = static_cast<char>(-(link.id2) / (j + 1));
        if (link.id1 * link.id2 > 0) {
            node[i].edge[m].end = j + 1;
            node[j].edge[n].end = i + 1;
        }
        else {
            node[i].edge[m].end = -(j + 1);
            node[j].edge[n].end = -(i + 1);
        }
        ++(node[i].numEdge);
        ++(node[j].numEdge);
    }
    fclose(graphLinkFP[this->targetLibraryIndex]);
    graphLinkFP[this->targetLibraryIndex] = NULL;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        GraphNode &tmp = node[nodeID];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);
    }
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
bool ScaffoldGraph::checkDeleteEdge(const GraphEdge &edge1, const GraphEdge &edge2, const GraphNode &node1, const GraphNode &node2)
{

    if (edge1.direction * edge2.direction < 0
        || edge1.length + node1.length <= edge2.length
        || edge2.length + node2.length <= edge1.length) {
        return 0;
    }

    if (edge1.direction > 0) {
        if (abs(edge1.length + node1.length - edge2.length) <= tolerence + this->getScaffoldOverlap(edge1.end, edge2.end)
         || abs(edge2.length + node2.length - edge1.length) <= tolerence + this->getScaffoldOverlap(edge2.end, edge1.end))
            return 0;
    }
    else {
        if (abs(edge1.length + node1.length - edge2.length) <= tolerence + this->getScaffoldOverlap(edge2.end, edge1.end)
         || abs(edge2.length + node2.length - edge1.length) <= tolerence + this->getScaffoldOverlap(edge1.end, edge2.end))
            return 0;
    }

    return 1;
}


void ScaffoldGraph::deleteThinEdge(void)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

//    std::cerr << "removing thin edges (NUM_LINK < min(" << minLink << ", EXPECTED_NUM_LINK))" <<  " ..." << std::endl;
    std::cerr << "removing thin edges (NUM_LINK < EXPECTED_NUM_LINK)" <<  " ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            GraphNode &node2 = this->node[abs(node1.edge[edgeID].end) - 1];

            const double currentLibraryCoverage1 = (*allLibraryMT)[targetLibraryIndex][0].getAverageCoverage() * this->calcNodeCoverage(node1) / this->averageCoverage;
            const double currentLibraryCoverage2 = (*allLibraryMT)[targetLibraryIndex][0].getAverageCoverage() * this->calcNodeCoverage(node2) / this->averageCoverage;
            const long expectedNumLink = this->calcExpectedLink(std::min(currentLibraryCoverage1, currentLibraryCoverage2), node1.length, node2.length, std::max(static_cast<double>((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize()), 0.0)) + 0.5;

//            if (node1.edge[edgeID].numLink < std::min(minLink, expectedNumLink)) {
            if (node1.edge[edgeID].numLink < std::round(expectedNumLink / 2.0)) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

void ScaffoldGraph::deleteThinEdgeConstant(const long linkThreshold)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing thin edges (NUM_LINK < " <<  linkThreshold << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (node1.edge[edgeID].numLink < linkThreshold) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::deleteErroneousEdge(const long numThread)
{
    vector<long> ids;
    long numDelete = 0;

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2)) continue;

                double errorRate1 = static_cast<double>(edge1.numLink)
                        / std::max(1.0, calcExpectedLinkNode(node[nodeID], node1, edge1.length));
                if (errorRate1 > EDGE_EXPECTED_RATE_UPPER_TH) continue;
                double errorRate2 = static_cast<double>(edge2.numLink)
                        / std::max(1.0, calcExpectedLinkNode(node[nodeID], node2, edge2.length));
                if (errorRate2 > EDGE_EXPECTED_RATE_UPPER_TH) continue;
                if (edge1.numLink < edge2.numLink && errorRate1 / errorRate2 <= EDGE_EXPECTED_RATE_TH) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (edge2.numLink < edge1.numLink && errorRate2 / errorRate1 <= EDGE_EXPECTED_RATE_TH) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge2.end);
                    }
                    ++numDelete;
                }
            }
        }
    }
    this->deleteEdges(ids);
    return numDelete;
}

long ScaffoldGraph::deleteErroneousEdgeNumLinkRate(const long numThread)
{
    vector<long> ids;
    long numDelete = 0;
	const double RATE_THRESHOLD = 0.125;

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2)) continue;

                if (edge1.numLink < RATE_THRESHOLD * edge2.numLink) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (edge2.numLink < RATE_THRESHOLD * edge1.numLink) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge2.end);
                    }
                    ++numDelete;
                }
            }
        }
    }
    this->deleteEdges(ids);
    return numDelete;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::crushBubbleIterative(const double bubbleThreshold, const double averageCoverage, const long numThread)
{
    unsigned long long numCrush;
    unsigned long long totalCrush = 0;
    if (bubbleThreshold == 0) return;

    cerr << "removing bubbles... (MAX_BUBBLE_IDENTITY = " << bubbleThreshold << ")" << endl;
    do {
        numCrush = this->crushBubble(bubbleThreshold, averageCoverage, numThread);
        totalCrush += numCrush;
        cerr << "NUM_REMOVED_BUBBLES=" << numCrush << endl;
    } while (numCrush > 0);

    cerr << "TOTAL_NUM_REMOVED_BUBBLES=" << totalCrush << endl;

}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::deleteErroneousEdgeIterative(const long numThread)
{
    cerr << "removing erroneous edges..." << endl << endl;;
    long totalDelete = 0;
    long numDelete;
    do {
        numDelete = this->deleteErroneousEdge(numThread);
        totalDelete += numDelete;
        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
    } while (numDelete > 0);

    cerr << "TOTAL_NUM_REMOVED_EDGES =" << totalDelete << endl << endl;
}

void ScaffoldGraph::deleteErroneousEdgeNumLinkRateIterative(const long numThread)
{
    cerr << "removing erroneous edges..." << endl << endl;;
    long totalDelete = 0;
    long numDelete;
    do {
        numDelete = this->deleteErroneousEdgeNumLinkRate(numThread);
        totalDelete += numDelete;
        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
    } while (numDelete > 0);

    cerr << "TOTAL_NUM_REMOVED_EDGES =" << totalDelete << endl << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::deleteRepeatEdge(void)
{
    vector<long> ids;

    cerr << "deleting edges from repeat contigs..." << endl;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        if (node[nodeID].numContig == 1) continue;
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2)) continue;

                for (long m = 0; m < node[nodeID].numContig; ++m) {
                    if (edge1.breakdown[m] < minLink || edge2.breakdown[m] < minLink) continue;

                    for (long n = 0; n < node[nodeID].numEdge; ++n) {
                        node[nodeID].edge[n].numLink -= node[nodeID].edge[n].breakdown[m];
                        node[nodeID].edge[n].breakdown[m] = 0;
                    }
                    contigPositionInScaffold[abs(node[nodeID].contig[m].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
                }
            }
        }
    }

    for (long i = 0; i < numNode; ++i) {
        for (long j = 0; j < node[i].numEdge; ++j) {
            if (node[i].edge[j].numLink < minLink) {
                ids.push_back(i + 1);
                ids.push_back(node[i].edge[j].end);
            }
        }
    }
    this->deleteEdges(ids);

}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::detectRepeat(const double averageCoverage)
{
    const double coverageThreshold = averageCoverage * 1.5;

    // # pragma omp for schedule(dynamic)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        if (node[nodeID].numContig == 1 && calcNodeCoverage(node[nodeID]) > coverageThreshold) {
            node[nodeID].state |= SC_REP;
            continue;
        }
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                long edgeEnd1, edgeEnd2;
                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);

                if (edge1.direction * edge2.direction < 0) continue;
                GraphNode &node1 = (node[abs(edge1.end) - 1]);
                if (edge1.length + node1.length <= edge2.length) continue;
                GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (edge2.length + node2.length <= edge1.length) continue;

                if (edge1.isForward()) {
                    edgeEnd1 = edge1.end;
                    edgeEnd2 = edge2.end;
                } else {
                    edgeEnd1 = edge2.end;
                    edgeEnd2 = edge1.end;
                }
                if (abs(edge1.length + node1.length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2.length + node2.length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1)) continue;

                node[nodeID].state |= SC_REP;
                edgeID1 = node[nodeID].numEdge;
                break;
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::deleteEdges(vector<long> &ids)
{
    long id1, id2;
    long idsSize = ids.size();
    GraphNode *nodePointer;
    for (long i = 0; i < idsSize; i += 2) {
        id1 = ids[i];
        id2 = ids[i + 1];
		if (id2 == 0)
			continue;

        nodePointer = &(node[id1 - 1]);
        for (long j = 0; j < nodePointer->numEdge; ++j) {
            if (nodePointer->edge[j].end == id2) {
                nodePointer->edge[j] = nodePointer->edge[nodePointer->numEdge - 1];
                --(nodePointer->numEdge);
                break;
            }
        }

        if (id2 > 0) {
            nodePointer = &(node[id2 - 1]);
            for (long j = 0; j < nodePointer->numEdge; ++j) {
                if (nodePointer->edge[j].end == id1) {
                    nodePointer->edge[j] = nodePointer->edge[nodePointer->numEdge - 1];
                    --(nodePointer->numEdge);
                    break;
                }
            }
        }
        else {
            nodePointer = &(node[-id2 - 1]);
            for (long j = 0; j < nodePointer->numEdge; ++j) {
                if (nodePointer->edge[j].end == -id1) {
                    nodePointer->edge[j] = nodePointer->edge[nodePointer->numEdge - 1];
                    --(nodePointer->numEdge);
                    break;
                }
            }
        }
    }
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        GraphNode &tmp = node[nodeID];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);
    }
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::split(void)
{
    long numNewNode = 0;
    long newContigPoolSize = 0;
    bool isSplit = false;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
    vector<char> breakpoint;
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        if (node[nodeID].numContig == 1) {
            fwrite(&(node[nodeID].numContig), sizeof(long), 1, scaffoldFP);
            fwrite(&(node[nodeID].contig[0]), sizeof(ScaffoldPart), 1, scaffoldFP);
            ++numNewNode;
            ++newContigPoolSize;
            continue;
        }
        isSplit = false;
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2)) continue;

                breakpoint.resize(node[nodeID].numContig + 1, 0);
                long m = 0;
                while (m < node[nodeID].numContig) {
                    while (m < node[nodeID].numContig && edge1.breakdown[m] < minLink && edge1.breakdown[m] < minLink)
                        ++m;
                    breakpoint[m] = 1;
                    while (m < node[nodeID].numContig && edge1.breakdown[m] >= minLink)
                        ++m;
                    breakpoint[m] = 1;
                    while (m < node[nodeID].numContig && edge2.breakdown[m] >= minLink)
                        ++m;
                    breakpoint[m] = 1;
                }
                isSplit = true;
            }
        }
        if (isSplit) {
            long j = 0;
            while (j < node[nodeID].numContig) {
                long start = node[nodeID].contig[j].start;
                long k = j;
                while (breakpoint[j + 1] == 0) {
                    node[nodeID].contig[j].start -= start;
                    node[nodeID].contig[j].end -= start;
                    ++j;
                }
                node[nodeID].contig[j].start -= start;
                node[nodeID].contig[j].end -= start;
                ++j;
                long tmp = j - k;
                fwrite(&tmp, sizeof(long), 1, scaffoldFP);
                for (tmp = k; tmp < j; ++tmp)
                    fwrite(&(node[nodeID].contig[tmp]), sizeof(ScaffoldPart), 1, scaffoldFP);
                ++numNewNode;
                newContigPoolSize += tmp;
            }
        } else {
            fwrite(&(node[nodeID].numContig), sizeof(long), 1, scaffoldFP);
            auto contigIterator = node[nodeID].contig.begin();
            auto contigEnd = node[nodeID].contig.end();
            for (; contigIterator != contigEnd; ++contigIterator)
                fwrite(&(*contigIterator), sizeof(ScaffoldPart), 1, scaffoldFP);
            ++numNewNode;
            newContigPoolSize += node[nodeID].numContig;
        }

    }

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::makeScaffold(void)
{
    long numNewContig = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    vector<GraphLayout> include, candidate;
    include.reserve(10000);
    candidate.reserve(10000);

    cerr << "scaffolding..." << endl;

    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long i = 0; i < numNode; ++i) {
        if (node[i].state & (SC_INC | SC_REP | SC_DEL)) continue;
        include.clear();
        candidate.clear();
        include.push_back(GraphLayout());
        include[0].start =  include[0].distance = 0;
        include[0].end = node[i].length;
        include[0].id = i + 1;
        numNewContig = node[i].numContig;
        node[i].state |= SC_INC;
        for (long j = 0; j < node[i].numEdge; ++j) {
            long tmpNodeID = abs(node[i].edge[j].end) - 1;
            if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

            candidate.push_back(GraphLayout());
            unsigned long long candidateEnd = candidate.size() - 1;
            if (node[i].edge[j].direction > 0) {
                candidate[candidateEnd].start = include[0].end + node[i].edge[j].length;
                candidate[candidateEnd].end = candidate[candidateEnd].start + node[tmpNodeID].length;
            } else {
                candidate[candidateEnd].end = -(node[i].edge[j].length);
                candidate[candidateEnd].start = candidate[candidateEnd].end - node[tmpNodeID].length;
            }
            candidate[candidateEnd].id = node[i].edge[j].end;
            candidate[candidateEnd].distance = 1;
            candidate[candidateEnd].numLink = node[i].edge[j].numLink;
        }


        while (candidate.size() > 0) {
            long minDistance = candidate[0].distance;
            long closest = candidate[0].start;
//            long maxCoverage = calcNodeCoverage(node[id2Index(candidate[0].id)]);
            long minCandidateID = 0;
            for (unsigned j = 1; j < candidate.size(); ++j) {
                if (candidate[j].distance < minDistance || (candidate[j].distance == minDistance && abs(candidate[j].start) < closest)) {
//                if (candidate[j].distance < minDistance  || (candidate[j].distance == minDistance && calcNodeCoverage(node[id2Index(candidate[j].id)]) > maxCoverage)) {
                    minDistance = candidate[j].distance;
//					maxCoverage = calcNodeCoverage(node[id2Index(candidate[j].id)]);
                    closest = abs(candidate[j].start);
                    minCandidateID = j;
                }
            }
/*
        while (candidate.size() > 0) {
            long minDistance = candidate[0].distance;
            long maxNumLink = candidate[0].numLink;
            long minCandidateID = 0;
            for (unsigned j = 1; j < candidate.size(); ++j) {
                if (candidate[j].distance < minDistance || (candidate[j].distance == minDistance && candidate[j].numLink > maxNumLink)) {
                    minDistance = candidate[j].distance;
                    maxNumLink = candidate[j].numLink;
                    minCandidateID = j;
                }
            }
*/


            long tmpNodeID = abs(candidate[minCandidateID].id) - 1;

            if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) {
                auto candidateIterator = candidate.begin();
                candidateIterator = candidate.erase(candidateIterator + minCandidateID);
                continue;
            }
            unsigned j = 0;
            for (; j < include.size(); ++j) {
				long overlapTolerence = std::min(tolerence, std::min(candidate[minCandidateID].end - candidate[minCandidateID].start, include[j].end - include[j].start) / 2);

/*
				if (candidate[minCandidateID].end <= include[j].start || candidate[minCandidateID].start >= include[j].end)
					continue;

				long overlapLength = -(calcGapSizeFromPositions(candidate[minCandidateID].start, candidate[minCandidateID].end, include[j].start, include[j].end));

				if (candidate[minCandidateID].start - include[j].end - (include[j].start - candidate[minCandidateID].end) < 0) {
				 	if (overlapLength <= overlapTolerence + this->getScaffoldOverlap(candidate[minCandidateID].id, include[j].id))
						continue;
				}
				else {
				 	if (overlapLength <= overlapTolerence + this->getScaffoldOverlap(include[j].id, candidate[minCandidateID].id))
						continue;
				}
*/
				if (candidate[minCandidateID].end <= include[j].start
				 || candidate[minCandidateID].start >= include[j].end
				 || abs(candidate[minCandidateID].start - include[j].end) <= overlapTolerence + this->getScaffoldOverlap(include[j].id, candidate[minCandidateID].id)
				 || abs(candidate[minCandidateID].end - include[j].start) <= overlapTolerence + this->getScaffoldOverlap(candidate[minCandidateID].id, include[j].id))
					continue;

                break;
            }
            if (j == include.size()) {
                include.push_back(candidate[minCandidateID]);

                GraphNode &newNode = node[abs(include[include.size() - 1].id) - 1];
                if (!(newNode.state & SC_REP)) {
                    for (long k = 0; k < newNode.numEdge; ++k) {
                        long tmpNodeID = abs(newNode.edge[k].end) - 1;
                        if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

                        GraphLayout tmpLayout;
                        if (include[include.size() - 1].id * newNode.edge[k].direction > 0) {
                            tmpLayout.start = include[include.size() - 1].end + newNode.edge[k].length;
                            tmpLayout.end = tmpLayout.start + node[tmpNodeID].length;
                        }
                        else {
                            tmpLayout.end = include[include.size() - 1].start - newNode.edge[k].length;
                            tmpLayout.start = tmpLayout.end - node[tmpNodeID].length;
                        }
                        tmpLayout.id = include[include.size() - 1].id > 0 ? newNode.edge[k].end : -(newNode.edge[k].end);
                        tmpLayout.distance = include[include.size() - 1].distance + 1;
                        tmpLayout.numLink = newNode.edge[k].numLink;
                        candidate.push_back(tmpLayout);
                    }
                }

                numNewContig += newNode.numContig;
                if (!(newNode.state & SC_REP))
                    newNode.state |= SC_INC;
            }

            auto candidateIterator = candidate.begin();
            candidateIterator = candidate.erase(candidateIterator + minCandidateID);
        }

        sort(include.begin(), include.end());
        long includeSizeMover = include.size();
        long j = 0;
        for (; node[abs(include[j].id) - 1].state & SC_REP; ++j)
            numNewContig -= node[abs(include[j].id) - 1].numContig;
        for (; node[abs(include[includeSizeMover - 1].id) - 1].state & SC_REP; --includeSizeMover)
            numNewContig -= node[abs(include[includeSizeMover - 1].id) - 1].numContig;
        fwrite(&numNewContig, sizeof(long), 1, scaffoldFP);
        long minStart = include[j].start;
        for (; j < includeSizeMover; ++j) {
            long tmpNodeID = abs(include[j].id) - 1;

            node[tmpNodeID].state |= SC_INC;

            include[j].start -= minStart;
            include[j].end -= minStart;

            if (include[j].start != 0) {
                long overlapLength = this->getScaffoldOverlap(include[j-1].id, include[j].id);
                if (overlapLength > 0 && overlapLength + include[j].start - include[j-1].end <= tolerence) {
                    overlapLength = include[j-1].end - include[j].start - overlapLength;
                    for (long k = j; k < includeSizeMover; ++k) {
                        include[k].end += overlapLength;
                        include[k].start += overlapLength;
                    }
                }
            }

            if (include[j].id > 0) {
                for (long k = 0; k < node[tmpNodeID].numContig; ++k) {
                    ScaffoldPart tmpScaffoldPart(node[tmpNodeID].contig[k].id, include[j].start + node[tmpNodeID].contig[k].start, include[j].start + node[tmpNodeID].contig[k].end);
                    fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                }
            }
            else {
                for (long k = node[tmpNodeID].numContig - 1; k >= 0; --k) {
                    ScaffoldPart tmpScaffoldPart(-(node[tmpNodeID].contig[k].id), include[j].start + node[tmpNodeID].length - node[tmpNodeID].contig[k].end, include[j].start + node[tmpNodeID].length - node[tmpNodeID].contig[k].start);
                    fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                }
            }
        }

        newContigPoolSize += numContig;
        ++numNewNode;
    }
    include.clear();
    candidate.clear();

    for (long i = 0; i < numNode; ++i) {
        if (!(node[i].state & SC_REP)) continue;
        long numContig = node[i].numContig;
        if (node[i].state & SC_INC) {
            for (long j = 0; j < numContig; ++j)
                contigPositionInScaffold[abs(node[i].contig[j].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
        } else {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            for (long j = 0; j < numContig; ++j)
                fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
            newContigPoolSize += node[i].numContig;
            ++numNewNode;
        }
    }
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::remake(const long numNewNode, const long newContigPoolSize, FILE *scaffoldFP)
{
    this->destroyGraph();
    numNode = numNewNode;
    node.resize(numNode);
    vector<long> numOccurrence(numContig, 0);

    //added by ouchi (contigPositionInScaffold)
    for (long i = 0; i < numContig; ++i) {
        contigPositionInScaffold[i].clear();
    }
    //

    rewind(scaffoldFP);
    for (long i = 0; i < numNode; ++i) {
        fread(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
        node[i].contig.resize(node[i].numContig);
        for (long j = 0; j < node[i].numContig; ++j)
            fread(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);

//        node[i].length = node[i].contig[node[i].numContig - 1].end;

		long start = node[i].contig.front().start;
		long end = node[i].contig.back().end;
		for (long j = 0; j < node[i].numContig; ++j) {
			if (start > node[i].contig[j].start)
				start = node[i].contig[j].start;
			if (end < node[i].contig[j].end)
				end = node[i].contig[j].end;
		}
		if (start != 0) {
			for (long j = 0; j < node[i].numContig; ++j) {
				node[i].contig[j].start -= start;
				node[i].contig[j].end -= start;
			}
		}
        node[i].length = end - start;

        node[i].state = 0;
        node[i].isHomo = false;
		node[i].oppositeBubbleNodeID = 0;
		node[i].numMappedTag.clear();
        for (long j = 0; j < node[i].numContig; ++j) {
            long tmp = abs(node[i].contig[j].id) - 1;
            ++numOccurrence[tmp];
//            contigPositionInScaffold[tmp].offset = j; //deleted by ouchi (contigPositionInScaffold)

            if (node[i].contig[j].id > 0)
                contigPositionInScaffold[tmp].emplace_back(j, i + 1); //edited by ouchi (contigPositionInScaffold)
            else
                contigPositionInScaffold[tmp].emplace_back(j, -(i + 1)); //edited by ouchi (contigPositionInScaffold)
        }
    }

    //deleted by ouchi (contigPositionInScaffold)
    /*for (long i = 0; i < numContig; ++i) {
        if (numOccurrence[i] != 1)
            contigPositionInScaffold[i].id = 0;
    }*/
    //
    //added by ouchi
    /*vector<long> deleteNode;
    cerr << "detect repeat contig..." << endl;
    long numRepeatContig = 0;
    for (long i = 0; i < numContig; ++i) {
        if (contigPositionInScaffold[i].size() > 1)
            ++numRepeatContig;
        for (long j = 1; j < contigPositionInScaffold[i].size(); ++j) {
            for (long k = 0; k < j; ++k) {
                if (!(node[id2Index(contigPositionInScaffold[i][k].id)].state & SC_DEL) &&
                    abs(contigPositionInScaffold[i][j].id) != abs(contigPositionInScaffold[i][k].id)) {
                    if (node[id2Index(contigPositionInScaffold[i][j].id)] == node[id2Index(contigPositionInScaffold[i][k].id)]) {
                        node[id2Index(contigPositionInScaffold[i][j].id)].state |= SC_DEL;
                        deleteNode.push_back(id2Index(contigPositionInScaffold[i][j].id));
                    }
                }
            }
        }
    }
    cerr << "TOTAL_NUM_REPEAT_CONTIG=" << numRepeatContig << endl;
    if (deleteNode.size() > 0) {
        cerr << "delete redundant node..." << endl;
        std::sort(deleteNode.begin(), deleteNode.end(), std::greater<long>());
        for (long i = 0; i < deleteNode.size(); ++i) {
            node.erase(node.begin() + deleteNode[i]);
            --numNode;
        }
        for (long i = 0; i < numContig; ++i) {
            contigPositionInScaffold[i].clear();
        }
        for (long i = 0; i < numNode; ++i) {
            for (long j = 0; j < node[i].numContig; ++j) {
                long tmp = abs(node[i].contig[j].id) - 1;
                if (node[i].contig[j].id > 0)
                    contigPositionInScaffold[tmp].emplace_back(j, i + 1);
                else
                    contigPositionInScaffold[tmp].emplace_back(j, -(i + 1));
            }
        }
        cerr << "TOTAL_NUM_DELETE=" << deleteNode.size() << endl;
    }
    // */

    this->classifyNode();
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long ScaffoldGraph::crushHeteroBubble(const double averageCoverage)
{
    long numCrush = 0;
    const double homoCoverageThreshold = averageCoverage * MAX_HOMO_RATE + 0.5;
    const double heteroCoverageThreshold = averageCoverage * MAX_HETERO_RATE + 0.5;
    vector<long> ids;
    vector<GraphLayout> work, layout1, layout2;
    this->detectRepeat(averageCoverage);

    if (bubbleThreshold == 0.0) return 0;

	if (bubbleFP == NULL)
		bubbleFP = platanus::makeTemporaryFile();

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = node[nodeID].edge[edgeID1];
                const GraphEdge &edge2 = node[nodeID].edge[edgeID2];

                work.clear();
                layout1.clear();
                layout2.clear();

                if (edge1.direction * edge2.direction < 0) continue;
                GraphNode *node1 = &(node[abs(edge1.end) - 1]);
                if ((node1->state & SC_DEL) || edge1.length + node1->length <= edge2.length) continue;
                GraphNode *node2 = &(node[abs(edge2.end) - 1]);
                if ((node2->state & SC_DEL) || edge2.length + node2->length <= edge1.length) continue;
                if (node1->isHomo && node2->isHomo) continue;
                long edgeEnd1, edgeEnd2;
                if (edge1.isForward()) {
                    edgeEnd1 = edge1.end;
                    edgeEnd2 = edge2.end;
                } else {
                    edgeEnd1 = edge2.end;
                    edgeEnd2 = edge1.end;
                }

                if (abs(edge1.length + node1->length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2.length + node2->length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1)) continue;

                this->layoutNodes(node1, layout1, work);
                this->layoutNodes(node2, layout2, work);
                long rightEdge = std::min(layout1.size(), layout2.size());
                long leftEdge = rightEdge;
                long layoutID = 0;
                for (; layoutID < leftEdge; ++layoutID) {
                    if (layout1[layoutID].id != layout2[layoutID].id)
                        break;
                }
                if (layoutID == 0 || layoutID == leftEdge) continue;
                leftEdge = layoutID - 1;

                if (calcNodeCoverage(node[abs(layout1[leftEdge].id) - 1]) >= homoCoverageThreshold) continue;
                for (layoutID = 1; layoutID <= rightEdge; ++layoutID) {
                    if (layout1[layout1.size() - layoutID].id != layout2[layout2.size() - layoutID].id)
                        break;
                }
                if (layoutID == 1) continue;
                rightEdge = layoutID - 1;

				if (abs(layout1[leftEdge].id) == abs(layout1[layout1.size() - rightEdge].id)) continue;

                if (calcNodeCoverage(node[abs(layout1[layout1.size() - rightEdge].id) - 1]) > homoCoverageThreshold) continue;
                double coverage1 = this->layoutAverageCoverage(layout1, leftEdge + 1, layout1.size() - rightEdge - leftEdge - 1);
                double coverage2 = this->layoutAverageCoverage(layout2, leftEdge + 1, layout2.size() - rightEdge - leftEdge - 1);
                const vector<GraphLayout> &layoutRef = coverage1 < coverage2 ? layout1 : layout2;
                if (rightEdge + leftEdge + 1 >= static_cast<long>(layoutRef.size()) || coverage1 > heteroCoverageThreshold || coverage2 > heteroCoverageThreshold) continue;

                for (layoutID = leftEdge + 1; layoutID < static_cast<long>(layoutRef.size()) - rightEdge; ++layoutID) {
                    node1 = &(node[abs(layoutRef[layoutID].id) - 1]);
                    for (long n = 0; n < node1->numEdge; ++n) {
                        ids.push_back(static_cast<long>(abs(layoutRef[layoutID].id) - 1) + 1);
                        ids.push_back(node1->edge[n].end);
                    }
                    for (long n = 0; n < node1->numContig; ++n)
                        contigPositionInScaffold[abs(node1->contig[n].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
                    node1->state |= SC_DEL;
                }
                vector<char> scaffoldSeq;

                this->layout2seq(layoutRef, leftEdge + 1, layoutRef.size() - rightEdge - leftEdge, scaffoldSeq);
                long scaffoldSize = scaffoldSeq.size();
                fwrite(&(scaffoldSize), sizeof(long), 1, bubbleFP);
                auto scaffoldIterator = scaffoldSeq.begin();
                auto scaffoldEnd = scaffoldSeq.end();
                for (; scaffoldIterator != scaffoldEnd; ++scaffoldIterator)
                    fwrite(&(*scaffoldIterator), sizeof(char), 1, bubbleFP);
                fwrite(coverage1 < coverage2 ? &coverage1 : &coverage2, sizeof(double), 1, bubbleFP);
                ++numCrush;
            }
        }
    }

    this->deleteEdges(ids);

    for (long i = 0; i < numNode; ++i)
        node[i].state &= ~SC_REP;

    cerr << "NUM_REMOVED_BUBBLES=" << numCrush << " (COVERAGE_THRESHOLD)" << endl;
    return numCrush;


}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::cutAndPrintSeq(const long minSeqLength, const unsigned long long readLength, const string &outFilename, const string &componentFilename)
{
    std::ofstream out(outFilename);
    std::ofstream com(componentFilename);
    vector<long> maxLeftOverlap(numContig, 0);
    vector<long> maxRightOverlap(numContig, 0);

    long minOverlap = this->minOverlap;
    this->minOverlap = 1;

    long maxNumContig = 0;
    long numScaffold = 0;
    long scaffoldLength = 0;
    for (unsigned i = 0; i < TABLE_DIVID; ++i) {
		auto overlapIterator = overlapTable[i].begin();
		auto overlapEnd = overlapTable[i].end();
		for (; overlapIterator != overlapEnd; ++overlapIterator) {
			if (overlapIterator->second.id1 > 0) {
				if (overlapIterator->second.length > maxRightOverlap[overlapIterator->second.id1 - 1])
					maxRightOverlap[overlapIterator->second.id1 - 1] = overlapIterator->second.length;
			} else if (overlapIterator->second.length > maxLeftOverlap[-(overlapIterator->second.id1) - 1])
				maxLeftOverlap[-(overlapIterator->second.id1) - 1] = overlapIterator->second.length;

			if (overlapIterator->second.id2 > 0) {
				if (overlapIterator->second.length > maxLeftOverlap[overlapIterator->second.id2 - 1])
					maxLeftOverlap[overlapIterator->second.id2 - 1] = overlapIterator->second.length;
			} else if (overlapIterator->second.length > maxRightOverlap[-(overlapIterator->second.id2) - 1])
				maxRightOverlap[-(overlapIterator->second.id2) - 1] = overlapIterator->second.length;
		}
    }
    for (long i = 0; i < numNode; ++i) {
        if (node[i].numContig > maxNumContig)
            maxNumContig = node[i].numContig;
    }
    vector<long> leftCut(maxNumContig, 0);
    vector<long> rightCut(maxNumContig, 0);
    vector<long> gap(maxNumContig, 0);

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        long contigID = 0;
        for (; contigID < node[nodeID].numContig; ++contigID) {
            if (contigPositionInScaffold[abs(node[nodeID].contig[contigID].id) - 1].size() == 1) //edited by ouchi (contigPositionInScaffold)
                break;
        }
        if (contigID == node[nodeID].numContig) continue;

        scaffoldLength = 0;
        if (node[nodeID].contig[0].id > 0)
            leftCut[0] = maxLeftOverlap[node[nodeID].contig[0].id - 1] / 2;
        else
            leftCut[0] = maxRightOverlap[-(node[nodeID].contig[0].id) - 1] / 2;

        for (contigID = 0; contigID < node[nodeID].numContig - 1; ++contigID) {
/*
            if (node[nodeID].contig[contigID].end > node[nodeID].contig[contigID + 1].start && this->getOverlap(node[nodeID].contig[contigID].id, node[nodeID].contig[contigID + 1].id) > minOverlap) {
                rightCut[contigID] = 0;
                leftCut[contigID + 1] = node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start;
                gap[contigID] = 0;
            } else if (node[nodeID].contig[contigID].end > node[nodeID].contig[contigID + 1].start) {
                rightCut[contigID] = leftCut[contigID + 1] = (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start) / 2;
                gap[contigID] = node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start;
            } else {
                rightCut[contigID] = leftCut[contigID + 1] = this->getSimilarOverlap(node[nodeID].contig[contigID].id, node[nodeID].contig[contigID + 1].id) / 2;
                gap[contigID] = node[nodeID].contig[contigID + 1].start - node[nodeID].contig[contigID].end + 2 * rightCut[contigID];
            }
*/
            if (node[nodeID].contig[contigID].end > node[nodeID].contig[contigID + 1].start) {
                rightCut[contigID] = leftCut[contigID + 1] = (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start) / 2;
                gap[contigID] = node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start;
            } else {
                rightCut[contigID] = leftCut[contigID + 1] = this->getSimilarOverlap(node[nodeID].contig[contigID].id, node[nodeID].contig[contigID + 1].id) / 2;
                gap[contigID] = node[nodeID].contig[contigID + 1].start - node[nodeID].contig[contigID].end + 2 * rightCut[contigID];
            }
            scaffoldLength += (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID].start) - leftCut[contigID] - rightCut[contigID] + gap[contigID];
        }

        if (node[nodeID].contig[contigID].id > 0)
            rightCut[contigID] = maxRightOverlap[node[nodeID].contig[contigID].id - 1] / 2;
        else
            rightCut[contigID] = maxLeftOverlap[-(node[nodeID].contig[contigID].id) - 1] / 2;

        scaffoldLength += (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID].start) - leftCut[contigID] - rightCut[contigID];
        gap[contigID] = 0;

        if (scaffoldLength < minSeqLength) continue;

        ++numScaffold;
        if (node[nodeID].numContig > 1) {
            out << ">scaffold" << numScaffold << "_len" << scaffoldLength << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeID]) + 0.5) << "_read" << readLength << endl;
            com << "scaffold" << numScaffold << "_len" << scaffoldLength << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeID]) + 0.5);
        } else {
            out << ">scaffold" << numScaffold << "_len" << scaffoldLength << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeID]) + 0.5) << "_read" << readLength << "_single" << endl;
            com << "scaffold" << numScaffold << "_len" << scaffoldLength << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeID]) + 0.5) << "_single";
        }
        long lineLength = 0;
        for (contigID = 0; contigID < node[nodeID].numContig; ++contigID) {
            if (node[nodeID].contig[contigID].id > 0) {
                for (long seqID = leftCut[contigID]; seqID < contig[node[nodeID].contig[contigID].id - 1].length - rightCut[contigID]; ++seqID) {
                    out.put(platanus::Bin2Char(contig[node[nodeID].contig[contigID].id - 1].base[seqID]));
                    ++lineLength;
                    if (lineLength % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                        out.put('\n');
                }
            } else {
                for (long seqID = contig[-(node[nodeID].contig[contigID].id) - 1].length - leftCut[contigID] - 1; seqID >= (long)rightCut[contigID]; --seqID) {
                    if (contig[-(node[nodeID].contig[contigID].id) - 1].base[seqID] != 4)
                        out.put(platanus::Bin2Char(0x3 ^ contig[-(node[nodeID].contig[contigID].id) - 1].base[seqID]));
                    else
                        out.put('N');
                    ++lineLength;
                    if (lineLength % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                        out.put('\n');
                }
            }

            for (long k = 0; k < gap[contigID]; ++k) {
                out.put('N');
                ++lineLength;
                if (lineLength % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                    out.put('\n');
            }
        }
        if (lineLength % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');

        for (long k = 0; k < node[nodeID].numContig; ++k)
            com << "\t" << node[nodeID].contig[k].id;
        com << std::endl;
    }

    out.close();
    com.close();

    this->minOverlap = minOverlap;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long ScaffoldGraph::crushBubble(const double bubbleThreshold, const double averageCoverage, const long numThread)
{
    long numCrush = 0;
    vector<long> ids;
    vector<GraphLayout> work, layout1, layout2;
    long re[4] = {0};

	if (bubbleFP == NULL)
		bubbleFP = platanus::makeTemporaryFile();

    this->detectRepeat(averageCoverage);
    omp_lock_t lock;
    omp_init_lock(&lock);
    omp_set_num_threads(numThread);
    // # pragma omp parallel for schedule(dynamic) private(work, layout1, layout2) reduction(+: numCrush)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = node[nodeID].edge[edgeID1];
                const GraphEdge &edge2 = node[nodeID].edge[edgeID2];

                work.clear();
                layout1.clear();
                layout2.clear();


                if (edge1.direction * edge2.direction < 0) continue;
                GraphNode *node1 = &(node[abs(edge1.end) - 1]);
                if ((node1->state & SC_DEL) || edge1.length + node1->length <= edge2.length) continue;
                GraphNode *node2 = &(node[abs(edge2.end) - 1]);
                if ((node2->state & SC_DEL) || edge2.length + node2->length <= edge1.length) continue;
                if (node1->isHomo && node2->isHomo) continue;

                long edgeEnd1, edgeEnd2;
                if (edge1.isForward()) {
                    edgeEnd1 = edge1.end;
                    edgeEnd2 = edge2.end;
                } else {
                    edgeEnd1 = edge2.end;
                    edgeEnd2 = edge1.end;
                }

                if (abs(edge1.length + node1->length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2.length + node2->length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1)) continue;

                // omp_set_lock(&lock);
                this->layoutNodes(node1, layout1, work);
                this->layoutNodes(node2, layout2, work);
                // omp_unset_lock(&lock);
                long rightEdge = std::min(layout1.size(), layout2.size());
                long leftEdge = rightEdge;
                long layoutID = 0;
                for (; layoutID < leftEdge; ++layoutID) {
                    if (layout1[layoutID].id != layout2[layoutID].id)
                        break;
                }
                if (layoutID == 0 || layoutID == leftEdge) {++re[3]; continue;}
                leftEdge = layoutID - 1;

                for (layoutID = 1; layoutID <= rightEdge; ++layoutID) {
                    if (layout1[layout1.size() - layoutID].id != layout2[layout2.size() - layoutID].id)
                        break;
                }
                if (layoutID == 1) {++re[3]; continue;}
                rightEdge = layoutID - 1;

				if (abs(layout1[leftEdge].id) == abs(layout1[layout1.size() - rightEdge].id)) continue;

                double coverage1 = this->layoutAverageCoverage(layout1, leftEdge + 1, layout1.size() - rightEdge - leftEdge - 1);
                double coverage2 = this->layoutAverageCoverage(layout2, leftEdge + 1, layout2.size() - rightEdge - leftEdge - 1);
                const double coverageThreshold = averageCoverage * 2.0;
                const vector<GraphLayout> &layoutRef = coverage1 < coverage2 ? layout1 : layout2;
                if (rightEdge + leftEdge + 1 >= (long)layoutRef.size() || coverage1 + coverage2 > coverageThreshold)
                {++re[0];
                continue;}

                vector<char> scaffoldSeq1;
                vector<char> scaffoldSeq2;
                vector<long> alternativeMatrix;

                this->layout2seq(layout1, leftEdge + 1, layout1.size() - rightEdge - leftEdge - 1, scaffoldSeq1);
                this->layout2seq(layout2, leftEdge + 1, layout2.size() - rightEdge - leftEdge - 1, scaffoldSeq2);
                long sizeDiff = scaffoldSeq1.size() - scaffoldSeq2.size();
                if (std::abs(sizeDiff) > std::max(scaffoldSeq1.size(), scaffoldSeq2.size()) * bubbleThreshold) {re[1]++; continue;}
                if (this->alignScaffold(scaffoldSeq1, scaffoldSeq2, alternativeMatrix, std::max(scaffoldSeq1.size(), scaffoldSeq2.size()) * bubbleThreshold) > std::max(scaffoldSeq1.size(), scaffoldSeq2.size()) * bubbleThreshold) {++re[2]; continue;}

                omp_set_lock(&lock);
                for (layoutID = leftEdge + 1; layoutID < (long)layoutRef.size() - rightEdge; ++layoutID) {
                    node1 = &(node[abs(layoutRef[layoutID].id) - 1]);
                    for (long n = 0; n < node1->numEdge; ++n) {
                        ids.push_back((long)(abs(layoutRef[layoutID].id) - 1) + 1);
                        ids.push_back(node1->edge[n].end);
                    }
                    for (long n = 0; n < node1->numContig; ++n)
                        contigPositionInScaffold[abs(node1->contig[n].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
                    node1->state |= SC_DEL;
                }

                this->layout2seq(layoutRef, leftEdge + 1, layoutRef.size() - rightEdge - leftEdge, scaffoldSeq1);
                long scaffoldSize = scaffoldSeq1.size();
                fwrite(&(scaffoldSize), sizeof(long), 1, bubbleFP);
                auto scaffoldIterator = scaffoldSeq1.begin();
                auto scaffoldEnd = scaffoldSeq1.end();
                for (; scaffoldIterator != scaffoldEnd; ++scaffoldIterator)
                    fwrite(&(*scaffoldIterator), sizeof(char), 1, bubbleFP);
                fwrite(coverage1 < coverage2 ? &coverage1 : &coverage2, sizeof(double), 1, bubbleFP);

                ++numCrush;
                omp_unset_lock(&lock);
            }
        }
    }
    omp_destroy_lock(&lock);
    this->deleteEdges(ids);

    for (long i = 0; i < numNode; ++i)
        node[i].state &= ~SC_REP;

    return numCrush;


}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::layoutNodes(GraphNode *newNode, vector<GraphLayout> &ret, vector<GraphLayout> &work)
{
    long retID;
    ret.resize(1);
    work.resize(0);
    ret[0].start =  ret[0].distance = 0;
    ret[0].end = newNode->length;
    ret[0].id = static_cast<long>(newNode - &(*(node.begin()))) + 1;
    newNode->state |= SC_INC;

    for (long edgeID = 0; edgeID < newNode->numEdge; ++edgeID) {
        long tmpNodeID = abs(newNode->edge[edgeID].end) - 1;
        if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;
        work.resize(work.size() + 1);
        long workSize = work.size();
        if (newNode->edge[edgeID].isForward()) {
            work[workSize - 1].start = ret[0].end + newNode->edge[edgeID].length;
            work[workSize - 1].end   = work[workSize - 1].start + node[tmpNodeID].length;
        }
        else {
            work[workSize - 1].end = -(newNode->edge[edgeID].length);
            work[workSize - 1].start = work[workSize - 1].end - node[tmpNodeID].length;
        }
        work[workSize - 1].id = newNode->edge[edgeID].end;
        work[workSize - 1].distance = 1;
        work[workSize - 1].numLink = newNode->edge[edgeID].numLink;
    }

    while (work.size() > 0) {
        long minDistance = work[0].distance;
        long minNumLink = work[0].start;
        long minDistanceID = 0;
        long workSize = work.size();
        for (long workID = 1; workID < workSize; ++workID) {
            if (work[workID].distance < minDistance) {
                minDistance = work[workID].distance;
                minNumLink = work[workID].start;
                minDistanceID = workID;
            } else if (work[workID].distance == minDistance && work[workID].start < minNumLink) {
                minNumLink = work[workID].start;
                minDistanceID = workID;
            }
        }

        long tmpNodeID = abs(work[minDistanceID].id) - 1;
        if ((node[tmpNodeID].state & SC_INC) && !((node[tmpNodeID].state & SC_INC) & SC_REP)) {
            auto workIterator = work.begin();
            workIterator = work.erase(workIterator + minDistanceID);
            continue;
        }

        for (retID = 0; retID < (long)ret.size(); ++retID) {
            if (work[minDistanceID].end <= ret[retID].start
                || work[minDistanceID].start >= ret[retID].end
                || abs(work[minDistanceID].start - ret[retID].end) <= tolerence + this->getScaffoldOverlap(ret[retID].id, work[minDistanceID].id)
                || abs(work[minDistanceID].end - ret[retID].start) <= tolerence + this->getScaffoldOverlap(work[minDistanceID].id, ret[retID].id))
                continue;
            break;
        }

        if (retID == static_cast<long>(ret.size())) {
            ret.resize(ret.size() + 1);
            ret[ret.size() - 1] = work[minDistanceID];

            newNode = &(node[abs(ret[ret.size() - 1].id) - 1]);
            if (~(newNode->state) & SC_REP) {
                for (long edgeID = 0; edgeID < newNode->numEdge; ++edgeID) {
                    long tmpNodeID = abs(newNode->edge[edgeID].end) - 1;
                    if ((node[tmpNodeID].state & SC_INC) && !(node[tmpNodeID].state & SC_REP)) continue;

                    work.resize(work.size() + 1);
                    workSize = work.size();
                    if (ret[ret.size() - 1].id * newNode->edge[edgeID].direction > 0) {
                        work[workSize - 1].start = ret[ret.size() - 1].end + newNode->edge[edgeID].length;
                        work[workSize - 1].end   = work[workSize - 1].start + node[tmpNodeID].length;
                    } else {
                        work[workSize - 1].end   = ret[ret.size() - 1].start - newNode->edge[edgeID].length;
                        work[workSize - 1].start = work[workSize - 1].end - node[tmpNodeID].length;
                    }
                    work[workSize - 1].id = ret[ret.size() - 1].id > 0 ? newNode->edge[edgeID].end : -(newNode->edge[edgeID].end);
                    work[workSize - 1].distance = ret[ret.size() - 1].distance + 1;
                    work[workSize - 1].numLink = newNode->edge[edgeID].numLink;
                }
            }

            if (!(newNode->state & SC_REP))
                newNode->state |= SC_INC;
        }

        auto workIterator = work.begin();
        workIterator = work.erase(workIterator + minDistanceID);
    }
    sort(ret.begin(), ret.end());

    long retSize = ret.size();
    for (retID = 1; retID < retSize; ++retID) {
        node[abs(ret[retID].id) - 1].state &= ~SC_INC;
        ret[retID].start -= ret[0].start;
        ret[retID].end   -= ret[0].start;
        if (ret[retID].start != 0) {
            long tmp = this->getScaffoldOverlap(ret[retID-1].id, ret[retID].id);
            if (tmp + ret[retID].start - ret[retID-1].end <= tolerence) {
                tmp = ret[retID-1].end - ret[retID].start - tmp;
                for (long retIDs = retID; retIDs < retSize; ++retIDs) {
                    ret[retIDs].end += tmp;
                    ret[retIDs].start += tmp;
                }
            } else if (ret[retID].start < ret[retID-1].end) {
                tmp = ret[retID-1].end - ret[retID].start + 1;
                for (long retIDs = retID; retIDs < retSize; ++retIDs) {
                    ret[retIDs].end += tmp;
                    ret[retIDs].start += tmp;
                }
            }
        }
    }
    node[abs(ret[0].id) - 1].state &= ~SC_INC;
    ret[0].end -= ret[0].start;
    ret[0].start = 0;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::layout2seq(const vector<GraphLayout> &leftOverlap, const long startPoint, const long numLeftOverlap, vector<char> &ret)
{
    GraphNode *tmpNodePointer;
    long k;
    ret.clear();
    platanus::SEQ *tmpContigPointer;
    for (long i = 0; i < numLeftOverlap; ++i) {
        if (leftOverlap[startPoint + i].id > 0) {
            tmpNodePointer = &(node[leftOverlap[startPoint + i].id - 1]);
            for (long j = 0; j < tmpNodePointer->numContig; ++j) {
                if (i == 0) {
                    k = 0;
                } else if (j == 0) {
                    k = leftOverlap[startPoint + i - 1].end - leftOverlap[startPoint + i].start;
                } else {
                    k = tmpNodePointer->contig[j - 1].end - tmpNodePointer->contig[j].start;
                }
                for (; k < 0; ++k)
                    ret.push_back(4);

                if (tmpNodePointer->contig[j].id > 0) {
                    tmpContigPointer = &(contig[tmpNodePointer->contig[j].id - 1]);
                    if (k < tmpContigPointer->length) {
                        ret.insert(ret.end(), tmpContigPointer->base.begin() + k, tmpContigPointer->base.begin() + tmpContigPointer->length);
                    }
                }
                else {
                    tmpContigPointer = &(contig[-(tmpNodePointer->contig[j].id) - 1]);
                    for (; k < tmpContigPointer->length; ++k) {
                        if (tmpContigPointer->base[tmpContigPointer->length - k - 1] >= 4)
                            ret.push_back(4);
                        else
                            ret.push_back(0x3 ^ (tmpContigPointer->base[tmpContigPointer->length - k - 1]));
                    }
                }
            }
        }
        else {
            tmpNodePointer = &(node[-(leftOverlap[startPoint + i].id) - 1]);
            for (long j = tmpNodePointer->numContig - 1; j >= 0; --j) {

                if (i == 0)
                    k = 0;
                else if (j == tmpNodePointer->numContig - 1)
                    k = leftOverlap[startPoint + i - 1].end - leftOverlap[startPoint + i].start;
                else
                    k = tmpNodePointer->contig[j].end - tmpNodePointer->contig[j + 1].start;

                for (; k < 0; ++k) {
                    ret.push_back(4);
                }

                if (tmpNodePointer->contig[j].id > 0) {
                    tmpContigPointer = &(contig[tmpNodePointer->contig[j].id - 1]);
                    for (; k < tmpContigPointer->length; ++k) {
                        if (tmpContigPointer->base[tmpContigPointer->length - k - 1] >= 4)
                            ret.push_back(4);
                        else
                            ret.push_back(0x3 ^ (tmpContigPointer->base[tmpContigPointer->length - k - 1]));
                    }
                }
                else {
                    tmpContigPointer = &(contig[-(tmpNodePointer->contig[j].id) - 1]);
                    if (k < tmpContigPointer->length) {
                        ret.insert(ret.end(), tmpContigPointer->base.begin() + k, tmpContigPointer->base.begin() + tmpContigPointer->length);
                    }
                }
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::alignScaffold(const vector<char> &scaffold1, const vector<char> &scaffold2, vector<long> &work, const long scoreThreshold) const
{
    long workSize = 2 * scoreThreshold + 1;
    work.resize(workSize * 2 + 1);
    long m, n;
    long min = 0;


    for (n = 0; n < scoreThreshold; ++n) {
        work[n] = scoreThreshold + 1;
    }
    for (n = scoreThreshold; n < workSize; ++n) {
        work[n] = n - scoreThreshold;
    }

    for (m = 0; m < static_cast<long>(scaffold1.size()); ++m) {

        if (m - scoreThreshold > 0 && scaffold1[m] == scaffold2[m - scoreThreshold - 1]) {
            work[(~m & 1) * workSize] = work[(m & 1) * workSize];
        } else {
            work[(~m & 1) * workSize] = work[(m & 1) * workSize] + 1;
            if (work[(~m & 1) * workSize] > work[(m & 1) * workSize + 1] + 1) {
                work[(~m & 1) * workSize] = work[(m & 1) * workSize + 1] + 1;
            }
        }
        min = work[(~m & 1) * workSize];

        for (n = 0; n < workSize - 1; ++n) {
            if (m - scoreThreshold + n >= 0 && m - scoreThreshold + n < static_cast<long>(scaffold2.size()) && scaffold1[m] == scaffold2[m - scoreThreshold + n]) {
                work[(~m & 1) * workSize + n + 1] = work[(m & 1) * workSize + n + 1];
            } else {
                work[(~m & 1) * workSize + n + 1] = work[(m & 1) * workSize + n + 1] + 1;
                if (work[(~m & 1) * workSize + n + 1] > work[(m & 1) * workSize + n + 2] + 1) {
                    work[(~m & 1) * workSize + n + 1] = work[(m & 1) * workSize + n + 2] + 1;
                }
                if (work[(~m & 1) * workSize + n + 1] > work[(~m & 1) * workSize + n] + 1) {
                    work[(~m & 1) * workSize + n + 1] = work[(~m & 1) * workSize + n] + 1;
                }
            }
            if (work[(~m & 1) * workSize + n + 1] < min) {
                min = work[(~m & 1) * workSize + n + 1];
            }
        }

        if (m - scoreThreshold + n >= 0 && m - scoreThreshold + n < static_cast<long>(scaffold2.size()) && scaffold1[m] == scaffold2[m - scoreThreshold + n]) {
            work[(~m & 1) * workSize + n + 1] = work[(m & 1) * workSize + n + 1];
        } else {
            work[(~m & 1) * workSize + n + 1] = work[(m & 1) * workSize + n + 1] + 1;
            if (work[(~m & 1) * workSize + n + 1] > work[(~m & 1) * workSize + n] + 1) {
                work[(~m & 1) * workSize + n + 1] = work[(~m & 1) * workSize + n] + 1;
            }
        }
        if (work[(~m & 1) * workSize + n + 1] < min) {
            min = work[(~m & 1) * workSize + n + 1];
        }

        if (min > scoreThreshold) return min;

    }
    return static_cast<long>(work[(m & 1) * workSize + scoreThreshold]);
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
double ScaffoldGraph::layoutAverageCoverage(const vector<GraphLayout> &leftOverlap, const long startPoint, const long leftOverlapSize) const
{
    long sum = 0;
    long num = 0;

    for (long i = 0; i < leftOverlapSize; ++i) {
        const GraphNode &nodeRef = node[abs(leftOverlap[i + startPoint].id) - 1];
        for (long j = 0; j < nodeRef.numContig; ++j) {
            num += contig[abs(nodeRef.contig[j].id) - 1].length;
            sum += coverage[abs(nodeRef.contig[j].id) - 1] * contig[abs(nodeRef.contig[j].id) - 1].length;
        }
    }

    if (num != 0)
        return static_cast<double>(sum) / num;
    else
        return 0.0;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::deleteHeteroEdge(void)
{
    long numDelete = 0;
    const double homoCoverageThreshold = static_cast<long>(averageCoverage * MAX_HOMO_RATE + 0.5);
    const double heteroCoverageThreshold = static_cast<long>(averageCoverage * MAX_HETERO_RATE + 0.5);
    vector<long> ids;
    if (bubbleThreshold == 0.0) return 0;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = node[nodeID].edge[edgeID1];
                const GraphEdge &edge2 = node[nodeID].edge[edgeID2];
                GraphNode *node1 = &(node[std::abs(edge1.end) - 1]);
                GraphNode *node2 = &(node[std::abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, *node1, *node2)) continue;
                if (node1->isHomo || node2->isHomo) continue;

                if (calcNodeCoverage(node[nodeID]) > homoCoverageThreshold) continue;

                double coverage1 = this->calcNodeCoverage(*node1);
                double coverage2 = this->calcNodeCoverage(*node2);
                long id = std::abs(edge1.end);
                if (node1->length > node2->length) {
                    node1 = node2;
                    coverage1 = coverage2;
                    id = std::abs(edge2.end);
                }
                if (std::max(coverage1, coverage2) > heteroCoverageThreshold) continue;

                ++numDelete;
                node1->state |= SC_DEL;
                for (long edgeID = 0; edgeID < node1->numEdge; ++edgeID) {
                    ids.push_back(id);
                    ids.push_back(node1->edge[edgeID].end);
                }
                for (long contigID = 0; contigID < node1->numContig; ++contigID) {
                    contigPositionInScaffold[abs(node1->contig[contigID].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
                }
//std::cout << edge1.end << ", " << edge2.end << "; " << node[std::abs(edge1.end) - 1].length << ", " << node[std::abs(edge2.end) - 1].length << "; " << this->calcNodeCoverage(node[std::abs(edge1.end) - 1]) << ", " << this->calcNodeCoverage(node[std::abs(edge2.end) - 1]) << endl;
            }
        }
    }

    this->deleteEdges(ids);
    cerr << "NUM_SPLIT_LINK (not originate from heterozygosity)=" << numDelete << " (COVERAGE_THRESHOLD)" << endl;

    return numDelete;
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
long ScaffoldGraph::getSimilarOverlap(const long id1, const long id2)
{
    const double maxMissRate = MAX_OVERLAP_IDENTITY_DIFF;
    long tolerenceMiss = 0;
    long numMiss = 0;
    long exclusiveOR1, exclusiveOR2;
    long position1, position2;
    long code1, code2;
    long overlap = 0;

    if (id1 > 0) {
        if (id2 > 0) {
            exclusiveOR1 = exclusiveOR2 = 0;
            position1 = contig[id1 - 1].length - 1;
            position2 = 0;
            code1 = -1;
            code2 = 1;
        } else {
            exclusiveOR1 = 0;
            exclusiveOR2 = 1;
            position1 = contig[id1 - 1].length - 1;
            position2 = contig[std::abs(id2) - 1].length - 1;
            code1 = -1;
            code2 = -1;
        }
    } else {
        if (id2 > 0) {
            exclusiveOR1 = 1;
            exclusiveOR2 = 0;
            position1 = 0;
            position2 = 0;
            code1 = 1;
            code2 = 1;
        } else {
            exclusiveOR1 = exclusiveOR2 = 0;
            position1 = 0;
            position2 = contig[std::abs(id2) - 1].length - 1;
            code1 = 1;
            code2 = -1;
        }
    }

    long maxOverlap = std::min(contig[std::abs(id1) - 1].length, contig[std::abs(id2) - 1].length);
    for (long j = maxOverlap; j >= this->minOverlap; --j) {
        tolerenceMiss = j * maxMissRate + 0.5;
        numMiss = 0;
        for (long k = 0; k < j; ++k) {
            if ((exclusiveOR1 ^ contig[std::abs(id1) - 1].base[position1 + (k * code1)]) == (exclusiveOR2 ^ contig[std::abs(id2) - 1].base[position2 + ((j - k - 1) * code2)])) continue;

            ++numMiss;
            if (numMiss > tolerenceMiss)
                break;
        }
        if (numMiss <= tolerenceMiss) {
            overlap = j;
            break;
        }

    }

    return overlap;
}



//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::removeHeteroOverlap(void)
{
    long numDelete = 0;
    vector<long> ids;
    const double heteroCoverageThreshold = static_cast<long>(this->averageCoverage * MAX_HETERO_RATE + 0.5);

    if (this->bubbleThreshold == 0.0) return;


    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge *edge1 = &(node[nodeID].edge[edgeID1]);
                const GraphEdge *edge2 = &(node[nodeID].edge[edgeID2]);

                if (edge1->direction * edge2->direction < 0 || edge1->length < -(tolerence) || edge2->length < -(tolerence)) continue;
                GraphNode *node1 = &(node[abs(edge1->end) - 1]);
                if ((node1->state & SC_DEL) || edge1->length + node1->length <= edge2->length) continue;
                GraphNode *node2 = &(node[abs(edge2->end) - 1]);
                if ((node2->state & SC_DEL) || edge2->length + node2->length <= edge1->length) continue;
                if (node2->isHomo) continue;

                long edgeEnd1, edgeEnd2;
                if (edge1->isForward()) {
                    edgeEnd1 = edge1->end;
                    edgeEnd2 = edge2->end;
                } else {
                    edgeEnd1 = edge2->end;
                    edgeEnd2 = edge1->end;
                }

                if (abs(edge1->length + node1->length - edge2->length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2->length + node2->length - edge1->length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1)) continue;

                if (node1->length < node2->length) {
                    std::swap(node1, node2);
                    std::swap(edge1, edge2);
                }
                double coverage1 = this->calcNodeCoverage(*node1);
                double coverage2 = this->calcNodeCoverage(*node2);
                const GraphEdge *edge3 = NULL;
                long edgeID = 0;
                for (; edgeID < node1->numEdge; ++edgeID) {
                    if (abs(node1->edge[edgeID].end) == abs(edge2->end)) {
                        edge3 = &(node1->edge[edgeID]);
                        break;
                    }
                }

                if (edgeID == node1->numEdge) continue;

                if ((node2->state & SC_DEL) || coverage1 < coverage2 || edge3->length > -(tolerence) || edge3->length < -(node1->length) || coverage2 > heteroCoverageThreshold) continue;

                for (long contigID = 0; contigID < node2->numContig; ++contigID)
                    contigPositionInScaffold[abs(node2->contig[contigID].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
                node2->state |= SC_DEL;

                ++numDelete;
            }
        }
    }

    cerr << "NUM_REMOVED_OVERLAP_CONTIGS=" << numDelete << " (CONTAINED_HETERO)" << endl;

    long newNumNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        if (node[nodeID].state & SC_DEL) continue;

        fwrite(&(node[nodeID].numContig), sizeof(long), 1, scaffoldFP);
        for (long contigID = 0; contigID < node[nodeID].numContig; ++contigID)
            fwrite(&(node[nodeID].contig[contigID]), sizeof(ScaffoldPart), 1, scaffoldFP);
        newContigPoolSize += node[nodeID].numContig;
        ++newNumNode;
    }
    this->remake(newNumNode, newContigPoolSize, scaffoldFP);

    fclose(scaffoldFP);
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
void ScaffoldGraph::printScaffoldBubble(const string &outFilename)
{

    vector<char> seq;
    std::ofstream out(outFilename);
	if (bubbleFP != NULL) {
		rewind(bubbleFP);
		unsigned long long i, seqLength;
		double coverage;
		long numSeq = 0;

		while (fread(&seqLength, sizeof(unsigned long long), 1, bubbleFP)) {
			seq.resize(seqLength, 0);
			for (i = 0; i < seqLength; ++i)
				fread(&(seq[i]), sizeof(char), 1, bubbleFP);
			fread(&coverage, sizeof(double), 1, bubbleFP);

			++numSeq;
			out << ">seq" << numSeq << "_len" << seqLength << "_cov" << static_cast<unsigned short>(coverage + 0.5) << endl;
			for (i = 0; i < seqLength; ++i) {
				out.put(platanus::Bin2Char(seq[i]));
				if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
					out.put('\n');
			}
			if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
				out.put('\n');
		}
		fclose(bubbleFP);
	}

    out.close();
}



void ScaffoldGraph::insertSizeDistribution(vector<SeqLib>& library, vector<long>& distribution, const long numThread)
{
    long insertLength;
    platanus::Position forwardResult;
    platanus::Position reverseResult;

    for (long threadID = 0; threadID < numThread; ++threadID) {
        rewind(library[threadID].mappedFP);
        while (fread(&forwardResult, sizeof(platanus::Position), 1, library[threadID].mappedFP)) {
            fread(&reverseResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);

            unsigned long forwardIndex = abs(forwardResult.id) - 1;
            if (contigPositionInScaffold[forwardIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
            forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex][0].id : -(contigPositionInScaffold[forwardIndex][0].id); //edited by ouchi (contigPositionInScaffold)
            forwardResult.offset = contigPositionInScaffold[forwardIndex][0].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
            forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

            unsigned long reverseIndex = abs(reverseResult.id) - 1;
            if (contigPositionInScaffold[reverseIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
            reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseIndex][0].id : -(contigPositionInScaffold[reverseIndex][0].id); //edited by ouchi (contigPositionInScaffold)
            reverseResult.offset = contigPositionInScaffold[reverseIndex][0].id > 0 ? reverseResult.offset : contig[reverseIndex].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
            reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[reverseIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

            if (forwardResult.id != -reverseResult.id) continue;

            if (forwardResult.id > 0 && forwardResult.offset < reverseResult.offset)
                insertLength = static_cast<long>(reverseResult.offset - forwardResult.offset + 1);
            else if (reverseResult.id > 0 && reverseResult.offset < forwardResult.offset)
                insertLength = static_cast<long>(forwardResult.offset - reverseResult.offset + 1);
            else
                continue;

            if (insertLength >= static_cast<long>(distribution.size())) {
                distribution.resize(insertLength + 1, static_cast<long>(0));
            }
            ++distribution[insertLength];
        }
    }
}



void ScaffoldGraph::scaffoldLengthList(vector<long>& list)
{
    list.resize(numNode);
    for (long i = 0; i < numNode; ++i)
        list[i] = node[i].length;
}


void ScaffoldGraph::clearEdges()
{
    for (int i = 0; i < numNode; ++i) {
        node[i].edge.clear();
        node[i].numEdge= 0;
    }
}


void ScaffoldGraph::node2seq(const GraphNode &node, vector<char> &ret)
{
    platanus::SEQ *tmpContigPointer;
    long k;

    ret.clear();

	for (long j = 0; j < node.numContig; ++j) {
		if (j == 0) {
			k = 0;
		} else {
			k = node.contig[j - 1].end - node.contig[j].start;
		}
		for (; k < 0; ++k)
			ret.push_back(4);

		if (node.contig[j].id > 0) {
			tmpContigPointer = &(contig[node.contig[j].id - 1]);
			if (k < tmpContigPointer->length) {
				ret.insert(ret.end(), tmpContigPointer->base.begin() + k, tmpContigPointer->base.begin() + tmpContigPointer->length);
			}
		}
		else {
			tmpContigPointer = &(contig[-(node.contig[j].id) - 1]);
			for (; k < tmpContigPointer->length; ++k) {
				if (tmpContigPointer->base[tmpContigPointer->length - k - 1] >= 4)
					ret.push_back(4);
				else
					ret.push_back(0x3 ^ (tmpContigPointer->base[tmpContigPointer->length - k - 1]));
			}
		}
	}
}


void ScaffoldGraph::setContigName(const vector<string> &name)
{
	contigName = name;
}


void ScaffoldGraph::setContigNameIndex(const std::unordered_map<string, unsigned> nameIndex)
{
	contigNameIndex = nameIndex;
}


void ScaffoldGraph::dumpAllEdges(const string &outputFilename)
{
    std::ofstream out(outputFilename.c_str());

	out <<"#"
		<< "node1_ID" << "\t"
		<< "node1_len" << "\t"
		<< "node1_state" << "\t"
		<< "node1_component" << "\t"
		<< "edge_len" << "\t"
		<< "edge_num_links" << "\t"
		<< "node2_ID" << "\t"
		<< "node2_len" << "\t"
		<< "node2_state" << "\t"
		<< "node2_component"
	<< endl;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &nodeRef = this->node[nodeIndex];
        for (long edgeIndex = 0; edgeIndex < nodeRef.numEdge; ++edgeIndex) {
			GraphEdge &edgeRef = nodeRef.edge[edgeIndex];

			if (edgeRef.direction > 0) {
				out << nodeIndex + 1;
				out << "\t" << nodeRef.length;
				if (nodeRef.state & SC_REP)
					out << "\trep\t";
				else
					out << "\tnon_rep\t";

				for (long contigIndex = 0; contigIndex < nodeRef.numContig; ++contigIndex) {
					if (nodeRef.contig[contigIndex].id < 0)
						out << "-";
					out << contigName[id2Index(nodeRef.contig[contigIndex].id)];
					if (contigIndex + 1 < nodeRef.numContig)
						out << ";";
				}

				out << "\t" << edgeRef.length << "\t" << edgeRef.numLink << "\t" << edgeRef.end;

				GraphNode &endNodeRef = this->node[id2Index(edgeRef.end)];
				out << "\t" << endNodeRef.length;
				if (endNodeRef.state & SC_REP)
					out << "\trep\t";
				else
					out << "\tnon_rep\t";


				if (edgeRef.end > 0) {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[contigIndex].id < 0)
							out << "-";
						out << contigName[id2Index(endNodeRef.contig[contigIndex].id)];
						if (contigIndex + 1 < endNodeRef.numContig)
							out << ";";
					}
				}
				else {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id > 0)
							out << "-";
						out << contigName[id2Index(endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id)];
						if (contigIndex + 1 < endNodeRef.numContig)
							out << ";";
					}
				}
			}
			else {
				GraphNode &endNodeRef = this->node[id2Index(edgeRef.end)];
				out << edgeRef.end;
				out << "\t" << endNodeRef.length;
				if (endNodeRef.state & SC_REP)
					out << "\trep\t";
				else
					out << "\tnon_rep\t";

				if (edgeRef.end > 0) {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[contigIndex].id < 0)
							out << "-";
						out << contigName[id2Index(endNodeRef.contig[contigIndex].id)];
						if (contigIndex + 1 < endNodeRef.numContig)
							out << ";";
					}
				}
				else {
					for (long contigIndex = 0; contigIndex < endNodeRef.numContig; ++contigIndex) {
						if (endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id > 0)
							out << "-";
						out << contigName[id2Index(endNodeRef.contig[endNodeRef.numContig - contigIndex - 1].id)];
						if (contigIndex + 1 < endNodeRef.numContig)
							out << ";";
					}
				}

				out << "\t" << edgeRef.length << "\t" <<  edgeRef.numLink << "\t" << nodeIndex + 1;

				out << "\t" << nodeRef.length;
				if (nodeRef.state & SC_REP)
					out << "\trep\t";
				else
					out << "\tnon_rep\t";

				for (long contigIndex = 0; contigIndex < nodeRef.numContig; ++contigIndex) {
					if (nodeRef.contig[contigIndex].id < 0)
						out << "-";
					out << contigName[id2Index(nodeRef.contig[contigIndex].id)];
					if (contigIndex + 1 < nodeRef.numContig)
						out << ";";
				}
			}
			out << std::endl;
		}
	}

    out.close();
}


void ScaffoldGraph::outputFastg(const string &outputFilename)
{
    std::ofstream out(outputFilename.c_str());

//added by ouchi for debug
	for (long nodeID = 0; nodeID < numNode; ++nodeID) {
		cerr << nodeID + 1 << endl;
		const GraphNode &node = this->node[nodeID];
		cerr << "EDGE:" << node.numEdge << endl;
		for (long edgeID = 0; edgeID < node.numEdge; ++edgeID) {
			const GraphEdge &edge = node.edge[edgeID];
			cerr << "d:" << edge.direction << ", ID:" << edge.end << ", l:" << edge.length << ", n:" << edge.numLink << endl;
		}
	}
//


	vector<string> nodeName;
    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node = this->node[nodeID];
		std::ostringstream oss;
		oss << "NODE_" << nodeID + 1 << "_"
			<< "length_" << node.length << "_"
			<< "cov_" << calcNodeCoverage(node);

		nodeName.push_back(oss.str());
	}

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
		vector<char> nodeSeq;
		unsigned long i;
		unsigned currentNumEdge;

		unsigned numLeftEdge = 0;
		unsigned numRightEdge = 0;
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
			GraphEdge &edge = (node[nodeID].edge[edgeID]);
			if (edge.direction > 0)
				++numRightEdge;
			else
				++numLeftEdge;
		}

		currentNumEdge = 0;
		out << '>' << nodeName[nodeID];
		if (numRightEdge > 0)
			out << ':';
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
			GraphEdge &edge = (node[nodeID].edge[edgeID]);
			if (edge.direction < 0)
				continue;

			if (currentNumEdge > 0)
				out << ',';
			out << nodeName[abs(node1.edge[edgeID].end) - 1];
			if (node1.edge[edgeID].end < 0)
				out << '\'';
			++currentNumEdge;
		}
		out << ";\n";

		this->node2seq(node1, nodeSeq);
		for (i = 0; i < nodeSeq.size(); ++i) {
			out << platanus::Bin2Char(nodeSeq[i]);
            if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
		}
        if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');

		if (numLeftEdge == 0)
			continue;

		currentNumEdge = 0;
		out << '>' << nodeName[nodeID];
		out << '\'' << ':';
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
			GraphEdge &edge = (node[nodeID].edge[edgeID]);
			if (edge.direction > 0)
				continue;

			if (currentNumEdge > 0)
				out << ',';
			out << nodeName[abs(node1.edge[edgeID].end) - 1];
			if (node1.edge[edgeID].end > 0)
				out << '\'';
			++currentNumEdge;
		}
		out << ";\n";

		this->node2seq(node1, nodeSeq);
		for (i = 0; i < nodeSeq.size(); ++i) {
			if (nodeSeq[nodeSeq.size() - i - 1] < 4)
				out << platanus::Bin2Char(nodeSeq[nodeSeq.size() - i - 1] ^ 3);
			else
				out << platanus::Bin2Char(nodeSeq[nodeSeq.size() - i - 1]);
            if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
		}
        if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }

    out.close();
}


void ScaffoldGraph::loadResultSeq(const long minSeqLength, const unsigned long long readLength)
{
	resultSeq.clear();
	resultSeq.resize(node.size());

    vector<long> maxLeftOverlap(numContig, 0);
    vector<long> maxRightOverlap(numContig, 0);

    long minOverlap = this->minOverlap;
    this->minOverlap = 1;

    long maxNumContig = 0;
    long numScaffold = 0;
    long scaffoldLength = 0;
    for (unsigned i = 0; i < TABLE_DIVID; ++i) {
		auto overlapIterator = overlapTable[i].begin();
		auto overlapEnd = overlapTable[i].end();
		for (; overlapIterator != overlapEnd; ++overlapIterator) {
			if (overlapIterator->second.id1 > 0) {
				if (overlapIterator->second.length > maxRightOverlap[overlapIterator->second.id1 - 1])
					maxRightOverlap[overlapIterator->second.id1 - 1] = overlapIterator->second.length;
			} else if (overlapIterator->second.length > maxLeftOverlap[-(overlapIterator->second.id1) - 1])
				maxLeftOverlap[-(overlapIterator->second.id1) - 1] = overlapIterator->second.length;

			if (overlapIterator->second.id2 > 0) {
				if (overlapIterator->second.length > maxLeftOverlap[overlapIterator->second.id2 - 1])
					maxLeftOverlap[overlapIterator->second.id2 - 1] = overlapIterator->second.length;
			} else if (overlapIterator->second.length > maxRightOverlap[-(overlapIterator->second.id2) - 1])
				maxRightOverlap[-(overlapIterator->second.id2) - 1] = overlapIterator->second.length;
		}
    }
    for (long i = 0; i < numNode; ++i) {
        if (node[i].numContig > maxNumContig)
            maxNumContig = node[i].numContig;
    }
    vector<long> leftCut(maxNumContig, 0);
    vector<long> rightCut(maxNumContig, 0);
    vector<long> gap(maxNumContig, 0);

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        long contigID = 0;
        for (; contigID < node[nodeID].numContig; ++contigID) {
            if (contigPositionInScaffold[abs(node[nodeID].contig[contigID].id) - 1].size() == 1) //edited by ouchi (contigPositionInScaffold)
                break;
        }
        if (contigID == node[nodeID].numContig) continue;

        scaffoldLength = 0;
        if (node[nodeID].contig[0].id > 0)
            leftCut[0] = maxLeftOverlap[node[nodeID].contig[0].id - 1] / 2;
        else
            leftCut[0] = maxRightOverlap[-(node[nodeID].contig[0].id) - 1] / 2;

        for (contigID = 0; contigID < node[nodeID].numContig - 1; ++contigID) {
/*
            if (node[nodeID].contig[contigID].end > node[nodeID].contig[contigID + 1].start && this->getOverlap(node[nodeID].contig[contigID].id, node[nodeID].contig[contigID + 1].id) > minOverlap) {
                rightCut[contigID] = 0;
                leftCut[contigID + 1] = node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start;
                gap[contigID] = 0;
            } else if (node[nodeID].contig[contigID].end > node[nodeID].contig[contigID + 1].start) {
                rightCut[contigID] = leftCut[contigID + 1] = (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start) / 2;
                gap[contigID] = node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start;
            } else {
                rightCut[contigID] = leftCut[contigID + 1] = this->getSimilarOverlap(node[nodeID].contig[contigID].id, node[nodeID].contig[contigID + 1].id) / 2;
                gap[contigID] = node[nodeID].contig[contigID + 1].start - node[nodeID].contig[contigID].end + 2 * rightCut[contigID];
            }
*/
            if (node[nodeID].contig[contigID].end > node[nodeID].contig[contigID + 1].start) {
                rightCut[contigID] = leftCut[contigID + 1] = (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start) / 2;
                gap[contigID] = node[nodeID].contig[contigID].end - node[nodeID].contig[contigID + 1].start;
            } else {
                rightCut[contigID] = leftCut[contigID + 1] = this->getSimilarOverlap(node[nodeID].contig[contigID].id, node[nodeID].contig[contigID + 1].id) / 2;
                gap[contigID] = node[nodeID].contig[contigID + 1].start - node[nodeID].contig[contigID].end + 2 * rightCut[contigID];
            }
            scaffoldLength += (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID].start) - leftCut[contigID] - rightCut[contigID] + gap[contigID];
        }

        if (node[nodeID].contig[contigID].id > 0)
            rightCut[contigID] = maxRightOverlap[node[nodeID].contig[contigID].id - 1] / 2;
        else
            rightCut[contigID] = maxLeftOverlap[-(node[nodeID].contig[contigID].id) - 1] / 2;
        scaffoldLength += (node[nodeID].contig[contigID].end - node[nodeID].contig[contigID].start) - leftCut[contigID] - rightCut[contigID];
        gap[contigID] = 0;

        if (scaffoldLength < minSeqLength) continue;

		std::ostringstream nameSS;
        ++numScaffold;
        if (node[nodeID].numContig > 1)
            nameSS << "scaffold" << numScaffold << "_len" << scaffoldLength << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeID]) + 0.5) << "_read" << readLength;
        else
            nameSS << "scaffold" << numScaffold << "_len" << scaffoldLength << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeID]) + 0.5) << "_read" << readLength << "_single";

		resultSeq[nodeID].name = nameSS.str();

		resultSeq[nodeID].seq.clear();
        for (contigID = 0; contigID < node[nodeID].numContig; ++contigID) {
            if (node[nodeID].contig[contigID].id > 0) {
                for (long seqID = leftCut[contigID]; seqID < contig[node[nodeID].contig[contigID].id - 1].length - rightCut[contigID]; ++seqID)
                    resultSeq[nodeID].seq.push_back(contig[node[nodeID].contig[contigID].id - 1].base[seqID]);
            }
			else {
                for (long seqID = contig[-(node[nodeID].contig[contigID].id) - 1].length - leftCut[contigID] - 1; seqID >= (long)rightCut[contigID]; --seqID) {
                    if (contig[-(node[nodeID].contig[contigID].id) - 1].base[seqID] != 4)
                        resultSeq[nodeID].seq.push_back(0x3 ^ contig[-(node[nodeID].contig[contigID].id) - 1].base[seqID]);
                    else
                        resultSeq[nodeID].seq.push_back(4);
                }
            }

            for (long k = 0; k < gap[contigID]; ++k)
				resultSeq[nodeID].seq.push_back(4);
        }
    }

    this->minOverlap = minOverlap;
}

void ScaffoldGraph::printResultSeq(const string &outFilename)
{
    std::ofstream out(outFilename);

    for (long i = 0; i < numNode; ++i) {
		if (resultSeq[i].seq.size() <= 0 || resultSeq[i].redundantFlag)
			continue;

		out << '>' << resultSeq[i].name << endl;
		unsigned long j;
		for (j = 0; j < resultSeq[i].seq.size(); ++j) {
			out << platanus::Bin2Char(resultSeq[i].seq[j]);
			if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
				out.put('\n');
		}
		if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
			out.put('\n');
	}

	out.close();
}

void ScaffoldGraph::writeNodeSeq(const GraphNode &targetNode, FILE *fp)
{
	vector<char> nodeSeq;
	long size;

	size = targetNode.contig.size();
	fwrite(&size, sizeof(long), 1, bubbleFP);
	for (auto itr = targetNode.contig.begin(); itr != targetNode.contig.end(); ++itr)
		fwrite(&(itr->id), sizeof(long), 1, fp);

	node2seq(targetNode, nodeSeq);
	size = nodeSeq.size();
	fwrite(&size, sizeof(long), 1, bubbleFP);
	for (auto itr = nodeSeq.begin(); itr != nodeSeq.end(); ++itr)
		fwrite(&(*itr), sizeof(char), 1, fp);
}

long ScaffoldGraph::readNodeSeq(vector<char> &nodeSeq, vector<long> &contigID, FILE *fp)
{
	long size;

	fread(&size, sizeof(long), 1, bubbleFP);
	contigID.resize(size);
	for (auto itr = contigID.begin(); itr != contigID.end(); ++itr)
		fread(&(*itr), sizeof(long), 1, fp);

	long ret = fread(&size, sizeof(long), 1, bubbleFP);
	nodeSeq.resize(size);
	for (auto itr = nodeSeq.begin(); itr != nodeSeq.end(); ++itr)
		fread(&(*itr), sizeof(char), 1, fp);

	return ret;
}

void ScaffoldGraph::mapBubbleToResultSeq()
{
	rewind(bubbleFP);
    rewind(bubbleOpositeFP);

	vector<char> bubbleSeq;
	vector<char> bubbleOpositeSeq;
	vector<long> bubbleContigID;
	vector<long> bubbleOpositeContigID;

	long numMappedBubble = 0;

	while (readNodeSeq(bubbleSeq, bubbleContigID, bubbleFP)) {
		readNodeSeq(bubbleOpositeSeq, bubbleOpositeContigID, bubbleOpositeFP);

		long mappedNodeID = 0;
		long mappedNodeOffset = 0;
		for (auto itr = bubbleOpositeContigID.begin(); itr != bubbleOpositeContigID.end(); ++itr) {
			if (*itr != 0 && contigPositionInScaffold[abs(*itr) - 1].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				mappedNodeID = sign(*itr) * contigPositionInScaffold[abs(*itr) - 1][0].id; //edited by ouchi (contigPositionInScaffold)
				mappedNodeOffset = node[abs(mappedNodeID) - 1].contig[contigPositionInScaffold[abs(*itr) - 1][0].offset].start; //edited by ouchi (contigPositionInScaffold)
				break;
			}
		}
		if (mappedNodeID == 0)
			continue;

//		std::cout << mappedNodeID << ", " << mappedNodeOffset <<  endl;
		if (mappedNodeID < 0) {
			reverseComplement(bubbleSeq);
			reverseComplement(bubbleOpositeSeq);
		}

//		vector<char> &targetSeq = resultSeq[abs(mappedNodeID) - 1].seq;
		string &targetSeq = resultSeq[abs(mappedNodeID) - 1].seq;
		auto rangeStart = targetSeq.begin() + std::max(mappedNodeOffset - (long)(bubbleOpositeSeq.size()), 0L);
		auto rangeEnd = targetSeq.begin() + std::min(mappedNodeOffset + (long)(bubbleOpositeSeq.size()), (long)(targetSeq.size()));
		auto bubbleOpositeItr = std::search(rangeStart, rangeEnd, bubbleOpositeSeq.begin(), bubbleOpositeSeq.end());
		if (bubbleOpositeItr == rangeEnd) {
//			std::cout << "not found" <<  endl;
			continue;
		}
//		else
//			std::cout << bubbleOpositeItr - targetSeq.begin() <<  endl;

		long mapStart = bubbleOpositeItr - targetSeq.begin();
		long mapEnd = mapStart + bubbleOpositeSeq.size();

		std::pair<long, long> commonPartLength = trimCommonEdgePart(bubbleSeq.begin(), bubbleSeq.end(), bubbleOpositeSeq.begin(), bubbleOpositeSeq.end());
		mapStart += commonPartLength.first;
		mapEnd -= commonPartLength.second;

		++numMappedBubble;

//		std::cout << bubbleSeq.size() << ", " << bubbleOpositeSeq.size() <<  endl;
		std::cout << ">bubble" << numMappedBubble << "_in_consensus" << endl;
		for (auto itr = bubbleSeq.begin() + commonPartLength.first; itr != bubbleSeq.end() - commonPartLength.second; ++itr)
			std::cout << platanus::Bin2Char(*itr);
		std::cout << endl;

		std::cout << ">bubble" << numMappedBubble << "_removed" << endl;
		for (auto itr = targetSeq.begin() + mapStart; itr != targetSeq.begin() + mapEnd; ++itr)
			std::cout << platanus::Bin2Char(*itr);
		std::cout << endl;
	}
}

long ScaffoldGraph::getNumLinkFromIDPair(long leftNodeID, long rightNodeID)
{
	GraphNode &leftNode = this->node[abs(leftNodeID) - 1];
	char directionFromLeft = sign(leftNodeID);
	char rightStrand = directionFromLeft * sign(rightNodeID);

	for (long i = 0; i < leftNode.numEdge; ++i) {
		if (abs(leftNode.edge[i].end) == abs(rightNodeID) &&
			leftNode.edge[i].direction == directionFromLeft &&
			sign(leftNode.edge[i].end) == rightStrand) {

			return leftNode.edge[i].numLink;
		}
	}

	return 0;
}

void ScaffoldGraph::getLinkedNode(const long sourceNodeIndex, const char targetDirection, vector<NodeIDWithGap> &nodeIDBuffer)
{
	nodeIDBuffer.clear();
	GraphNode &node = this->node[sourceNodeIndex];
	for (long edgeIndex = 0; edgeIndex < node.numEdge; ++edgeIndex) {
		GraphEdge &edge = node.edge[edgeIndex];
		if (edge.direction == targetDirection && edge.numLink >= minLink)
			nodeIDBuffer.emplace_back(edge.end, edge.length);
	}
}

void ScaffoldGraph::getLinkedUniqueNode(const long sourceNodeIndex, const char targetDirection, vector<NodeIDWithGap> &nodeIDBuffer)
{
	long numLinkedNode = 0;

	nodeIDBuffer.clear();
	GraphNode &node = this->node[sourceNodeIndex];
	for (long edgeIndex = 0; edgeIndex < node.numEdge; ++edgeIndex) {
		GraphEdge &edge = node.edge[edgeIndex];
		if (edge.direction == targetDirection) {
			if (numLinkedNode > 0) {
				nodeIDBuffer.clear();
				return;
			}
			++numLinkedNode;
			if (edge.numLink >= minLink)
				nodeIDBuffer.emplace_back(edge.end, edge.length);
		}
	}
}

long ScaffoldGraph::remakeGraphAccordingToGappedPath(vector<GraphPathGapped> &pathBuffer)
{
	long numNewNode = 0;
	long newContigPoolSize = 0;
	long numJoin = 0;

    FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); ++pathItr) {
		//deleted by ouchi
		//if ((node[abs(pathItr->node.front().id) - 1].state & SC_INC) || (node[abs(pathItr->node.back().id) - 1].state & SC_INC))
		//	continue;
		//

		//added by ouchi
		auto nodeItr = pathItr->node.begin();
		for (; nodeItr != pathItr->node.end(); ++nodeItr) {
			if (node[abs(nodeItr->id) - 1].state & SC_INC)
				break;
		}
		if (nodeItr != pathItr->node.end())
			continue;
		//

		numJoin += pathItr->node.size() - 1; //edited by ouchi (pahtItr->node.size();)
		newContigPoolSize += writeAndMarkGappedNodes(pathItr->node, newContigFP);
		++numNewNode;
	}

	this->writeSingletonNode(numNewNode, newContigPoolSize, newContigFP);
    this->remake(numNewNode, newContigPoolSize, newContigFP);
    fclose(newContigFP);

	return numJoin;
}

long ScaffoldGraph::remakeGraphAccordingToGappedPathPair(vector<GraphPathGapped> &pathBuffer)
{
	long numNewNode = 0;
	long newContigPoolSize = 0;
	long numJoin = 0;
	vector<NodeIDWithGap> nodeForWriting;
	vector<NodeIDWithGap> nodeForWriting2;

    FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); pathItr += 2) {
		nodeForWriting.clear();
		auto nodeItr = pathItr->node.begin();
		for (; nodeItr != pathItr->node.end(); ++nodeItr) {
			if (node[abs(nodeItr->id) - 1].state & SC_INC)
				break;
			else
				nodeForWriting.push_back(*nodeItr);
		}
		if (nodeItr != pathItr->node.end())
			continue;

		nodeForWriting2.clear();
		auto pathItr2 = pathItr + 1;
		nodeItr = pathItr2->node.begin();
		for (; nodeItr != pathItr2->node.end(); ++nodeItr) {
			if (node[abs(nodeItr->id) - 1].state & SC_INC)
				break;
			else
				nodeForWriting2.push_back(*nodeItr);
		}
		if (nodeItr != pathItr2->node.end())
			continue;

		numJoin += (pathItr->node.size() + pathItr2->node.size());
		newContigPoolSize += writeAndMarkGappedNodes(nodeForWriting, newContigFP);
		newContigPoolSize += writeAndMarkGappedNodes(nodeForWriting2, newContigFP);
		numNewNode += 2;
	}

	this->writeSingletonNode(numNewNode, newContigPoolSize, newContigFP);
    this->remake(numNewNode, newContigPoolSize, newContigFP);
    fclose(newContigFP);

	return numJoin;
}


long ScaffoldGraph::writeAndMarkGappedNodes(const vector<NodeIDWithGap> &gappedPath, FILE *storeFP)
{
	long numContig = 0;

	for (auto it = gappedPath.begin(); it != gappedPath.end(); ++it)
		numContig += node[abs(it->id) - 1].numContig;

	fwrite(&numContig, sizeof(long), 1, storeFP);

	long start = 0;
	for (auto it = gappedPath.begin(); it != gappedPath.end(); ++it) {
		start += it->gap;

		long index;
		if (it->id > 0) {
			index = it->id - 1;
			node[index].state |= SC_INC;
			for (long i = 0; i < node[index].numContig; ++i) {
				ScaffoldPart scaffoldPartForWriting(node[index].contig[i].id, start + node[index].contig[i].start, start + node[index].contig[i].end);
				fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
			}
		}
		else {
			index = -(it->id) - 1;
			node[index].state |= SC_INC;
			for (long i = node[index].numContig - 1; i >= 0; --i) {
				ScaffoldPart scaffoldPartForWriting(-(node[index].contig[i].id), start + node[index].length - node[index].contig[i].end, start + node[index].length - node[index].contig[i].start);
				fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
			}
		}

		start += node[index].length;
	}

	return numContig;
}

void ScaffoldGraph::writeSingletonNode(long &numNewNode, long &newContigPoolSize, FILE *storeFP)
{
    for (long i = 0; i < numNode; ++i) {
        long numContig = node[i].numContig;
        if (!(node[i].state & SC_INC) && !(node[i].state & SC_DEL)) {
            fwrite(&(node[i].numContig), sizeof(long), 1, storeFP);
            for (long j = 0; j < numContig; ++j)
                fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, storeFP);
            newContigPoolSize += node[i].numContig;
            ++numNewNode;
        }
    }
}

void ScaffoldGraph::updateInsertLengthFP(vector<SeqLib>& lib, const long numThread)
{
    platanus::Position forwardResult;
    platanus::Position reverseResult;

	fseek(lib[0].insertLengthFP, 0, SEEK_END);

    for (long threadID = 0; threadID < numThread; ++threadID) {
        rewind(lib[threadID].mappedFP);
        while (fread(&forwardResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP)) {
            fread(&reverseResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);

            unsigned long forwardIndex = abs(forwardResult.id) - 1;
            if (contigPositionInScaffold[forwardIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)

            unsigned long reverseIndex = abs(reverseResult.id) - 1;
            if (contigPositionInScaffold[reverseIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)

			if (forwardIndex == reverseIndex) continue;

            forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex][0].id : -(contigPositionInScaffold[forwardIndex][0].id); //edited by ouchi (contigPositionInScaffold)
            forwardResult.offset = contigPositionInScaffold[forwardIndex][0].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
            forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

            reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseIndex][0].id : -(contigPositionInScaffold[reverseIndex][0].id); //edited by ouchi (contigPositionInScaffold)
            reverseResult.offset = contigPositionInScaffold[reverseIndex][0].id > 0 ? reverseResult.offset : contig[reverseIndex].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
            reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[reverseIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

            if (forwardResult.id != -reverseResult.id) continue;

            long insertLength = abs(forwardResult.offset - reverseResult.offset) ;
			fwrite(&insertLength, sizeof(long), 1, lib[0].insertLengthFP);
        }
    }
}

long ScaffoldGraph::joinUnambiguousNodePairGapped(const long numThread)
{
	vector<vector<GraphPathGapped> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<NodeIDWithGap> nodeBuffer;
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & SC_DEL)
				continue;

			for (long i = 0; i < 2; ++i) {
				char initialDirection = 2*i - 1;
//				getLinkedNode(initialNodeIndex, initialDirection, nodeBuffer);
				getLinkedUniqueNode(initialNodeIndex, initialDirection, nodeBuffer);
				if (nodeBuffer.size() != 1)
					continue;

				long altNodeID = nodeBuffer[0].id;

				getLinkedNode(abs(altNodeID) - 1, -(initialDirection * sign(altNodeID)), nodeBuffer);
				if (nodeBuffer.size() != 1)
					continue;

				pathBuffer[t].emplace_back(initialNodeIndex, 2);
				if (i == 0) {
					pathBuffer[t].back().node[0].id = altNodeID;
					pathBuffer[t].back().node[1].id = initialNodeIndex + 1;
				}
				else {
					pathBuffer[t].back().node[0].id = initialNodeIndex + 1;
					pathBuffer[t].back().node[1].id = altNodeID;
				}
				pathBuffer[t].back().node[0].gap = 0;
				pathBuffer[t].back().node[1].gap = nodeBuffer[0].gap;

				break;
			}
		}
	}

	vector<GraphPathGapped> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathGappedSelfIDLess());

	long numJoin = remakeGraphAccordingToGappedPath(mergedPathBuffer);

	return numJoin;
}

void ScaffoldGraph::joinUnambiguousNodePairGappedIterative(const long numThread)
{
//	const double MAX_HETERO_COVERAGE_FACTOR = 1.5;

	long total = 0;
	long num;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "joining unambiguous pair of nodes in a scaffold graph.." << endl;
	do {
		setMinLink(1);
		makeGraph(numThread);
		setMinLink(currentMinLink);
//		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
//		deleteNonOverlapHomoEdge();
		num = joinUnambiguousNodePairGapped(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);
	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}


void ScaffoldGraph::deleteEdgeFromShortNode(const long lengthThreshold)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from short nodes (LENGTH < " <<  lengthThreshold << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (node1.length < lengthThreshold) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

void ScaffoldGraph::deleteEdgeFromShortAndAbnormalCoverageNode(const long minLength, const double minCoverage, const double maxCoverage)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from short and abnormal-coverage nodes (LENGTH < " <<  minLength << ", COVERAGE < " << minCoverage <<  ", " << maxCoverage << " < COVERAGE) ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (node1.length < minLength) {
				double coverage = this->calcNodeCoverage(node1);
				if (coverage < minCoverage || coverage > maxCoverage) {
					ids.push_back(nodeID + 1);
					ids.push_back(node1.edge[edgeID].end);
					++numDelete;
				}
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

void ScaffoldGraph::setContigUnlinkFlags(const string fileName)
{
	std::pair<string, string> contigNamePair;
	string line;

    std::ifstream ifs(fileName);

    while (ifs && getline(ifs, line)) {
		std::stringstream ss;
		ss << line;
		ss >> contigNamePair.first >> contigNamePair.second;

		auto itr1 = this->contigNameIndex.find(contigNamePair.first);
		auto itr2 = this->contigNameIndex.find(contigNamePair.second);
		if (itr1 == this->contigNameIndex.end() || itr2 == this->contigNameIndex.end())
			continue;

		this->contigUnlinkSet.insert(itr1->second < itr2->second ? std::make_pair(itr1->second, itr2->second) : std::make_pair(itr2->second, itr1->second));
	}

    ifs.close();
}
