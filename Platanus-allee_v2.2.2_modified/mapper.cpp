/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus-allee.

Platanus-allee is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus-allee is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus-allee; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "mapper.h"
#include <omp.h>

using std::vector;
using std::cerr;
using std::endl;


const double Mapper::MIN_IDENTITY_FOR_GAP_CLOSE = 0.98;
const double Mapper::MIN_IDENTITY_FOR_SCAFFOLD = 0.95;
const double Mapper::MIN_IDENTITY_TO_CHECK_MAPPING = 0.95;
const int Mapper::DEFAULT_MAPKEY_LIMIT = 256;
const double HeteroMapper::MIN_IDENTITY_FOR_SCAFFOLD = 0.95;
const long Mapper::DISPLAY_NUM_READ_UNIT = 10000000;


Mapper::Mapper(): keyLength(0), mapKeyUpperLimit(DEFAULT_MAPKEY_LIMIT), mask(0), seedLength(0), indexLength(0), seqPoolSize(0), numSeq(0), seq(), table(), positionPool(nullptr)
{
}

Mapper::Mapper(const long seed, const long key): keyLength(key), mapKeyUpperLimit(DEFAULT_MAPKEY_LIMIT), mask(0), seedLength(seed), indexLength(0), seqPoolSize(0), numSeq(0), seq(), table(), positionPool(nullptr)
{
	keyLength = key > 32 ? 32 : key;
	mask = keyLength < 32 ? ~(~0ull << (2 * keyLength)) : ~0ull;
}


//////////////////////////////////////////////////////////////////////////////////////
// initilize mapper (not constructor)
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::setContig(platanus::Contig &contig)
{
    numSeq = contig.numSeq;
    seq = std::move(contig.seq);
    for (auto it = seq.begin(), end = seq.end(); it != end; ++it)
        seqPoolSize += it->base.size();

    nameIndex = contig.nameIndex;
}





//////////////////////////////////////////////////////////////////////////////////////
// insert contigs into kmer hash table
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::insertHashKey(void)
{
    unsigned long long seqLength = 0;
    unsigned long long hashKey = 0;
    unsigned long long start = 0;
    unsigned long long total = 0;
    bool keyInit = false;
    unsigned j = 0;
    FILE *tmpFP = platanus::makeTemporaryFile();
    // make this struct
    platanus::Position position;
    table.clear();

    for (int i = 0; i < numSeq; ++i) {
        hashKey = 0;
        seqLength = seq[i].length;
        if (static_cast<long>(seqLength) < keyLength) continue;
        start = 0;
        keyInit = true;
        while (start < seqLength - keyLength + 1) {
            // initialize hash key
            if (keyInit) {
                hashKey = 0;
                for (j = 0; j < static_cast<unsigned>(keyLength - 1); ++j) {
                    // if sequence base contain "N", don't make hash key
                    if (seq[i].base[start + j] == 4) break;
                    hashKey = (hashKey << 2) | static_cast<unsigned long long>(seq[i].base[start + j]);
                }
                if (j == static_cast<unsigned>(keyLength - 1)) {
                    keyInit = false;
                } else {
                    start += j + 1;
                    continue;
                }
            }
            if (seq[i].base[start + keyLength - 1] == 4) {
                start += keyLength;
                keyInit = true;
                continue;
            }

            hashKey = ((hashKey << 2) & mask) | static_cast<unsigned long long>(seq[i].base[start + keyLength - 1]);
            position.id = i + 1;
            position.offset = start;
		//this line add key to table,not but .position
            ++(table[hashKey].num);
            fwrite(&hashKey, sizeof(unsigned long long), 1, tmpFP);
            fwrite(&position, sizeof(platanus::Position), 1, tmpFP);
            ++total;
            ++start;
        }
    }
    platanus::Position *pos = reservePositionPool(total);

    fillPositionPool(tmpFP, pos);
}



//////////////////////////////////////////////////////////////////////////////////////
// new positionPool using smart pointer
//////////////////////////////////////////////////////////////////////////////////////
platanus::Position *Mapper::reservePositionPool(const unsigned long long total)
{
    positionPool.reset(new platanus::Position[total]);
    return positionPool.get();
}



//////////////////////////////////////////////////////////////////////////////////////
// fill positionPool and complete map table
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::fillPositionPool(FILE *tmpFP, platanus::Position *pos)
{
    unsigned long long hashKey;
    platanus::Position position;
        rewind(tmpFP);
        while (fread(&hashKey, sizeof(unsigned long long), 1, tmpFP)) {
            fread(&position, sizeof(platanus::Position), 1, tmpFP);
            auto it = table.find(hashKey);
            if (it->second.position == NULL) {
                it->second.position = pos;
                pos += it->second.num;
                it->second.num = 1;
                *(it->second.position) = position;
            } else {
                *(it->second.position + it->second.num) = position;
                ++(it->second.num);
            }
        }
        fclose(tmpFP);
}


//////////////////////////////////////////////////////////////////////////////////////
// make kmer table and calculate max occurrence value
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::makeKmerTable(void)
{
    cerr << "K=" << keyLength << ", making hash table..." << endl;
    this->insertHashKey();
}



//////////////////////////////////////////////////////////////////////////////////////
// mapping kmer (however, only k <= 32)
//////////////////////////////////////////////////////////////////////////////////////
long Mapper::mapKey(const Kmer31 &kmer, vector<platanus::Position> &positionBuffer) const
{
    positionBuffer.clear();
    auto tableIt = table.find(kmer.forward);
    if (tableIt != table.end()) {
	  positionBuffer.resize(tableIt->second.num);
	  std::copy(tableIt->second.position, tableIt->second.position + tableIt->second.num, positionBuffer.begin());
    }
    unsigned forwardSize = positionBuffer.size();
    tableIt = table.find(kmer.reverse);
    if (tableIt != table.end()) {
	  positionBuffer.insert(positionBuffer.end(), tableIt->second.position, tableIt->second.position + tableIt->second.num);
    }

    for (unsigned i = forwardSize; i < positionBuffer.size(); ++i)
	  positionBuffer[i].id *= -1;
    return positionBuffer.size();

}


long Mapper::mapKeyWithLimit(const Kmer31 &kmer, vector<platanus::Position> &positionBuffer) const
{
    positionBuffer.clear();
    auto tableIt = table.find(kmer.forward);
    if (tableIt != table.end() && tableIt->second.num <= this->mapKeyUpperLimit) {
	  positionBuffer.resize(tableIt->second.num);
	  std::copy(tableIt->second.position, tableIt->second.position + tableIt->second.num, positionBuffer.begin());
    }
    unsigned forwardSize = positionBuffer.size();
    tableIt = table.find(kmer.reverse);
    if (tableIt != table.end() && tableIt->second.num <= this->mapKeyUpperLimit) {
	  positionBuffer.insert(positionBuffer.end(), tableIt->second.position, tableIt->second.position + tableIt->second.num);
    }

    for (unsigned i = forwardSize; i < positionBuffer.size(); ++i)
	  positionBuffer[i].id *= -1;
    return positionBuffer.size();

}


//////////////////////////////////////////////////////////////////////////////////////
// Merge bubble and contig by mapping
//////////////////////////////////////////////////////////////////////////////////////
void HeteroMapper::mergeBubble(const long numThread)
{
    long threadID = 0;
    long seqID = 0;
    long bufferSize = 0;
    long maxLeftLength = 0;
    long maxRightLength = 0;
    vector<platanus::Position> positionBuffer;
    platanus::Position leftResult;
    platanus::Position rightResult;
    if (bubbleMap.getNumSeq() == 0) return;

    cerr << "mapping bubbles on contigs..." << endl;
    if (bubblePosition.size() != 0)
        bubblePosition.clear();

    bubblePosition.resize(bubbleMap.getNumSeq());
    omp_set_num_threads(numThread);

    #    pragma omp parallel for  schedule(static, 1) \
        private(seqID, positionBuffer, bufferSize, maxLeftLength, maxRightLength) \
        firstprivate(leftResult, rightResult)

    for (threadID = 0; threadID < numThread; ++threadID) {
        long bubbleNumSeq = bubbleMap.getNumSeq();
        for (seqID = threadID; seqID < bubbleNumSeq; seqID += numThread) {
            Kmer31 leftKmer(keyLength);
            Kmer31 rightKmer(keyLength);
            const platanus::SEQ &bubbleSeq = bubbleMap.getSeq(seqID);
            if (bubbleSeq.length < 2 * keyLength) continue;

            if (leftKmer.setKmer(bubbleSeq, 0, keyLength, seedLength) || rightKmer.setKmer(bubbleSeq, bubbleSeq.length - keyLength, keyLength, seedLength)) continue;
            // mapping bubble left seed in contig
            bufferSize = contigMap.mapKey(leftKmer, positionBuffer);
            maxLeftLength = 0;
            for (long i = 0; i < bufferSize; ++i) {
                long j;
                // position.id > 0 means seed is mapped forward contig
                if (positionBuffer[i].id > 0) {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + j >= contigMap.getSeqLength(positionBuffer[i].id - 1)
                            || contigMap.getBase(positionBuffer[i].id - 1, positionBuffer[i].offset + j) != bubbleSeq.base[j])
                            break;
                    }
                }
                // position.id < 0 means seed is mapped reverse contig
                else {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + keyLength - j <= 0
                            || contigMap.getBase(-(positionBuffer[i].id) - 1, positionBuffer[i].offset + keyLength - j - 1) != (0x3^bubbleSeq.base[j]))
                            break;
                    }
                    positionBuffer[i].offset += keyLength - 1;
                }
                if (j > maxLeftLength) {
                    maxLeftLength = j;
                    leftResult =  positionBuffer[i];
                } else if (j == maxLeftLength)
                    leftResult.id = 0;
            }

            if (maxLeftLength < seedLength) continue;

            // mapping bubble right seed in contig
            bufferSize = contigMap.mapKey(rightKmer, positionBuffer);
            maxRightLength = 0;
            for (long i = 0; i < bufferSize; ++i) {
                long j;
                // position.id > 0 means seed is mapped forward contig
                if (positionBuffer[i].id > 0) {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + keyLength - j <= 0
                            || contigMap.getBase(positionBuffer[i].id - 1, positionBuffer[i].offset + keyLength - j - 1) != bubbleSeq.base[bubbleSeq.length - j - 1])
                            break;
                    }
                    positionBuffer[i].offset += keyLength - 1;
                }
                // position.id < 0 means seed is mapped reverse contig
                else {
                    for (j = keyLength; j < bubbleSeq.length; ++j) {
                        if (positionBuffer[i].offset + j >= contigMap.getSeqLength(-(positionBuffer[i].id) - 1)
                            || contigMap.getBase(-(positionBuffer[i].id) - 1, positionBuffer[i].offset + j) != (0x3^bubbleSeq.base[bubbleSeq.length - j - 1]))
                            break;
                    }
                }
                if (j > maxRightLength) {
                    maxRightLength = j;
                    rightResult =  positionBuffer[i];
                } else if (j == maxRightLength)
                    rightResult.id = 0;
            }

            if (maxRightLength < seedLength || rightResult.id != leftResult.id || leftResult.id == 0 || rightResult.id == 0) continue;

            bubblePosition[seqID].id = leftResult.id;
            bubblePosition[seqID].offset = (leftResult.offset + rightResult.offset) / 2;
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// decide the contig where the variable read is mapped
// this version, I only translate C language.
// it can be more acceleratable using kmer mapping over k > 32
//////////////////////////////////////////////////////////////////////////////////////
platanus::Position Mapper::mapRead(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer, const long wordLength)
{
    long numResult = 0;
    // vector<platanus::Position> result(platanus::ConstParam::MAX_READ_LEN);
    vector<platanus::Position> result(read.length/wordLength + 2);
    Kmer31 kmer;
    long bufferSize;
    long firstMax, secondMax;
    int i, j, k = 0;
    for (i = read.length - wordLength; i > -(wordLength); i -= wordLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j -1]);
        }

        if (j != keyLength) continue;

        // first, mapping kmer <= 32
        bufferSize = mapKeyWithLimit(kmer, positionBuffer);
        result[numResult].id = 0;
        // extend mapping kmer >= 32
        for (j = 0; j < bufferSize; ++j) {
            // position.id > 0 means seed is mapped forward contig
            if (positionBuffer[j].id > 0) {
                if (positionBuffer[j].offset > seq[positionBuffer[j].id - 1].length - wordLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (seq[positionBuffer[j].id - 1].base[positionBuffer[j].offset + k] != read.base[i + k])
                        break;
                }
                positionBuffer[j].offset -= i;
            }
            // position.id < 0 means seed is mapped reverse contig
            else {
                if (positionBuffer[j].offset < wordLength - keyLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (seq[-(positionBuffer[j].id) - 1].base[positionBuffer[j].offset + keyLength - k - 1] != (3^read.base[i + k]))
                        break;
                }
                positionBuffer[j].offset += i + keyLength - 1;
            }

            if (k == wordLength) {
                if (result[numResult].id == 0)
                    result[numResult] = positionBuffer[j];
                else
                    break;
            }
        }

        // check whther exist exact match kmer read or not
	  //seed not multi hit and not unmapped
        if (j == bufferSize && result[numResult].id != 0)
            ++numResult;
		//num of valid seed
    }
    if (numResult == 0) {
        result[0].id = 0;
        return result[0];
    }
    result.resize(numResult);
    sort(result.begin(), result.end());

    // select the highest overlapped kmer one
    result.resize(numResult + 1);
    i = firstMax = secondMax = result[numResult].id = 0;
    while (result[i].id != 0) {
        j = i;
        while (result[i].id == result[i+1].id && result[i].offset == result[i+1].offset)
            ++i;
        ++i;
        if (i - j >= firstMax) {
            secondMax = firstMax;
            firstMax = i - j;
            k = j;
        }
    }

    // if cannot decide top one, assign unmapped
    if (firstMax == secondMax) {
        result[0].id = 0;
        return result[0];
    }

    return result[k];
}

void Mapper::mapReadMultiReports(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer, const long wordLength, const long minContigLength)
{
	positionBuffer.clear();
    long numResult = 0;
    // vector<platanus::Position> result(platanus::ConstParam::MAX_READ_LEN);
    vector<platanus::Position> result(read.length/wordLength + 2);
    Kmer31 kmer;
    long bufferSize;
    int i, j, k = 0;

    for (i = read.length - wordLength; i > -(wordLength); i -= wordLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j -1]);
        }

        if (j != keyLength) continue;

        // first, mapping kmer <= 32
        bufferSize = mapKeyWithLimit(kmer, positionBuffer);
        result[numResult].id = 0;
        // extend mapping kmer >= 32
        for (j = 0; j < bufferSize; ++j) {
			if (seq[abs(positionBuffer[j].id) - 1].length < minContigLength)
				continue;

            // position.id > 0 means seed is mapped forward contig
            if (positionBuffer[j].id > 0) {
                if (positionBuffer[j].offset > seq[positionBuffer[j].id - 1].length - wordLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (seq[positionBuffer[j].id - 1].base[positionBuffer[j].offset + k] != read.base[i + k])
                        break;
                }
                positionBuffer[j].offset -= i;
            }
            // position.id < 0 means seed is mapped reverse contig
            else {
                if (positionBuffer[j].offset < wordLength - keyLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (seq[-(positionBuffer[j].id) - 1].base[positionBuffer[j].offset + keyLength - k - 1] != (3^read.base[i + k]))
                        break;
                }
                positionBuffer[j].offset += i + keyLength - 1;
            }

            if (k == wordLength) {
                if (result[numResult].id == 0)
                    result[numResult] = positionBuffer[j];
                else
                    break;
            }
        }

        if (j == bufferSize && result[numResult].id != 0)
            ++numResult;
    }
    if (numResult == 0) {
		positionBuffer.clear();
        return;
    }
    result.resize(numResult);
    sort(result.begin(), result.end());

    result.resize(numResult + 1);
    result[numResult].id = 0;

	positionBuffer.clear();
    i = 0;
    while (result[i].id != 0) {
        j = i;
		long sum = result[i].offset;
        while (result[i].id == result[i+1].id) {
            ++i;
			sum += result[i].offset;
		}
        ++i;
		positionBuffer.resize(positionBuffer.size() + 1);
		positionBuffer.back().id = result[j].id;
		positionBuffer.back().offset = std::round(static_cast<double>(sum) / (i - j));
    }
}

platanus::Position Mapper::mapReadMultiSeedFiltered(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer)
{
	platanus::Position position;
	for (auto seedItr = this->multiSeedLength.begin(); seedItr != this->multiSeedLength.end(); ++seedItr) {
		position = this->mapRead(read, positionBuffer, *seedItr);
//		if (!(this->judgeMappedPositionUngapAlignment(read, position, MIN_IDENTITY_TO_CHECK_MAPPING)))
//			position.id = 0;

		if (position.id != 0)
			break;
	}

	return position;
}

void Mapper::mapReadMultiReportsMultiSeed(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer, const long minContigLength)
{
	for (auto seedItr = this->multiSeedLength.begin(); seedItr != this->multiSeedLength.end(); ++seedItr) {
		this->mapReadMultiReports(read, positionBuffer, *seedItr, minContigLength);
		if (positionBuffer.size() > 0)
			break;
	}
}

void Mapper::mapReadMultiReportsMultiSeedFiltered(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer, const long minContigLength)
{
	for (auto seedItr = this->multiSeedLength.begin(); seedItr != this->multiSeedLength.end(); ++seedItr) {
		this->mapReadMultiReports(read, positionBuffer, *seedItr, minContigLength);
//		this->filterMappedPositionUngapAlignment(read, positionBuffer, MIN_IDENTITY_TO_CHECK_MAPPING);
		if (positionBuffer.size() > 0)
			break;
	}
}

platanus::Position Mapper::mapReadUngapAlignment(const platanus::SEQ &read, const double minIdentity, vector<platanus::Position> &positionBuffer, double *returnIdentity)
{
	const long MATCH_SCORE= 1;
	const long MISMATCH_SCORE = -50;

    long i;
    long j;
    long k = 0;
    long numHit;
    long numMis;
    long misThreshold;
    long alignmentStart;
    long alignmentEnd;
    long bufferSize;
	long score;
    long maxScore;
	platanus::Position maxScorePosition;
    Kmer31 kmer;

	numHit = 0;
	maxScore = LONG_MIN;
	maxScorePosition.id = 0;

    for (i = read.length - seedLength; i > -(seedLength); i -= seedLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j -1]);
        }

        if (j != keyLength) continue;

        bufferSize = mapKey(kmer, positionBuffer);
        for (j = 0; j < bufferSize; ++j) {
			numMis = 0;
			platanus::Position &seedPosition = positionBuffer[j];
            if (seedPosition.id > 0) {
				platanus::SEQ &target = seq[seedPosition.id - 1];

				if (seedPosition.id == maxScorePosition.id && seedPosition.offset - i == maxScorePosition.offset)
					continue;

                if (seedPosition.offset > target.length - seedLength)
                    continue;

                for (k = keyLength; k < seedLength; ++k) {
                    if (target.base[seedPosition.offset + k] != read.base[i + k])
                        break;
                }
				if (k != seedLength)
					continue;

				alignmentStart = -std::min(static_cast<long>(seedPosition.offset), i);
				alignmentEnd = std::min(read.length - i , target.length - seedPosition.offset);
				misThreshold = (int64_t)((1.0 - minIdentity) * (alignmentEnd - alignmentStart));

				for (k = alignmentStart; k < 0; ++k) {
					if (target.base[seedPosition.offset + k] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + k] != read.base[i + k]) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != 0)
					continue;

				for (k = seedLength; k < alignmentEnd; ++k) {
					if (target.base[seedPosition.offset + k] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + k] != read.base[i + k]) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != alignmentEnd)
					continue;

				seedPosition.offset -= i;
            }
            else {
				platanus::SEQ &target = seq[-(seedPosition.id) - 1];

				if (seedPosition.id == maxScorePosition.id && seedPosition.offset + i + keyLength - 1 == maxScorePosition.offset)
					continue;

				if (seedPosition.offset < seedLength - keyLength)
                    continue;

                for (k = keyLength; k < seedLength; ++k) {
                    if (target.base[seedPosition.offset + keyLength - k - 1] != (3^read.base[i + k]))
                        break;
                }
				if (k != seedLength)
					continue;

				alignmentStart = -std::min(i, target.length - seedPosition.offset - keyLength);
				alignmentEnd = std::min(read.length - i , static_cast<long>(seedPosition.offset));
				misThreshold = static_cast<long>((1.0 - minIdentity) * (alignmentEnd - alignmentStart));

				for (k = alignmentStart; k < 0; ++k) {
					if (target.base[seedPosition.offset + keyLength - k - 1] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + keyLength - k - 1] != (3^read.base[i + k])) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != 0)
					continue;;

				for (k = seedLength; k < alignmentEnd; ++k) {
					if (target.base[seedPosition.offset + keyLength - k - 1] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + keyLength - k - 1] != (3^read.base[i + k])) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != alignmentEnd)
					continue;

				seedPosition.offset += i + keyLength - 1;
            }

			score = MATCH_SCORE*(alignmentEnd - alignmentStart - numMis) +  MISMATCH_SCORE*numMis;

			if (score > maxScore) {
				numHit = 1;
				maxScore = score;
				maxScorePosition = seedPosition;
				if (returnIdentity != NULL)
					*returnIdentity = static_cast<double>((alignmentEnd - alignmentStart - numMis)) / (alignmentEnd - alignmentStart);
			}
			else if (score == maxScore) {
				++numHit;
			}
		}
    }

	if (numHit != 1) {
		maxScorePosition.id = 0;
		positionBuffer.clear();
	}
	else
		positionBuffer.assign(1, maxScorePosition);

	return maxScorePosition;
}

void Mapper::mapReadUngapAlignmentMulti(const platanus::SEQ &read, const double minIdentity, vector<platanus::Position> &positionBuffer, vector<platanus::Position> &maxScorePositionBuffer)
{
	const long MATCH_SCORE= 1;
	const long MISMATCH_SCORE = -50;

    long i;
    long j;
    long k = 0;
    long numMis;
    long misThreshold;
    long alignmentStart;
    long alignmentEnd;
    long bufferSize;
	long score;
    long maxScore;
    Kmer31 kmer;

	maxScore = LONG_MIN;

	maxScorePositionBuffer.clear();

    for (i = read.length - seedLength; i > -(seedLength); i -= seedLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j -1]);
        }

        if (j != keyLength) continue;

        bufferSize = mapKey(kmer, positionBuffer);
        for (j = 0; j < bufferSize; ++j) {
			numMis = 0;
			platanus::Position &seedPosition = positionBuffer[j];
            if (seedPosition.id > 0) {
				platanus::SEQ &target = seq[seedPosition.id - 1];

				auto itr = maxScorePositionBuffer.begin();
				for (; itr != maxScorePositionBuffer.end(); ++itr) {
					if (seedPosition.id == itr->id && seedPosition.offset - i == itr->offset)
						break;
				}
				if (itr != maxScorePositionBuffer.end())
					continue;

                if (seedPosition.offset > target.length - seedLength)
                    continue;

                for (k = keyLength; k < seedLength; ++k) {
                    if (target.base[seedPosition.offset + k] != read.base[i + k])
                        break;
                }
				if (k != seedLength)
					continue;

				alignmentStart = -std::min(static_cast<long>(seedPosition.offset), i);
				alignmentEnd = std::min(read.length - i , target.length - seedPosition.offset);
				misThreshold = (int64_t)((1.0 - minIdentity) * (alignmentEnd - alignmentStart));

				for (k = alignmentStart; k < 0; ++k) {
					if (target.base[seedPosition.offset + k] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + k] != read.base[i + k]) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != 0)
					continue;

				for (k = seedLength; k < alignmentEnd; ++k) {
					if (target.base[seedPosition.offset + k] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + k] != read.base[i + k]) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != alignmentEnd)
					continue;

				seedPosition.offset -= i;
            }
            else {
				platanus::SEQ &target = seq[-(seedPosition.id) - 1];

				auto itr = maxScorePositionBuffer.begin();
				for (; itr != maxScorePositionBuffer.end(); ++itr) {
					if (seedPosition.id == itr->id && seedPosition.offset - i == itr->offset)
					if (seedPosition.id == itr->id && seedPosition.offset + i + keyLength - 1 == itr->offset)
						break;
				}
				if (itr != maxScorePositionBuffer.end())
					continue;

				if (seedPosition.offset < seedLength - keyLength)
                    continue;

                for (k = keyLength; k < seedLength; ++k) {
                    if (target.base[seedPosition.offset + keyLength - k - 1] != (3^read.base[i + k]))
                        break;
                }
				if (k != seedLength)
					continue;

				alignmentStart = -std::min(i, target.length - seedPosition.offset - keyLength);
				alignmentEnd = std::min(read.length - i , static_cast<long>(seedPosition.offset));
				misThreshold = static_cast<long>((1.0 - minIdentity) * (alignmentEnd - alignmentStart));

				for (k = alignmentStart; k < 0; ++k) {
					if (target.base[seedPosition.offset + keyLength - k - 1] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + keyLength - k - 1] != (3^read.base[i + k])) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != 0)
					continue;;

				for (k = seedLength; k < alignmentEnd; ++k) {
					if (target.base[seedPosition.offset + keyLength - k - 1] != 4 && read.base[i + k] != 4 && target.base[seedPosition.offset + keyLength - k - 1] != (3^read.base[i + k])) {
						++numMis;
						if (numMis > misThreshold)
							break;
					}
				}
				if (k != alignmentEnd)
					continue;

				seedPosition.offset += i + keyLength - 1;
            }

			score = MATCH_SCORE*(alignmentEnd - alignmentStart - numMis) +  MISMATCH_SCORE*numMis;

			if (score > maxScore) {
				maxScore = score;
				maxScorePositionBuffer.assign(1, seedPosition);
			}
			else if (score == maxScore) {
				maxScorePositionBuffer.push_back(seedPosition);
			}
		}
    }

	return;
}

//////////////////////////////////////////////////////////////////////////////////////
// map each read and count the number of read mapped same contig and different one
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapPairMT(vector<SeqLib> &library, const long minInsertion, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    long numMappedDiffContig = 0;
    long sum = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;

    cerr << "mapping reads..." << endl;

    omp_set_num_threads(numThread);


    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: sum, totalPair, numMappedSameContig, numMappedDiffContig)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;


            if (forward.length < seedLength || reverse.length < seedLength) continue;

            // find the contig where this read is mapped
//			double identity;
//            forwardPosition = this->mapReadUngapAlignment(forward, MIN_IDENTITY_FOR_SCAFFOLD, positionBuffer, identity);
//            if (forwardPosition.id == 0) {
                forwardPosition = this->mapReadMultiSeedFiltered(forward, positionBuffer);
                if (forwardPosition.id == 0) continue;
//            }
//            reversePosition = this->mapReadUngapAlignment(reverse, MIN_IDENTITY_FOR_SCAFFOLD, positionBuffer, identity);
//            if (reversePosition.id == 0) {
                reversePosition = this->mapReadMultiSeedFiltered(reverse, positionBuffer);
                if (reversePosition.id == 0) continue;
//            }
            // if forward read and reverse one are mapped same contig
            if (forwardPosition.id == -reversePosition.id) {
                if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset)
                    insertLength = static_cast<long>(reversePosition.offset - forwardPosition.offset + 1);
                else if (reversePosition.id > 0 && reversePosition.offset < forwardPosition.offset)
                    insertLength = static_cast<long>(forwardPosition.offset - reversePosition.offset + 1);
                else
                    continue;

                if (insertLength < std::min(forward.length, reverse.length)) continue;

                fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
                ++numMappedSameContig;
            }
            else if (forwardPosition.id != reversePosition.id) {
                fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                ++numMappedDiffContig;
            }

            sum += forward.length + reverse.length;
        }
    }
/*
    if (numMappedSameContig == 0) {
        throw platanus::MapError("No read mapped in the same contig!!");
    }
    if (numMappedDiffContig == 0) {
        throw platanus::MapError("no read links contigs!!");
    }
*/
    if (numMappedDiffContig != 0)
		library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());
	else
		library[0].setAverageCoverage(1);

    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMappedSameContig + numMappedDiffContig << " (" << static_cast<double>(numMappedSameContig + numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_DIFFERENT_CONTIGS = " << numMappedDiffContig << " (" << static_cast<double>(numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}

void Mapper::mapTagPairMT(vector<SeqLib> &library, const long numThread)
{
    long totalRead = 0;
    long numMapped = 0;
    long sum = 0;

    cerr << "mapping tagged reads..." << endl;
    cerr << "Each dot below indicates " << DISPLAY_NUM_READ_UNIT/1000000 << "M reads processed." << endl;

    omp_set_num_threads(numThread);

    # pragma omp parallel for schedule(static, 1) reduction(+: sum, totalRead, numMapped)
    for (long i = 0; i < numThread; ++i) {
        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();

		platanus::SEQ forward;
		platanus::SEQ reverse;
		platanus::Position forwardPosition;
		platanus::Position reversePosition;
		int forwardTagID;
		int reverseTagID;
		vector<platanus::Position> positionBuffer;

		long NextDisplayNumPair = DISPLAY_NUM_READ_UNIT;
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
			fread(&forwardTagID, sizeof(int), 1, library[i].pairFP);
            reverse.readTemporaryFile(library[i].pairFP);
			fread(&reverseTagID, sizeof(int), 1, library[i].pairFP);
			totalRead += 2;

			if (i == 0 && totalRead*numThread >= NextDisplayNumPair) {
				cerr << '.' << std::flush;
				NextDisplayNumPair += DISPLAY_NUM_READ_UNIT;
			}

            if (forward.length >= seedLength) {
				forwardPosition = this->mapReadMultiSeedFiltered(forward, positionBuffer);
				if (forwardPosition.id != 0) {
					++numMapped;
					sum += forward.length;
					int contigIndex = abs(forwardPosition.id) - 1;
					fwrite(&forwardTagID, sizeof(int), 1, library[i].mappedFP);
					fwrite(&contigIndex, sizeof(int), 1, library[i].mappedFP);
				}
			}

            if (reverse.length >= seedLength) {
				reversePosition = this->mapReadMultiSeedFiltered(reverse, positionBuffer);
				if (reversePosition.id != 0) {
					++numMapped;
					sum += reverse.length;
					int contigIndex = abs(reversePosition.id) - 1;
					fwrite(&reverseTagID, sizeof(int), 1, library[i].mappedFP);
					fwrite(&contigIndex, sizeof(int), 1, library[i].mappedFP);
				}
			}
        }
    }

    if (numMapped != 0)
		library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());
	else
		library[0].setAverageCoverage(1);

    cerr << endl;
    cerr << "TOTAL_READ = " << totalRead << endl;
    cerr << "MAPPED_READ = " << numMapped << " (" << static_cast<double>(numMapped) / totalRead << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}

//added by ouchi
void Mapper::mapHiCPairMT(vector<SeqLib> &library, const long numThread)
{
    long totalPair = 0;
    long numMapped = 0;
    long sum = 0;

    cerr << "mapping HiC reads..." << endl;
    cerr << "Each dot below indicates " << DISPLAY_NUM_READ_UNIT/1000000 << "M reads processed." << endl;

    omp_set_num_threads(numThread);

    # pragma omp parallel for schedule(static, 1) reduction(+: sum, totalPair, numMapped)
    for (long i = 0; i < numThread; ++i) {
        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();

		platanus::SEQ forward;
		platanus::SEQ reverse;
		platanus::Position forwardPosition;
		platanus::Position reversePosition;
		vector<platanus::Position> positionBuffer;

		long NextDisplayNumPair = DISPLAY_NUM_READ_UNIT;
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
			++totalPair;

			if (i == 0 && totalPair*numThread >= NextDisplayNumPair) {
				cerr << '.' << std::flush;
				NextDisplayNumPair += DISPLAY_NUM_READ_UNIT;
			}

            if (forward.length < seedLength || reverse.length < seedLength) continue;

            forwardPosition = this->mapReadMultiSeedFiltered(forward, positionBuffer);
            if (forwardPosition.id == 0) continue;
            reversePosition = this->mapReadMultiSeedFiltered(reverse, positionBuffer);
            if (reversePosition.id == 0) continue;

            fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            ++numMapped;
            sum += forward.length + reverse.length;
        }
    }

    if (numMapped != 0)
		library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());
	else
		library[0].setAverageCoverage(1);

    cerr << endl;
    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMapped << " (" << static_cast<double>(numMapped) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}
//


void Mapper::mapPairAndSaveReadLink(vector<SeqLib> &library, const long minInsertion, const long minContigLength, const long numThread)
{
    vector<platanus::Position> forwardPositionBuffer;
    vector<platanus::Position> reversePositionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMapped = 0;
    long sum = 0;
    char tmpChar;

    cerr << "mapping reads..." << endl;
    cerr << "Each dot below indicates " << DISPLAY_NUM_READ_UNIT/1000000 << "M reads processed." << endl;

    omp_set_num_threads(numThread);

	unsigned i = 0;
    #    pragma omp parallel for  schedule(static, 1) private(i, forwardPositionBuffer, reversePositionBuffer,forward, reverse, insertLength) reduction(+: sum, totalPair, numMapped)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedSameContigFP != NULL)
            fclose(library[i].mappedSameContigFP);
        library[i].mappedSameContigFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();

        if (library[i].mappedReadFP != NULL)
            fclose(library[i].mappedReadFP);
        library[i].mappedReadFP = platanus::makeTemporaryFile();

		long NextDisplayNumPair = DISPLAY_NUM_READ_UNIT/2;
        rewind(library[i].pairFP);

		//std::cerr << "test" << std::endl; //added by ouchi for test

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;

		//std::cerr << totalPair << "\t" << std::endl; //added by ouchi for test

			if (i == 0 && totalPair*numThread >= NextDisplayNumPair) {
				cerr << '.' << std::flush;
				NextDisplayNumPair += DISPLAY_NUM_READ_UNIT/2;
			}

            if (forward.length < seedLength || reverse.length < seedLength) continue;

			this->mapReadMultiReportsMultiSeedFiltered(forward, forwardPositionBuffer, minContigLength);
//			if (forwardPositionBuffer.size() == 0)
//				this->mapReadUngapAlignment(forward, MIN_IDENTITY_TO_CHECK_MAPPING, forwardPositionBuffer, NULL);

			this->mapReadMultiReportsMultiSeedFiltered(reverse, reversePositionBuffer, minContigLength);
//			if (reversePositionBuffer.size() == 0)
//				this->mapReadUngapAlignment(reverse, MIN_IDENTITY_TO_CHECK_MAPPING, reversePositionBuffer, NULL);

			unsigned numPairLink = 0;
			unsigned numReadLink = 0;
			vector<platanus::Position> finalPositionBuffer;

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> mappedPairIDFlag;
			std::pair<int, int> pairID;

			for (unsigned j = 0; j < forwardPositionBuffer.size(); ++j) {
				for (unsigned k = 0; k < reversePositionBuffer.size(); ++k) {
					if (abs(forwardPositionBuffer[j].id) != abs(reversePositionBuffer[k].id)) {

						if (abs(forwardPositionBuffer[j].id) < abs(reversePositionBuffer[k].id))
							pairID = std::make_pair(abs(forwardPositionBuffer[j].id), -sign(forwardPositionBuffer[j].id) * reversePositionBuffer[k].id);
						else
							pairID = std::make_pair(abs(reversePositionBuffer[k].id), -sign(reversePositionBuffer[k].id) * forwardPositionBuffer[j].id);

						if (mappedPairIDFlag.find(pairID) == mappedPairIDFlag.end()) {
							finalPositionBuffer.push_back(forwardPositionBuffer[j]);
							finalPositionBuffer.push_back(reversePositionBuffer[k]);
							++numPairLink;
						}
					}
					else if (forwardPositionBuffer[j].id == -reversePositionBuffer[k].id) {
						if (forwardPositionBuffer[j].id > 0 && forwardPositionBuffer[j].offset < reversePositionBuffer[k].offset) {
							insertLength = static_cast<long>(reversePositionBuffer[k].offset - forwardPositionBuffer[j].offset + 1);
							if (insertLength >= std::min(forward.length, reverse.length)) {
								fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);

								fwrite(&(forwardPositionBuffer[j]), sizeof(platanus::Position), 1, library[i].mappedSameContigFP);
								fwrite(&insertLength, sizeof(long), 1, library[i].mappedSameContigFP);
							}
						}
						else if (reversePositionBuffer[k].id > 0 && reversePositionBuffer[k].offset < forwardPositionBuffer[j].offset) {
							insertLength = static_cast<long>(forwardPositionBuffer[j].offset - reversePositionBuffer[k].offset + 1);
							if (insertLength >= std::min(forward.length, reverse.length)) {
								fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);

								fwrite(&(reversePositionBuffer[k]), sizeof(platanus::Position), 1, library[i].mappedSameContigFP);
								fwrite(&insertLength, sizeof(long), 1, library[i].mappedSameContigFP);
							}
						}
					}
				}
			}


			if (forwardPositionBuffer.size() > 1) {
				for (unsigned j = 0; j < forwardPositionBuffer.size() - 1; ++j) {
					for (unsigned k = j + 1; k < forwardPositionBuffer.size(); ++k) {
						if (abs(forwardPositionBuffer[j].id) == abs(forwardPositionBuffer[k].id))
							continue;

						if (abs(forwardPositionBuffer[j].id) < abs(forwardPositionBuffer[k].id))
							pairID = std::make_pair(abs(forwardPositionBuffer[j].id), sign(forwardPositionBuffer[j].id) * forwardPositionBuffer[k].id);
						else
							pairID = std::make_pair(abs(forwardPositionBuffer[k].id), sign(forwardPositionBuffer[k].id) * forwardPositionBuffer[j].id);

						if (mappedPairIDFlag.find(pairID) == mappedPairIDFlag.end()) {
							mappedPairIDFlag.insert(pairID);
							finalPositionBuffer.push_back(forwardPositionBuffer[j]);
							finalPositionBuffer.push_back(forwardPositionBuffer[k]);
							++numReadLink;
						}
					}
				}
			}

			if (reversePositionBuffer.size() > 1) {
				for (unsigned j = 0; j < reversePositionBuffer.size() - 1; ++j) {
					for (unsigned k = j + 1; k < reversePositionBuffer.size(); ++k) {
						if (abs(reversePositionBuffer[j].id) == abs(reversePositionBuffer[k].id))
							continue;

						if (abs(reversePositionBuffer[j].id) < abs(reversePositionBuffer[k].id))
							pairID = std::make_pair(abs(reversePositionBuffer[j].id), sign(reversePositionBuffer[j].id) * reversePositionBuffer[k].id);
						else
							pairID = std::make_pair(abs(reversePositionBuffer[k].id), sign(reversePositionBuffer[k].id) * reversePositionBuffer[j].id);

						if (mappedPairIDFlag.find(pairID) == mappedPairIDFlag.end()) {
							finalPositionBuffer.push_back(reversePositionBuffer[j]);
							finalPositionBuffer.push_back(reversePositionBuffer[k]);
							++numReadLink;
						}
					}
				}
			}


			fwrite(&numPairLink, sizeof(unsigned), 1, library[i].mappedFP);
			fwrite(&numReadLink, sizeof(unsigned), 1, library[i].mappedFP);
			for (unsigned j = 0; j < finalPositionBuffer.size(); ++j)
				fwrite(&(finalPositionBuffer[j]), sizeof(platanus::Position), 1, library[i].mappedFP);

			if (!(forwardPositionBuffer.empty() && reversePositionBuffer.empty()))
				++numMapped;

            sum += forward.length + reverse.length;
        }
    }

	library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());

    for (unsigned i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);

//        rewind(library[i].mappedSameContigFP);
//        while (fread(&tmpChar, sizeof(char), 1, library[i].mappedSameContigFP))
//            putc(tmpChar, library[0].mappedSameContigFP);
//        fclose(library[i].mappedSameContigFP);
    }

    cerr << endl;
    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMapped << " (" << static_cast<double>(numMapped) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// map single read which overlap gaps
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapSmallGap(vector<SeqLib> &library, FILE *gapFP, const long numThread)
{
    FILE *tmpFP[numThread];

    std::cerr << "mapping reads that cover small gaps..." << std::endl;
    cerr << "Each dot below indicates " << DISPLAY_NUM_READ_UNIT/1000000 << "M reads processed." << endl;

    omp_set_num_threads(numThread);

    # pragma omp parallel for schedule (static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
        char gapSeq[platanus::ConstParam::MAX_READ_LEN];
        long gapSeqLength;
		long totalRead = 0;
        vector<platanus::Position> leftBuffer;
        vector<platanus::Position> rightBuffer;
        platanus::Position leftResult;
        platanus::Position rightResult;
        platanus::SEQ read;
        Kmer31 leftKmer;
        Kmer31 rightKmer;

        tmpFP[threadID] = platanus::makeTemporaryFile();

		long NextDisplayNumPair = DISPLAY_NUM_READ_UNIT;
        rewind(library[threadID].pairFP);

        while (read.readTemporaryFile(library[threadID].pairFP)) {
			++totalRead;
			if (threadID == 0 && totalRead*numThread >= NextDisplayNumPair) {
				cerr << '.' << std::flush;
				NextDisplayNumPair += DISPLAY_NUM_READ_UNIT;
			}

            if (read.length < 2 * this->seedLength) continue;
            long i;
            long start = 0;
            long end = 0;
            for (i = 0; i < this->keyLength; ++i) {
                if (read.base[i] == 4 || read.base[read.length - this->seedLength + i] == 4) break;
                leftKmer.forward = (leftKmer.forward << 2) | static_cast<unsigned long long>(read.base[i]);
                leftKmer.reverse = (leftKmer.reverse << 2) | static_cast<unsigned long long>((0x3^read.base[this->keyLength - i - 1]));
                rightKmer.forward = (rightKmer.forward << 2) | static_cast<unsigned long long>(read.base[read.length - this->seedLength + i]);
                rightKmer.reverse = (rightKmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[read.length - this->seedLength + this->keyLength - i - 1]);
            }
            if (i != this->keyLength) continue;

            this->mapReadGapCloseLeft(read, leftKmer, leftBuffer);
            this->mapReadGapCloseRight(read, rightKmer, rightBuffer);


            leftResult.id = 0;
            rightResult.id = 0;

            for (i = 0; i < static_cast<long>(leftBuffer.size()); ++i) {
                unsigned long long j = 0;
                for (; j < rightBuffer.size(); ++j) {
                    if (leftBuffer[i].id == 0
                        || rightBuffer[j].id == 0
                        || leftBuffer[i].id != rightBuffer[j].id
                        || (leftBuffer[i].id > 0 && (leftBuffer[i].offset >= rightBuffer[j].offset || rightBuffer[j].offset - leftBuffer[i].offset + this->seedLength > 2 * read.length))
                        || (leftBuffer[i].id < 0 && (leftBuffer[i].offset <= rightBuffer[j].offset || leftBuffer[i].offset - rightBuffer[j].offset + this->seedLength > 2 * read.length))) continue;

                    if (leftResult.id != 0) break;
                    leftResult = leftBuffer[i];
                    rightResult = rightBuffer[j];
                }
                if (j != rightBuffer.size()) break;
            }
            if (static_cast<unsigned long long>(i) != leftBuffer.size() || leftResult.id == 0) continue;
            if (leftResult.id > 0) {
                start = searchGapStart(leftResult, rightResult, leftResult);
                if (start == 0) continue;
                end = searchGapEnd(leftResult, rightResult, leftResult, read.length);

                if (!checkGapContinuious(leftResult, rightResult, leftResult, start, end, read.length)) continue;

                gapSeqLength = 0;
                for (i = start; i < end; ++i) {
                    gapSeq[gapSeqLength] = read.base[i];
                    ++gapSeqLength;
                }
                leftResult.offset += start;
            } else {
                start = searchGapStart(rightResult, leftResult, leftResult);
                if (start == 0) continue;
                end = searchGapEnd(rightResult, leftResult, leftResult, read.length);

                if (!checkGapContinuious(rightResult, leftResult, leftResult, start, end, read.length)) continue;

                gapSeqLength = 0;
                for (i = read.length - start - 1; i >= read.length - end; --i) {
                    gapSeq[gapSeqLength] = 0x3 ^ read.base[i];
                    ++gapSeqLength;
                }
                leftResult.offset = rightResult.offset + start;
                leftResult.id *= -1;
            }

            gapSeqLength = end - start;
            fwrite(&leftResult, sizeof(platanus::Position), 1, tmpFP[threadID]);
            fwrite(&gapSeqLength, sizeof(long), 1, tmpFP[threadID]);
            if (gapSeqLength > 0)
                fwrite(gapSeq, sizeof(char), gapSeqLength, tmpFP[threadID]);
        }
    }
    cerr << endl;

    char c;
    fseek(gapFP, 0, SEEK_END);
    for (long threadID = 0; threadID < numThread; ++threadID) {
        rewind(tmpFP[threadID]);
        while (fread(&c, sizeof(char), 1, tmpFP[threadID])) {
            putc(c, gapFP);
        }
        fclose(tmpFP[threadID]);
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// wrapper function
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapReadGapCloseLeft(const platanus::SEQ &read, const Kmer31 &kmer, vector<platanus::Position> &buffer)
{
    this->mapReadGapClose(read, kmer, buffer, true);
}


//////////////////////////////////////////////////////////////////////////////////////
// wrapper function
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapReadGapCloseRight(const platanus::SEQ &read, const Kmer31 &kmer, vector<platanus::Position> &buffer)
{
    this->mapReadGapClose(read, kmer, buffer, false);
}


//////////////////////////////////////////////////////////////////////////////////////
// map single read actual for gap_close
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::mapReadGapClose(const platanus::SEQ &read, const Kmer31 &kmer, vector<platanus::Position> &buffer, const bool isLeft=true)
{
    long diff;
    long bufferSize = this->mapKeyWithLimit(kmer, buffer);
    if (isLeft) {
        diff = 0;
    } else {
        diff = read.length - this->seedLength;
    }

    for (long i = 0; i < bufferSize; ++i) {
        long j = this->keyLength;
        if (buffer[i].id > 0) {
            if (buffer[i].offset > this->seq[buffer[i].id - 1].length - this->seedLength) continue;
            for (; j < this->seedLength; ++j) {
                if (this->seq[buffer[i].id - 1].base[buffer[i].offset + j] != read.base[diff + j]) break;
            }
        } else {
            if (buffer[i].offset < this->seedLength - this->keyLength) continue;
            for (; j < this->seedLength; ++j) {
                if (this->seq[-(buffer[i].id) - 1].base[buffer[i].offset + this->keyLength - j - 1] != (0x3^read.base[diff + j])) break;
            }
            buffer[i].offset -= this->seedLength - this->keyLength;
        }
        if (j != this->seedLength) {
            buffer[i].id = 0;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
long Mapper::searchGapStart(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult)
{
    long start = 0;
    for (long i = result1.offset + this->seedLength; i < result2.offset; ++i) {
        if (this->seq[abs(leftResult.id) - 1].base[i] == 4) {
            start = i - result1.offset;
            break;
        }
    }
    return start;
}



//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
long Mapper::searchGapEnd(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult, const long length)
{
    long end = 0;
    for (long i = result2.offset - 1; i >= result1.offset + this->seedLength; --i) {
        if (this->seq[abs(leftResult.id) - 1].base[i] == 4) {
            end = length - (result2.offset + this->seedLength - 1 - i);
            break;
        }
    }
    return end;
}



//////////////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////////////
bool Mapper::checkGapContinuious(const platanus::Position &result1, const platanus::Position &result2, const platanus::Position &leftResult, const long start, const long end, const long length)
{
    long i;
    for (i = result1.offset + start; i < result2.offset - (length - end - this->seedLength); ++i) {
        if (this->seq[abs(leftResult.id) - 1].base[i] != 4) break;
    }
    return i == result2.offset - (length - end - this->seedLength);
}



//////////////////////////////////////////////////////////////////////////////////////
// map each read and count the number of read mapped same contig
// for gap close function
//////////////////////////////////////////////////////////////////////////////////////
void Mapper::gatherPairReadMappedSameContig(vector<SeqLib> &library, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;

    cerr << "mapping reads..." << endl;
    cerr << "Each dot below indicates " << DISPLAY_NUM_READ_UNIT/1000000 << "M reads processed." << endl;

    omp_set_num_threads(numThread);


    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: totalPair, numMappedSameContig)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);
		long NextDisplayNumPair = DISPLAY_NUM_READ_UNIT/2;

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;

			if (i == 0 && totalPair*numThread >= NextDisplayNumPair) {
				cerr << '.' << std::flush;
				NextDisplayNumPair += DISPLAY_NUM_READ_UNIT/2;
			}

            if (forward.length < seedLength || reverse.length < seedLength) continue;

            forwardPosition = this->mapReadMultiSeedFiltered(forward, positionBuffer);
            reversePosition = this->mapReadMultiSeedFiltered(reverse, positionBuffer);
//			double identity;
//            forwardPosition = this->mapReadUngapAlignment(forward, MIN_IDENTITY_FOR_GAP_CLOSE, positionBuffer, identity);
//            reversePosition = this->mapReadUngapAlignment(reverse, MIN_IDENTITY_FOR_GAP_CLOSE, positionBuffer, identity);
            if (forwardPosition.id == 0 && reversePosition.id == 0) continue;

            fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            for (int j = 0; j < forward.numUnknown; ++j)
                forward.base[forward.positionUnknown[j]] = 4;
            fwrite(&(forward.length), sizeof(long), 1, library[i].mappedFP);
            fwrite(forward.base.c_str(), sizeof(char), forward.length, library[i].mappedFP);

            fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            for (int j = 0; j < reverse.numUnknown; ++j)
                reverse.base[reverse.positionUnknown[j]] = 4;
            fwrite(&(reverse.length), sizeof(long), 1, library[i].mappedFP);
            fwrite(reverse.base.c_str(), sizeof(char), reverse.length, library[i].mappedFP);

            if (forwardPosition.id == -reversePosition.id) {
                if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset) {
                    insertLength = static_cast<long>(reversePosition.offset) - forwardPosition.offset;
                } else if (reversePosition.id > 0 && forwardPosition.offset > reversePosition.offset) {
                    insertLength = static_cast<long>(forwardPosition.offset) - reversePosition.offset;
                } else
                    continue;

                if (insertLength < std::min(forward.length, reverse.length)) continue;
                fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
                ++numMappedSameContig;
            }
        }
    }

/*
    if (numMappedSameContig == 0) {
        throw platanus::MapError("no read mapped in the same contig!!");
    }
*/

    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << endl;
    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// decide the contig where the variable read is mapped
// contig contains bubble seq
//////////////////////////////////////////////////////////////////////////////////////
platanus::Position HeteroMapper::mapRead(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer, const long wordLength)
{
    long numResult = 0;
    vector<platanus::Position> result(read.length/wordLength + 2);
    Kmer31 kmer;
    long firstMax, secondMax;
    int i, j, k = 0;
    for (i = read.length - wordLength; i > -(wordLength); i -= wordLength) {
        if (i < 0)
            i = 0;
        kmer.forward = kmer.reverse = 0;
        for (j = 0; j < keyLength; ++j) {
            if (read.base[i + j] == 4) break;
            kmer.forward = (kmer.forward << 2) | static_cast<unsigned long long>(read.base[i + j]);
            kmer.reverse = (kmer.reverse << 2) | static_cast<unsigned long long>(0x3^read.base[i + keyLength - j - 1]);
        }
        if (j != keyLength) continue;
        contigMap.mapKeyWithLimit(kmer, positionBuffer);
        result[numResult].id = 0;
        auto it = positionBuffer.begin();
        auto end = positionBuffer.end();
        for (; it != end; ++it) {
            if (it->id > 0) {
                if (it->offset > contigMap.getSeqLength(it->id - 1) - wordLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (contigMap.getBase(it->id - 1, it->offset + k) != read.base[i + k])
                        break;
                }
                it->offset -= i;
            }
            else {
                if (it->offset < wordLength - keyLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (contigMap.getBase(-(it->id) - 1, it->offset + keyLength - k - 1) != (3^read.base[i + k]))
                        break;
                }
                it->offset += i + keyLength - 1;
            }

            if (k == wordLength) {
                if (result[numResult].id == 0)
                    result[numResult] = *it;
                else
                    break;
            }
        }
        if (it == end) {
            if (result[numResult].id != 0) {
                ++numResult;
                continue;
            }
        }
        else
            continue;

        bubbleMap.mapKeyWithLimit(kmer, positionBuffer);
        result[numResult].id = 0;

        it = positionBuffer.begin();
        end = positionBuffer.end();

        for (; it != end; ++it) {
            if (it->id > 0) {
                if (it->offset > bubbleMap.getSeqLength(it->id - 1) - wordLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (bubbleMap.getBase(it->id - 1, it->offset + k) != read.base[i + k])
                        break;
                }
                it->offset -= i;
                if (bubblePosition[it->id - 1].id > 0) {
                    it->offset = bubblePosition[it->id - 1].offset + it->offset - bubbleMap.getSeqLength(it->id - 1) / 2;
                    it->id = bubblePosition[it->id - 1].id;
                }
                else {
                    it->offset = bubblePosition[it->id - 1].offset - it->offset + bubbleMap.getSeqLength(it->id - 1) / 2;
                    it->id = bubblePosition[it->id - 1].id;
                }
            }
            else {
                if (it->offset < wordLength - keyLength)
                    continue;
                for (k = keyLength; k < wordLength; ++k) {
                    if (bubbleMap.getBase(-(it->id) - 1, it->offset + keyLength - k - 1) != (0x3^read.base[i + k]))
                        break;
                }
                it->offset += i + keyLength - 1;
                if (bubblePosition[-(it->id) - 1].id > 0) {
                    it->offset = bubblePosition[-(it->id) - 1].offset + it->offset - bubbleMap.getSeqLength(-(it->id) - 1) / 2;
                    it->id = -(bubblePosition[-(it->id) - 1].id);
                }
                else {
                    it->offset = bubblePosition[-(it->id) - 1].offset - it->offset + bubbleMap.getSeqLength(-(it->id) - 1) / 2;
                    it->id = -(bubblePosition[-(it->id) - 1].id);
                }
            }

            if (k == wordLength) {
                if (result[numResult].id == 0)
                    result[numResult] = *it;
                else
                    break;
            }
        }
        if (it == end && result[numResult].id != 0) {
            ++numResult;
            break;
        }
    }
    if (numResult == 0) {
        result[0].id = 0;
        return result[0];
    }
    result.resize(numResult);
    sort(result.begin(), result.end());

    i = firstMax = secondMax = 0;
    result.emplace_back(platanus::Position());
    while (result[i].id != 0) {
        j = i;
        while (result[i].id == result[i+1].id && result[i].offset == result[i+1].offset)
            ++i;
        ++i;
        if (i - j >= firstMax) {
            secondMax = firstMax;
            firstMax = i - j;
            k = j;
        }
    }

    if (firstMax == secondMax) {
        result[0].id = 0;
        return result[0];
    }

    return result[k];
}



//////////////////////////////////////////////////////////////////////////////////////
// map each read and count the number of read mapped same contig (contain bubbles) and different one
//////////////////////////////////////////////////////////////////////////////////////
void HeteroMapper::mapPairMT(vector<SeqLib> &library, const long minInsertion, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    long numMappedDiffContig = 0;
    long sum = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;
    if (bubbleMap.getNumSeq() == 0) {
        contigMap.mapPairMT(library, minInsertion, numThread);
        return;
    }

    cerr << "mapping reads..." << endl;

    omp_set_num_threads(numThread);
    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: sum, totalPair, numMappedSameContig, numMappedDiffContig)
    for (i = 0; i < numThread; ++i) {

        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;
            if (forward.length < seedLength || reverse.length < seedLength) continue;
            forwardPosition = contigMap.mapRead(forward, positionBuffer, this->seedLength);
//			double identity;
//            forwardPosition = contigMap.mapReadUngapAlignment(forward, MIN_IDENTITY_FOR_SCAFFOLD, positionBuffer, identity);
            if (forwardPosition.id == 0) {
                forwardPosition = this->mapRead(forward, positionBuffer, this->seedLength);
                if (forwardPosition.id == 0) continue;
            }

            reversePosition = contigMap.mapRead(reverse, positionBuffer, this->seedLength);
//            reversePosition = contigMap.mapReadUngapAlignment(reverse, MIN_IDENTITY_FOR_SCAFFOLD, positionBuffer, identity);
            if (reversePosition.id == 0) {
                reversePosition = this->mapRead(reverse, positionBuffer, this->seedLength);
                if (reversePosition.id == 0) continue;
            }


            if (forwardPosition.id == 0 || reversePosition.id == 0) continue;
            if (forwardPosition.id == -reversePosition.id) {
                if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset)
                    insertLength = static_cast<long>(reversePosition.offset - forwardPosition.offset + 1);
                else if (reversePosition.id > 0 && reversePosition.offset < forwardPosition.offset)
                    insertLength = static_cast<long>(forwardPosition.offset - reversePosition.offset + 1);
                else
                    continue;

                if (insertLength < minInsertion || insertLength < std::min(forward.length, reverse.length))
                    continue;

                fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
                ++numMappedSameContig;
            }
            else if (forwardPosition.id != reversePosition.id) {
                fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
                ++numMappedDiffContig;
            }
            sum += forward.length + reverse.length;
        }
    }

    if (numMappedSameContig == 0) {
        throw platanus::MapError("no read mapped in the same contig!!");
    }
    if (numMappedDiffContig == 0) {
        throw platanus::MapError("no read links contigs!!");
    }

    library[0].setAverageCoverage(static_cast<double>(sum) / contigMap.getSeqPoolSize());
    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }
    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMappedSameContig + numMappedDiffContig << " (" << static_cast<double>(numMappedSameContig + numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_DIFFERENT_CONTIGS = " << numMappedDiffContig << " (" << static_cast<double>(numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}


void Mapper::mapPairToCalculateCoverage(vector<SeqLib> &library, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ forward;
    platanus::SEQ reverse;
    long totalPair = 0;
    long insertLength = 0;
    long numMappedSameContig = 0;
    long numMappedDiffContig = 0;
    long sum = 0;
    platanus::Position forwardPosition;
    platanus::Position reversePosition;
    char tmpChar;
    int i;

    cerr << "mapping reads..." << endl;

    omp_set_num_threads(numThread);

    #    pragma omp parallel for  schedule(static, 1) private(positionBuffer, forward, reverse, forwardPosition, reversePosition, insertLength) reduction(+: sum, totalPair, numMappedSameContig, numMappedDiffContig)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedFP != NULL)
            fclose(library[i].mappedFP);
        library[i].mappedFP = platanus::makeTemporaryFile();
        rewind(library[i].pairFP);

        while (forward.readTemporaryFile(library[i].pairFP)) {
            reverse.readTemporaryFile(library[i].pairFP);
            ++totalPair;

			double forwardIdentity;
            if (forward.length >= seedLength) {
				forwardPosition = this->mapReadUngapAlignment(forward, 0.0, positionBuffer, &forwardIdentity);
				if (forwardPosition.id != 0) {
					fwrite(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP);
					for (int j = 0; j < forward.numUnknown; ++j)
						forward.base[forward.positionUnknown[j]] = 4;
					fwrite(&(forward.length), sizeof(long), 1, library[i].mappedFP);
					fwrite(forward.base.c_str(), sizeof(char), forward.length, library[i].mappedFP);
					fwrite(&forwardIdentity, sizeof(double), 1, library[i].mappedFP);
				}
			}

			double reverseIdentity;
            if (reverse.length >= seedLength) {
				reversePosition = this->mapReadUngapAlignment(reverse, 0.0, positionBuffer, &reverseIdentity);
				if (reversePosition.id != 0) {
					fwrite(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
					for (int j = 0; j < reverse.numUnknown; ++j)
						reverse.base[reverse.positionUnknown[j]] = 4;
					fwrite(&(reverse.length), sizeof(long), 1, library[i].mappedFP);
					fwrite(reverse.base.c_str(), sizeof(char), reverse.length, library[i].mappedFP);
					fwrite(&reverseIdentity, sizeof(double), 1, library[i].mappedFP);
				}
			}

            if (forwardPosition.id != 0 && reversePosition.id != 0) {
				if (forwardPosition.id == -reversePosition.id) {
					if (forwardPosition.id > 0 && forwardPosition.offset < reversePosition.offset)
						insertLength = static_cast<long>(reversePosition.offset - forwardPosition.offset + 1);
					else if (reversePosition.id > 0 && reversePosition.offset < forwardPosition.offset)
						insertLength = static_cast<long>(forwardPosition.offset - reversePosition.offset + 1);
					else
						continue;

					if (insertLength < std::min(forward.length, reverse.length)) continue;

					fwrite(&insertLength, sizeof(long), 1, library[i].insertLengthFP);
					++numMappedSameContig;
				}
				else if (forwardPosition.id != reversePosition.id) {
					++numMappedDiffContig;
				}
			}

            sum += forward.length + reverse.length;
        }
    }

    if (numMappedDiffContig != 0)
		library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());
	else
		library[0].setAverageCoverage(1);

    for (i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << "TOTAL_PAIR = " << totalPair << endl;
    cerr << "MAPPED_PAIR = " << numMappedSameContig + numMappedDiffContig << " (" << static_cast<double>(numMappedSameContig + numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_DIFFERENT_CONTIGS = " << numMappedDiffContig << " (" << static_cast<double>(numMappedDiffContig) / totalPair << ")" << endl;
    cerr << "MAPPED_IN_SAME_CONTIG = " << numMappedSameContig << " (" << static_cast<double>(numMappedSameContig) / totalPair << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}

void Mapper::filterMappedPositionUngapAlignment(const platanus::SEQ &read, vector<platanus::Position> &positionBuffer, const double minIdentity)
{
    unsigned long i;
    unsigned long numHit = 0;

	for (i = 0; i < positionBuffer.size(); ++i) {
		if (judgeMappedPositionUngapAlignment(read, positionBuffer[i], minIdentity)) {
			if (i > numHit)
				positionBuffer[numHit] = positionBuffer[i];
			++numHit;
		}
	}

	positionBuffer.resize(numHit);
}

bool Mapper::judgeMappedPositionUngapAlignment(const platanus::SEQ &read, const platanus::Position readLeftPosition, const double minIdentity)
{
	if (readLeftPosition.id == 0)
		return false;

    long k = 0;
	long numMis = 0;
    long alignmentStart;
    long alignmentEnd;

	if (readLeftPosition.id > 0) {
		platanus::SEQ &target = seq[readLeftPosition.id - 1];
		alignmentStart = -std::min(static_cast<long>(readLeftPosition.offset), 0L);
		alignmentEnd = std::min(read.length, target.length - readLeftPosition.offset);
		long misThreshold = (int64_t)((1.0 - minIdentity) * (alignmentEnd - alignmentStart));

		for (k = alignmentStart; k < alignmentEnd; ++k) {
			if (target.base[readLeftPosition.offset + k] != 4 && read.base[k] != 4 && target.base[readLeftPosition.offset + k] != read.base[k]) {
				++numMis;
				if (numMis > misThreshold)
					break;
			}
		}
	}
	else {
		platanus::SEQ &target = seq[-(readLeftPosition.id) - 1];
		alignmentStart = -std::min(0L, target.length - readLeftPosition.offset);
		alignmentEnd = std::min(read.length, static_cast<long>(readLeftPosition.offset));
		long misThreshold = static_cast<long>((1.0 - minIdentity) * (alignmentEnd - alignmentStart));

		for (k = alignmentStart; k < alignmentEnd; ++k) {
			if (target.base[readLeftPosition.offset - k] != 4 && read.base[k] != 4 && target.base[readLeftPosition.offset - k] != (3^read.base[k])) {
				++numMis;
				if (numMis > misThreshold)
					break;
			}
		}
	}

	if (k == alignmentEnd)
		return true;
	else
		return false;
}

void Mapper::mapLongReadAndSaveReadLink(vector<SeqLib> &library, const long minContigLength, const long numThread)
{
    vector<platanus::Position> positionBuffer;
    platanus::SEQ read;
    long totalRead = 0;
    long numMapped = 0;
    long sum = 0;
	const long match = this->seedLength;

    cerr << "mapping long reads..." << endl;

    omp_set_num_threads(numThread);

	unsigned i = 0;
    #    pragma omp parallel for  schedule(static, 1) private(i, positionBuffer, read) reduction(+: sum, totalRead, numMapped)
    for (i = 0; i < numThread; ++i) {
        if (library[i].insertLengthFP != NULL)
            fclose(library[i].insertLengthFP);
        library[i].insertLengthFP = platanus::makeTemporaryFile();

        if (library[i].mappedReadFP != NULL)
            fclose(library[i].mappedReadFP);
        library[i].mappedReadFP = platanus::makeTemporaryFile();

        rewind(library[i].pairFP);

        while (read.readTemporaryFile(library[i].pairFP)) {
            ++totalRead;
            if (read.length < seedLength)
				continue;

			this->mapReadMultiReportsMultiSeed(read, positionBuffer, minContigLength);
			if (positionBuffer.size() == 0)
				continue;

			unsigned numPosition = positionBuffer.size();
			fwrite(&numPosition, sizeof(unsigned), 1, library[i].mappedReadFP);
			for (unsigned j = 0; j < positionBuffer.size(); ++j) {
				fwrite(&(positionBuffer[j]), sizeof(platanus::Position), 1, library[i].mappedReadFP);
				fwrite(&match, sizeof(long), 1, library[i].mappedReadFP);
			}

			fwrite(&(read.length), sizeof(long), 1, library[i].insertLengthFP);

			++numMapped;
            sum += read.length;
        }
    }

	library[0].setAverageCoverage(static_cast<double>(sum) / getSeqPoolSize());

    char tmpChar;
    for (unsigned i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << "TOTAL_READ = " << totalRead << endl;
    cerr << "MAPPED_READ = " << numMapped << " (" << static_cast<double>(numMapped) / totalRead << ")" << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}

void Mapper::reduceAlignmentsGreedy(vector<LongReadAlignment> &alignments, const long readLength, const long tolerenceLength)
{
	std::sort(alignments.begin(), alignments.end(), Mapper::LongReadAlignment::LongReadAlignmentGreater());

	unsigned long numRetained = alignments.size();
	for (unsigned long i = 1; i < alignments.size(); ++i) {
		for (unsigned long j = 0; j < i; ++j) {
			if (alignments[j].score != LONG_MIN && std::min(alignments[i].readEnd - alignments[j].readStart, alignments[j].readEnd - alignments[i].readStart) > tolerenceLength) {
				alignments[i].score = LONG_MIN;
				--numRetained;
				break;
			}
		}
	}

	std::partial_sort(alignments.begin(), alignments.begin() + numRetained, alignments.end(), Mapper::LongReadAlignment::LongReadAlignmentGreater());
	alignments.resize(numRetained);
}

void Mapper::readLongReadPAFfileAndSaveLink(const std::string PAFFilename, vector<SeqLib> &library, const long minAlignmentLength, const double minCoverage, const double minIdentity, const long tolerenceLength, const long numThread)
{
    cerr << "reading result of aligner ..." << endl;

	std::string oneLine, preQName(""), qName, tName, strand;
	long qLength, qStart, qEnd, tLength, tStart, tEnd, match, readLength=0;
	long threadIndex = 0, sumQAlignment = 0, sumTAlignment = 0;

	std::ifstream ifs;

	ifs.open(PAFFilename);
    while (ifs && getline(ifs, oneLine)) {
		std::istringstream ss(oneLine);
		ss >> qName >> qLength >> qStart >> qEnd >> strand >> tName >> tLength >> tStart >> tEnd >> match;
		sumQAlignment += qEnd - qStart;
		sumTAlignment += tEnd - tStart;
	}
	ifs.close();

	double readInsertionRate = (double)sumQAlignment / sumTAlignment;


	vector<LongReadAlignment> alignmentBuffer;
    long numMapped = 0;
    long sumMappedReadLength = 0;

	vector<FILE*> bufferFP(numThread);
    for (unsigned i = 0; i < numThread; ++i)
        bufferFP[i] = platanus::makeTemporaryFile();

	ifs.open(PAFFilename);
	while (1) {
		getline(ifs, oneLine);
		std::istringstream ss(oneLine);
		ss >> qName >> qLength >> qStart >> qEnd >> strand >> tName >> tLength >> tStart >> tEnd >> match;

		if (ifs.eof() || preQName != qName) {
			if (!alignmentBuffer.empty()) {
				fwrite(&readLength, sizeof(long), 1, bufferFP[threadIndex]);

				size_t bufferSize = alignmentBuffer.size();
				fwrite(&bufferSize, sizeof(size_t), 1, bufferFP[threadIndex]);
				for (unsigned j = 0; j < alignmentBuffer.size(); ++j)
					fwrite(&(alignmentBuffer[j]), sizeof(LongReadAlignment), 1, bufferFP[threadIndex]);

				threadIndex = (threadIndex + 1) %  numThread;
			}

			alignmentBuffer.clear();
			alignmentBuffer.clear();
			++numMapped;
			sumMappedReadLength += readLength;

			if (ifs.eof())
				break;
		}

		int alignmentLength = std::max(qEnd - qStart, tEnd - tStart);

		if ((double)match / alignmentLength >= minIdentity && (alignmentLength >= minAlignmentLength || (double)alignmentLength / std::min(qLength, tLength) >= minCoverage)) {
			unsigned contigIndex = this->nameIndex[tName];
			alignmentBuffer.resize(alignmentBuffer.size() + 1);

			alignmentBuffer.back().score = match;
			alignmentBuffer.back().readStart = qStart;
			alignmentBuffer.back().readEnd = qEnd;
			if (strand == "+") {
				alignmentBuffer.back().targetPosition.id = contigIndex + 1;
				alignmentBuffer.back().targetPosition.offset = tStart - qStart/readInsertionRate;
				alignmentBuffer.back().targetStart = tStart;
				alignmentBuffer.back().targetEnd = tEnd;
			}
			else {
				alignmentBuffer.back().targetPosition.id = -(contigIndex + 1);
				alignmentBuffer.back().targetPosition.offset = (tEnd - 1) + qStart/readInsertionRate;
				alignmentBuffer.back().targetStart = tEnd - 1;
				alignmentBuffer.back().targetEnd = tStart - 1;
			}

		}

		readLength = qLength;
		preQName = qName;
	}
	ifs.close();


    omp_set_num_threads(numThread);
    #pragma omp parallel for schedule(static, 1) private(alignmentBuffer, readLength)
	for (unsigned threadIndex = 0; threadIndex < numThread; ++threadIndex) {
        if (library[threadIndex].insertLengthFP != NULL)
            fclose(library[threadIndex].insertLengthFP);
        library[threadIndex].insertLengthFP = platanus::makeTemporaryFile();

        if (library[threadIndex].mappedReadFP != NULL)
            fclose(library[threadIndex].mappedReadFP);
        library[threadIndex].mappedReadFP = platanus::makeTemporaryFile();

		rewind(bufferFP[threadIndex]);
		while (fread(&readLength, sizeof(long), 1, bufferFP[threadIndex])) {
			size_t bufferSize;
			fread(&bufferSize, sizeof(size_t), 1, bufferFP[threadIndex]);
			alignmentBuffer.resize(bufferSize);
			for (unsigned j = 0; j < alignmentBuffer.size(); ++j)
				fread(&(alignmentBuffer[j]), sizeof(LongReadAlignment), 1, bufferFP[threadIndex]);

			reduceAlignmentsGreedy(alignmentBuffer, readLength, tolerenceLength);

			unsigned numAlignment = alignmentBuffer.size();
			fwrite(&numAlignment, sizeof(unsigned), 1, library[threadIndex].mappedReadFP);
			for (unsigned j = 0; j < alignmentBuffer.size(); ++j) {
				fwrite(&(alignmentBuffer[j].targetPosition), sizeof(platanus::Position), 1, library[threadIndex].mappedReadFP);
				fwrite(&(alignmentBuffer[j].targetStart), sizeof(int), 1, library[threadIndex].mappedReadFP);
				fwrite(&(alignmentBuffer[j].targetEnd), sizeof(int), 1, library[threadIndex].mappedReadFP);
				fwrite(&(alignmentBuffer[j].score), sizeof(long), 1, library[threadIndex].mappedReadFP);
			}

			fwrite(&readLength, sizeof(long), 1, library[threadIndex].insertLengthFP);
		}

		fclose(bufferFP[threadIndex]);
	}


	library[0].setAverageCoverage(static_cast<double>(sumMappedReadLength) / getSeqPoolSize());
    char tmpChar;
    for (unsigned i = 1; i < numThread; ++i) {
        rewind(library[i].insertLengthFP);
        while (fread(&tmpChar, sizeof(char), 1, library[i].insertLengthFP))
            putc(tmpChar, library[0].insertLengthFP);
        fclose(library[i].insertLengthFP);
    }

    cerr << "MAPPED_READ = " << numMapped << endl;
    cerr << "AVERAGE_COVERAGE = " << library[0].getAverageCoverage() << endl;
}
