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

#include "pairedDBG.h"
#include <omp.h>
#include <climits>
#include <cfloat>
#include <string>
#include <sstream>
#include <queue>
#include <array>
#include <numeric>

using std::vector;
using std::string;
using std::pair;
using std::cerr;
using std::endl;


// for node-state
const unsigned PairedDBG::DBG_HETERO = 0x8;
const unsigned PairedDBG::DBG_PRIMARY_BUBBLE = 0x10;
const unsigned PairedDBG::DBG_SECONDARY_BUBBLE = 0x20;

// for edge-state
const char PairedDBG::DBG_OVERLAP = 0x1;

// for contig-state
const char PairedDBG::DBG_CONTIG_BUBBLE_JUNCTION = 0x1;
const char PairedDBG::DBG_CONTIG_PRIMARY_BUBBLE = 0x2;
const char PairedDBG::DBG_CONTIG_SECONDARY_BUBBLE = 0x4;

// related to scaffolding modes
const unsigned PairedDBG::OVERLAP_MODE = 0x1;
const unsigned PairedDBG::PAIRED_END_LINK_MODE = 0x2;
const unsigned PairedDBG::LENGTH_CUTOFF_MODE = 0x4;
const unsigned PairedDBG::NON_DBG_OVERLAP_MODE = 0x8;
const unsigned PairedDBG::LONG_READ_LINK_MODE = 0x10;
const unsigned PairedDBG::BUBBLE_AWARE_MODE = 0x20;
const unsigned PairedDBG::SECONDARY_BUBBLE_REMOVAL_MODE = 0x40;
const unsigned PairedDBG::PREVIOUS_DIVISION_AWARE_MODE = 0x80;

const long PairedDBG::MAX_ITERATION_OF_CROSS_SOLUTION = 5;

const double PairedDBG::NON_PAIRED_END_TOLERENCE_FACTOR = 0.1;
const double PairedDBG::HETERO_COVERAGE_THRESHOLD_FACTOR = 1.75;
const double PairedDBG::HETERO_FORK_COVERAGE_THRESHOLD_FACTOR = 1.25;
const double PairedDBG::CROSS_LINK_RATE_THRESHOLD = 0.25;
const double PairedDBG::CROSS_SCORE_RATE_THRESHOLD = 0.5;
const double PairedDBG::MIN_BUBBLE_COUNT_FACTOR = 10000.0;

//added by ouchi
const double PairedDBG::HiC_CONTACTMAP_S_THRESHOLD = 1.3;
const double PairedDBG::HiC_CONTACTMAP_T_THRESHOLD = 1.8;
//

void PairedDBG::storeGraphLinkFromOverlap(vector<GraphLinkWithFlag> &graphLinkPool)
{
    # pragma omp parallel for schedule(dynamic)
    for (unsigned d = 0; d < TABLE_DIVID; ++d) {
		auto overlapIterator = overlapTable[d].begin();
		auto overlapEnd = overlapTable[d].end();
		for (; overlapIterator != overlapEnd; ++overlapIterator) {
			GraphLinkWithFlag overlapLink;
			overlapLink.overlapFlag = true;
			overlapLink.id1 = overlapIterator->second.id1;
			overlapLink.id2 = overlapIterator->second.id2;

            long i = id2Index(overlapLink.id1);
            if (contigPositionInScaffold[i].size() == 1) //added by ouchi (contigPositionInScaffold)
	            overlapLink.id1 = overlapLink.id1 > 0 ? contigPositionInScaffold[i][0].id : -(contigPositionInScaffold[i][0].id); //edited by ouchi (contigPositionInScaffold)
			else //added by ouchi (contigPositionInScaffold)
				overlapLink.id1 = 0; //added by ouchi (contigPositionInScaffold)

            long j = id2Index(overlapLink.id2);
            if (contigPositionInScaffold[j].size() == 1) //added by ouchi (contigPositionInScaffold)
	            overlapLink.id2 = overlapLink.id2 > 0 ? contigPositionInScaffold[j][0].id : -(contigPositionInScaffold[j][0].id); //edited by ouchi (contigPositionInScaffold)
			else //added by ouchi (contigPositionInScaffold)
				overlapLink.id2 = 0; //added by ouchi (contigPositionInScaffold)

			if (overlapLink.id1 * overlapLink.id2 == 0 || abs(overlapLink.id1) == abs(overlapLink.id2))
				continue;

			overlapLink.gap = -(getScaffoldOverlap(overlapLink.id1, overlapLink.id2));

			if (overlapLink.gap == -(this->minOverlap)) {
				if (abs(overlapLink.id1) > abs(overlapLink.id2)) {
					std::swap(overlapLink.id1, overlapLink.id2);
					overlapLink.id1 *= -1;
					overlapLink.id2 *= -1;
				}

				if (!(this->contigUnlinkSet.empty()) && this->contigUnlinkSet.find(i < j ? std::make_pair(i, j) : std::make_pair(j, i)) != this->contigUnlinkSet.end())
					continue;

				#pragma omp critical
				{
					graphLinkPool.push_back(overlapLink);
				}
			}
		}
    }
}

void PairedDBG::storeGraphLinkFromMappedPair(vector<GraphLinkWithFlag> &graphLinkPool, long numThread)
{
	if (!(this->mode & PAIRED_END_LINK_MODE) || allLibraryMT == NULL)
		return;

	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		unsigned numPairLink;
		unsigned numReadLink;

		rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				//added by ouchi
				long i = id2Index(forwardResult.id);
				long j = id2Index(reverseResult.id);
				long i2 = 0; long j2 = 0;
				/*if (contigPositionInScaffold[i].size() * contigPositionInScaffold[j].size() != 1) {
					platanus::Position scaffoldResultF, scaffoldResultR;
					long numLink = 0;
					for (long k = 0; k < contigPositionInScaffold[i].size(); ++k) {
						for (long l = 0; l < contigPositionInScaffold[j].size(); ++l) {
							scaffoldResultF.id = forwardResult.id > 0 ? contigPositionInScaffold[i][k].id : -(contigPositionInScaffold[i][k].id);
							scaffoldResultF.offset = contigPositionInScaffold[i][k].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1;
							scaffoldResultF.offset += node[id2Index(scaffoldResultF.id)].contig[contigPositionInScaffold[i][k].offset].start;

							scaffoldResultR.id = reverseResult.id > 0 ? contigPositionInScaffold[j][l].id : -(contigPositionInScaffold[j][l].id);
							scaffoldResultR.offset = contigPositionInScaffold[j][l].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1;
							scaffoldResultR.offset += node[id2Index(scaffoldResultR.id)].contig[contigPositionInScaffold[j][l].offset].start;

							if (abs(scaffoldResultF.id) == abs(scaffoldResultR.id)) {
								if (scaffoldResultF.id != -scaffoldResultR.id)
									continue;
								long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
								if (scaffoldResultF.id > 0) {
									if (scaffoldResultR.offset - scaffoldResultF.offset < 0 || scaffoldResultR.offset - scaffoldResultF.offset - averageInsSize > tolerence)
										continue;
								} else {
									if (scaffoldResultF.offset - scaffoldResultR.offset < 0 || scaffoldResultF.offset - scaffoldResultR.offset - averageInsSize > tolerence)
										continue;
								}
							} else {
								long gap = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
								if (scaffoldResultF.id > 0) {
									gap -= node[scaffoldResultF.id - 1].length - scaffoldResultF.offset;
								} else {
									gap -= scaffoldResultF.offset + 1;
								}
								if (scaffoldResultR.id > 0) {
									gap -= node[scaffoldResultR.id - 1].length -scaffoldResultR.offset;
								} else {
									gap -= scaffoldResultR.offset + 1;
								}

								if (-gap > std::min(tolerence, std::min(node[id2Index(scaffoldResultF.id)].length, node[id2Index(scaffoldResultR.id)].length)) + this->getScaffoldOverlap(scaffoldResultF.id, scaffoldResultR.id))
									continue;
							}

							++numLink;
							i2 = k; j2 = l;
							if (numLink > 1) {
								k = contigPositionInScaffold[i].size();
								break;
							}
						}
					}
					if (numLink != 1) continue;
				}
				// */

				//long i = id2Index(forwardResult.id); //deleted by ouchi
				if (contigPositionInScaffold[i].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i][i2].id : -(contigPositionInScaffold[i][i2].id); //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset = contigPositionInScaffold[i][i2].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i][i2].offset].start; //edited by ouchi (contigPositionInScaffold)

				//long j = id2Index(reverseResult.id); //deleted by ouci
				if (contigPositionInScaffold[j].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j][j2].id : -(contigPositionInScaffold[j][j2].id); //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset = contigPositionInScaffold[j][j2].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j][j2].offset].start; //edited by ouchi (contigPositionInScaffold)

				if (abs(forwardResult.id) == abs(reverseResult.id))
					continue;

				GraphLinkWithFlag graphLink;
				graphLink.overlapFlag = false;
				if (abs(forwardResult.id) < abs(reverseResult.id)) {
					graphLink.id1 = forwardResult.id;
					graphLink.offset1 = contigPositionInScaffold[i][i2].offset; //edited by ouchi (contigPositionInScaffold)
					graphLink.id2 = -(reverseResult.id);
					graphLink.offset2 = contigPositionInScaffold[j][j2].offset; //edited by ouchi (contigPositionInScaffold)
				}
				else {
					graphLink.id1 = reverseResult.id;
					graphLink.offset1 = contigPositionInScaffold[j][j2].offset; //edited by ouchi (contigPositionInScaffold)
					graphLink.id2 = -(forwardResult.id);
					graphLink.offset2 = contigPositionInScaffold[i][i2].offset; //edited by ouchi (contigPositionInScaffold)
				}

				std::pair<int, int> redundancyCheckKey = std::make_pair(graphLink.id1, graphLink.id2);
				if (redundancyCheckSet.find(redundancyCheckKey) != redundancyCheckSet.end())
					continue;

				redundancyCheckSet.insert(redundancyCheckKey);


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
				if (node[id2Index(forwardResult.id)].length < cutoffLength) continue;
				graphLink.gap -= gapMinus;


				if (reverseResult.id > 0) {
	//				if (node[reverseResult.id - 1].length < cutoffLength || abs(forwardResult.id) == abs(reverseResult.id)) continue;
					if (abs(forwardResult.id) == abs(reverseResult.id)) continue;
					graphLink.gap -= node[reverseResult.id-1].length - reverseResult.offset;
				} else {
	//				if (node[(-1 * reverseResult.id) - 1].length < cutoffLength || abs(forwardResult.id) == abs(reverseResult.id)) continue;
					if (abs(forwardResult.id) == abs(reverseResult.id)) continue;
					graphLink.gap -= reverseResult.offset + 1;
				}

				if (-graphLink.gap > std::min(tolerence, std::min(node[id2Index(graphLink.id1)].length, node[id2Index(graphLink.id2)].length)) + this->getScaffoldOverlap(graphLink.id1, graphLink.id2))
					continue;

				if (!(this->contigUnlinkSet.empty()) && this->contigUnlinkSet.find(i < j ? std::make_pair(i, j) : std::make_pair(j, i)) != this->contigUnlinkSet.end())
					continue;

				# pragma omp critical (push)
				{
					graphLinkPool.push_back(graphLink);
				}
			}


			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				long i = id2Index(forwardResult.id);
				if (contigPositionInScaffold[i].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i][0].id : -(contigPositionInScaffold[i][0].id); //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset = contigPositionInScaffold[i][0].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				long j = id2Index(reverseResult.id);
				if (contigPositionInScaffold[j].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j][0].id : -(contigPositionInScaffold[j][0].id); //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset = contigPositionInScaffold[j][0].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				if (abs(forwardResult.id) == abs(reverseResult.id))
					continue;

				GraphLinkWithFlag graphLink;
				graphLink.overlapFlag = false;

				long forwardLeft;
				long forwardRight;
				if (forwardResult.id > 0) {
					forwardLeft = -(forwardResult.offset);
					forwardRight = node[forwardResult.id - 1].length - forwardResult.offset - 1;
				}
				else {
					forwardLeft = -(node[-(forwardResult.id) - 1].length - forwardResult.offset - 1);
					forwardRight = forwardResult.offset;
				}

				long reverseLeft;
				long reverseRight;
				if (reverseResult.id > 0) {
					reverseLeft = -(reverseResult.offset);
					reverseRight = node[reverseResult.id - 1].length - reverseResult.offset - 1;
				}
				else {
					reverseLeft = -(node[-(reverseResult.id) - 1].length - reverseResult.offset - 1);
					reverseRight = reverseResult.offset;
				}

				if (forwardLeft <= reverseLeft) {
					if (forwardRight > reverseRight)
						continue;
					graphLink.gap = -(forwardRight - reverseLeft + 1);
				}
				else {
					if (reverseRight > forwardRight)
						continue;
					graphLink.gap = -(reverseRight - forwardLeft + 1);
				}

				if (abs(forwardResult.id) < abs(reverseResult.id)) {
					if (forwardRight < reverseRight) {
						graphLink.id1 = forwardResult.id;
						graphLink.offset1 = contigPositionInScaffold[i][0].offset; //edited by ouchi (contigPositionInScaffold)
						graphLink.id2 = reverseResult.id;
						graphLink.offset2 = contigPositionInScaffold[j][0].offset; //edited by ouchi (contigPositionInScaffold)
					}
					else {
						graphLink.id1 = -(forwardResult.id);
						graphLink.offset1 = contigPositionInScaffold[i][0].offset; //edited by ouchi (contigPositionInScaffold)
						graphLink.id2 = -(reverseResult.id);
						graphLink.offset2 = contigPositionInScaffold[j][0].offset; //edited by ouchi (contigPositionInScaffold)
					}
				}
				else {
					if (forwardRight < reverseRight) {
						graphLink.id1 = -(reverseResult.id);
						graphLink.offset1 = contigPositionInScaffold[j][0].offset; //edited by ouchi (contigPositionInScaffold)
						graphLink.id2 = -(forwardResult.id);
						graphLink.offset2 = contigPositionInScaffold[i][0].offset; //edited by ouchi (contigPositionInScaffold)
					}
					else {
						graphLink.id1 = reverseResult.id;
						graphLink.offset1 = contigPositionInScaffold[j][0].offset; //edited by ouchi (contigPositionInScaffold)
						graphLink.id2 = forwardResult.id;
						graphLink.offset2 = contigPositionInScaffold[i][0].offset; //edited by ouchi (contigPositionInScaffold)
					}
				}

				std::pair<int, int> redundancyCheckKey = std::make_pair(graphLink.id1, graphLink.id2);
				if (redundancyCheckSet.find(redundancyCheckKey) != redundancyCheckSet.end())
					continue;

				redundancyCheckSet.insert(redundancyCheckKey);


/*
				// check contigs aren't overlapeed too large area
				if (-(graphLink.gap) > this->minOverlap) continue;
				# pragma omp critical (push)
				{
					graphLinkPool.push_back(graphLink);
				}
*/
			}
		}
	}
}

void PairedDBG::storeGraphLinkFromMappedLongRead(vector<GraphLinkWithFlag> &graphLinkPool, long numThread)
{
	if (!(this->mode & LONG_READ_LINK_MODE) || longReadLibraryMT == NULL)
		return;

	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numResult;
		platanus::Position mapResult;
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		long score;

		rewind((*longReadLibraryMT)[threadID].mappedReadFP);
		while (fread(&numResult, sizeof(unsigned), 1, (*longReadLibraryMT)[threadID].mappedReadFP)) {

			vector<MapInfoForGraph> mapInfoBuffer;
			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&mapResult, sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fseek((*longReadLibraryMT)[threadID].mappedReadFP, 2 * sizeof(int), SEEK_CUR);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);

				long contigIndex = id2Index(mapResult.id);
				if (contigPositionInScaffold[contigIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				mapResult.id = mapResult.id > 0 ? contigPositionInScaffold[contigIndex][0].id : -(contigPositionInScaffold[contigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				mapResult.offset = contigPositionInScaffold[contigIndex][0].id > 0 ? mapResult.offset : contig[contigIndex].length - mapResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				mapResult.offset += node[id2Index(mapResult.id)].contig[contigPositionInScaffold[contigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				mapInfoBuffer.emplace_back(mapResult, contigIndex, score);
			}

			if (mapInfoBuffer.empty())
				continue;

			std::stable_sort(mapInfoBuffer.begin(), mapInfoBuffer.end(), MapInfoForGraph::PositionIDLessScoreGreater());

			vector<MapInfoForGraph> mergedMapResultBuffer;
			auto maxScoreItr = mapInfoBuffer.begin();
			while (maxScoreItr != mapInfoBuffer.end()) {
				auto boundItr = std::upper_bound(maxScoreItr, mapInfoBuffer.end(), *(maxScoreItr), MapInfoForGraph::PositionIDLess());

				score = 0;
				for (auto itr = maxScoreItr; itr != boundItr; ++itr)
					score += itr->score;

				mergedMapResultBuffer.emplace_back(maxScoreItr->position, maxScoreItr->contigIndex, score);

				maxScoreItr = boundItr;
			}

			for (auto mapItr1 = mergedMapResultBuffer.begin(); mapItr1 != mergedMapResultBuffer.end() - 1; ++mapItr1) {
				for (auto mapItr2 = mapItr1 + 1; mapItr2 != mergedMapResultBuffer.end(); ++mapItr2) {
					if (abs(mapItr1->position.id) == abs(mapItr2->position.id))
						continue;

					GraphLinkWithFlag graphLink;
					graphLink.overlapFlag = false;

					long forwardLeft;
					long forwardRight;
					if (mapItr1->position.id > 0) {
						forwardLeft = -(mapItr1->position.offset);
						forwardRight = node[mapItr1->position.id - 1].length - mapItr1->position.offset - 1;
					}
					else {
						forwardLeft = -(node[-(mapItr1->position.id) - 1].length - mapItr1->position.offset - 1);
						forwardRight = mapItr1->position.offset;
					}

					long reverseLeft;
					long reverseRight;
					if (mapItr2->position.id > 0) {
						reverseLeft = -(mapItr2->position.offset);
						reverseRight = node[mapItr2->position.id - 1].length - mapItr2->position.offset - 1;
					}
					else {
						reverseLeft = -(node[-(mapItr2->position.id) - 1].length - mapItr2->position.offset - 1);
						reverseRight = mapItr2->position.offset;
					}

					if (forwardLeft <= reverseLeft) {
						if (forwardRight > reverseRight)
							continue;
						graphLink.gap = -(forwardRight - reverseLeft + 1);
					}
					else {
						if (reverseRight > forwardRight)
							continue;
						graphLink.gap = -(reverseRight - forwardLeft + 1);
					}

					if (abs(mapItr1->position.id) < abs(mapItr2->position.id)) {
						if (forwardRight < reverseRight) {
							graphLink.id1 = mapItr1->position.id;
							graphLink.offset1 = contigPositionInScaffold[mapItr1->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
							graphLink.id2 = mapItr2->position.id;
							graphLink.offset2 = contigPositionInScaffold[mapItr2->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
						}
						else {
							graphLink.id1 = -(mapItr1->position.id);
							graphLink.offset1 = contigPositionInScaffold[mapItr1->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
							graphLink.id2 = -(mapItr2->position.id);
							graphLink.offset2 = contigPositionInScaffold[mapItr2->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
						}
					}
					else {
						if (forwardRight < reverseRight) {
							graphLink.id1 = -(mapItr2->position.id);
							graphLink.offset1 = contigPositionInScaffold[mapItr2->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
							graphLink.id2 = -(mapItr1->position.id);
							graphLink.offset2 = contigPositionInScaffold[mapItr1->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
						}
						else {
							graphLink.id1 = mapItr2->position.id;
							graphLink.offset1 = contigPositionInScaffold[mapItr2->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
							graphLink.id2 = mapItr1->position.id;
							graphLink.offset2 = contigPositionInScaffold[mapItr1->contigIndex][0].offset; //edited by ouchi (contigPositionInScaffold)
						}
					}

					if (-(graphLink.gap) > this->minOverlap)
						continue;

					graphLink.score = mapItr1->score + mapItr2->score;

					if (!(this->contigUnlinkSet.empty()) && this->contigUnlinkSet.find(mapItr1->contigIndex < mapItr2->contigIndex ? std::make_pair(mapItr1->contigIndex, mapItr2->contigIndex) : std::make_pair(mapItr2->contigIndex, mapItr1->contigIndex)) != this->contigUnlinkSet.end())
						continue;

					# pragma omp critical (push)
					{
						graphLinkPool.push_back(graphLink);
					}
				}
			}
		}
	}

}

void PairedDBG::calcLink(const long libraryIndex, const long linkThreshold, const long numThread)
{
    if (graphLinkFP[libraryIndex] != NULL)
        fclose(graphLinkFP[libraryIndex]);
    graphLinkFP[libraryIndex] = platanus::makeTemporaryFile();

//    cerr << "linking scaffolds (MIN_LINK = " << minLink << ")" << endl;
    cerr << "linking scaffolds..." << endl;

    vector<GraphLinkWithFlag> graphLinkPool;

	if (this->mode & OVERLAP_MODE)
		storeGraphLinkFromOverlap(graphLinkPool);

	if (this->mode & PAIRED_END_LINK_MODE)
		storeGraphLinkFromMappedPair(graphLinkPool, numThread);

	if (this->mode & LONG_READ_LINK_MODE)
		storeGraphLinkFromMappedLongRead(graphLinkPool, numThread);

    long totalLink = graphLinkPool.size();

	if (totalLink == 0)
		return;

    cerr << "sorting links in contigID order..." << endl;
    std::stable_sort(graphLinkPool.begin(), graphLinkPool.end());
    graphLinkPool.emplace_back();
    graphLinkPool[totalLink].id1 = 0;

    cerr << "estimating contig distances..." << endl;

    vector<GraphLinkWithFlagPoolIndex> indexes(1, 0);
	++indexes.back().numLink;
	if (graphLinkPool[0].overlapFlag)
		indexes.back().overlapFlag = true;
    for (long idx = 1; idx < totalLink; ++idx) {
        if (graphLinkPool[idx - 1].id1 != graphLinkPool[idx].id1 || graphLinkPool[idx - 1].id2 != graphLinkPool[idx].id2) {
            if (!(indexes.back().overlapFlag) && indexes.back().numLink < linkThreshold) {
                indexes.pop_back();
            }
            indexes.emplace_back(idx);
        }

        ++indexes.back().numLink;
		if (graphLinkPool[idx].overlapFlag) {
			indexes.back().overlapFlag = true;
			indexes.back().gap = graphLinkPool[idx].gap;
		}
    }
    if (!(indexes.back().overlapFlag) && indexes.back().numLink < linkThreshold) {
        indexes.pop_back();
    }
    sort(indexes.begin(), indexes.end(), GraphLinkWithFlagPoolIndexGreater());

	if (numThread == 1) {
		for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
			calcLinkAndWriteGraphLinkWithFlagFile(graphLinkPool, indexes[idxIndex], libraryIndex, false);
		}
	}
	else {
		// pragma omp parallel for schedule(dynamic) num_threads(numThread)
		for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
			calcLinkAndWriteGraphLinkWithFlagFile(graphLinkPool, indexes[idxIndex], libraryIndex, true);
		}
	}
}

void PairedDBG::calcLinkAndWriteGraphLinkWithFlagFile(const std::vector<GraphLinkWithFlag>& graphLinkPool, const GraphLinkWithFlagPoolIndex& index, const long libraryIndex, const bool multiThreadFlag)
{
    unsigned long indexBegin = index.index;
    long numLink = index.numLink;
	long score = 0;

    GraphLinkWithFlag graphLink = graphLinkPool[indexBegin];
	graphLink.overlapFlag = index.overlapFlag;
    vector<long> breakdown1(node[id2Index(graphLink.id1)].numContig, 0);
    vector<long> breakdown2(node[id2Index(graphLink.id2)].numContig, 0);
    for (long idxLink = indexBegin, endLink = idxLink + numLink; idxLink < endLink; ++idxLink) {
        ++(breakdown1[graphLinkPool[idxLink].offset1]);
        ++(breakdown2[graphLinkPool[idxLink].offset2]);
    }

	if (index.overlapFlag) {
		graphLink.gap = index.gap;
	}
	else {
		vector<GraphLink> targetLinkPool(numLink);
		for (long i = 0; i < numLink; ++i) {
			targetLinkPool[i] = static_cast<GraphLink>(graphLinkPool[i + indexBegin]);
			score += graphLinkPool[i + indexBegin].score;
		}
//		graphLink.gap = estimateGapSize(targetLinkPool, 0, numLink);
		graphLink.gap = estimateGapSizeAverage(targetLinkPool, 0, numLink);
	}
    if (graphLink.gap == LONG_MIN) return;

    long breakdownSize = breakdown1.size() + breakdown2.size();
    std::unique_ptr<long[]> tmp(new long[breakdownSize]());
    std::copy(breakdown1.begin(), breakdown1.end(), tmp.get());
    std::copy(breakdown2.begin(), breakdown2.end(), &(tmp[breakdown1.size()]));

	if (multiThreadFlag) {
		#pragma omp critical (write_link)
		{
			fwrite(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex]);
			fwrite(&score, sizeof(long), 1, graphLinkFP[libraryIndex]);
			fwrite(&graphLink, sizeof(GraphLinkWithFlag), 1, graphLinkFP[libraryIndex]);
			fwrite(tmp.get(), sizeof(long), breakdownSize, graphLinkFP[libraryIndex]);
		}
	}
	else {
		fwrite(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex]);
		fwrite(&score, sizeof(long), 1, graphLinkFP[libraryIndex]);
		fwrite(&graphLink, sizeof(GraphLinkWithFlag), 1, graphLinkFP[libraryIndex]);
		fwrite(tmp.get(), sizeof(long), breakdownSize, graphLinkFP[libraryIndex]);
	}
}

//added by ouchi
void PairedDBG::makeMPLink(const long numThread, const long pairLibraryBegin, const long pairLibraryEnd)
{
	if (allLibraryMT == NULL)
		return;

    for (int i = 0; i < numNode; ++i) {
		node[i].MPLinks.clear();
		node[i].MPLinks.resize((*allLibraryMT).size(), {});
	}

	for (unsigned libraryID = pairLibraryBegin; libraryID < pairLibraryEnd; ++libraryID) {
		# pragma omp parallel for schedule(static, 1)
		for (long threadID = 0; threadID < numThread; ++threadID) {
			platanus::Position forwardResult;
			platanus::Position reverseResult;
			unsigned numPairLink;
			unsigned numReadLink;

			rewind((*allLibraryMT)[libraryID][threadID].mappedFP);
			while (fread(&numPairLink, sizeof(unsigned), 1, (*allLibraryMT)[libraryID][threadID].mappedFP)) {
				fread(&numReadLink, sizeof(unsigned), 1, (*allLibraryMT)[libraryID][threadID].mappedFP);
				std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

				for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
					fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[libraryID][threadID].mappedFP);
					fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[libraryID][threadID].mappedFP);

					long i = id2Index(forwardResult.id);
					long j = id2Index(reverseResult.id);
					long i2 = 0; long j2 = 0;

					if (contigPositionInScaffold[i].size() != 1) continue;
					forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i][i2].id : -(contigPositionInScaffold[i][i2].id);
					forwardResult.offset = contigPositionInScaffold[i][i2].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1;
					forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i][i2].offset].start;

					if (contigPositionInScaffold[j].size() != 1) continue;
					reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j][j2].id : -(contigPositionInScaffold[j][j2].id);
					reverseResult.offset = contigPositionInScaffold[j][j2].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1;
					reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j][j2].offset].start;

					if (abs(forwardResult.id) == abs(reverseResult.id))
						continue;

					std::pair<int, int> redundancyCheckKey;
					if (abs(forwardResult.id) < abs(reverseResult.id)) {
						redundancyCheckKey = std::make_pair(forwardResult.id, -(reverseResult.id));
					} else {
						redundancyCheckKey = std::make_pair(reverseResult.id, -(forwardResult.id));
					}
					if (redundancyCheckSet.find(redundancyCheckKey) != redundancyCheckSet.end())
						continue;
					redundancyCheckSet.insert(redundancyCheckKey);

					/*long gap = (*allLibraryMT)[libraryID][0].getAverageInsSize();
					if (forwardResult.id > 0) {
						gap -= (node[forwardResult.id - 1].length - forwardResult.offset);
					} else {
						gap -= (forwardResult.offset + 1);
					}
					if (reverseResult.id > 0) {
						gap -= (node[reverseResult.id - 1].length - reverseResult.offset);
					} else {
						gap -= (reverseResult.offset + 1);
					}
					if (-gap > std::min(tolerence, std::min(node[id2Index(forwardResult.id)].length, node[id2Index(reverseResult.id)].length)) + this->getScaffoldOverlap(forwardResult.id, -(reverseResult.id)))
						continue;*/

					MPLink mpLink1, mpLink2;
					mpLink1.direction = sign(forwardResult.id);
					if (forwardResult.id > 0) {
						mpLink1.end = -(reverseResult.id);
					} else {
						mpLink1.end = reverseResult.id;
					}
					mpLink1.offset1 = forwardResult.offset;
					mpLink1.offset2 = reverseResult.offset;
					mpLink2.direction = sign(reverseResult.id);
					if (forwardResult.id > 0) {
						mpLink2.end = -(forwardResult.id);
					} else {
						mpLink2.end = forwardResult.id;
					}
					mpLink2.offset1 = reverseResult.offset;
					mpLink2.offset2 = forwardResult.offset;
					# pragma omp critical (push)
					{
						node[id2Index(forwardResult.id)].MPLinks[libraryID].push_back(mpLink1);
						node[id2Index(reverseResult.id)].MPLinks[libraryID].push_back(mpLink2);
					}

				}

				for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
					fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[libraryID][threadID].mappedFP);
					fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[libraryID][threadID].mappedFP);
				}

			}
		}
	}

}
//

//added by ouchi
void PairedDBG::makeHiCLink(const long numThread, const long HiCLinkMode)
{
    if (HiCLibraryMT == NULL)
        return;

    for (int i = 0; i < numNode; ++i) {
		node[i].HiCLinks.clear();
	}

	# pragma omp parallel for schedule(static, 1)
	for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;

		rewind((*HiCLibraryMT)[threadID].mappedFP);
		while (fread(&forwardResult, sizeof(platanus::Position), 1, (*HiCLibraryMT)[threadID].mappedFP)) {
			fread(&reverseResult, sizeof(platanus::Position), 1, (*HiCLibraryMT)[threadID].mappedFP);

			if (HiCLinkMode == 0) {
				if (contigPositionInScaffold[id2Index(forwardResult.id)].size() != 1)
					continue;
				if (contigPositionInScaffold[id2Index(reverseResult.id)].size() != 1)
					continue;
			}

			if (contigPositionInScaffold[id2Index(forwardResult.id)].size() > 2)
				continue;
			if (contigPositionInScaffold[id2Index(reverseResult.id)].size() > 2)
				continue;

			long k1 = 0;
			if (contigPositionInScaffold[id2Index(forwardResult.id)].size() == 2) {
				long k = 0;
				long consensusID1 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(forwardResult.id)][k].id)][0].second.id;
				long tmpOffset = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(forwardResult.id)][k].id)][0].second.offset;
				char h1 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(forwardResult.id)][k].id)][0].first;
				long offset1 = consensusnode[id2Index(consensusID1)].node[h1][tmpOffset].start;
				if (consensusID1 > 0)
					offset1 += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[id2Index(forwardResult.id)][k].offset].start;
				else
					offset1 += node[id2Index(forwardResult.id)].length - node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[id2Index(forwardResult.id)][k].offset].end;
				if (consensusID1 * contigPositionInScaffold[id2Index(forwardResult.id)][k].id > 0)
					offset1 += forwardResult.offset;
				else
					offset1 += contig[id2Index(forwardResult.id)].length - forwardResult.offset - 1;
				k = 1;
				long consensusID2 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(forwardResult.id)][k].id)][0].second.id;
				tmpOffset = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(forwardResult.id)][k].id)][0].second.offset;
				char h2 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(forwardResult.id)][k].id)][0].first;
				long offset2 = consensusnode[id2Index(consensusID2)].node[h2][tmpOffset].start;
				if (consensusID2 > 0)
					offset2 += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[id2Index(forwardResult.id)][k].offset].start;
				else
					offset2 += node[id2Index(forwardResult.id)].length - node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[id2Index(forwardResult.id)][k].offset].end;
				if (consensusID2 * contigPositionInScaffold[id2Index(forwardResult.id)][k].id > 0)
					offset2 += forwardResult.offset;
				else
					offset2 += contig[id2Index(forwardResult.id)].length - forwardResult.offset - 1;
				if (!(abs(consensusID1) == abs(consensusID2) && h1 != h2 && abs(offset1 - offset2) < tolerence))
					continue;
				if (HiCLinkMode == 2)
					k1 = 1;
			}

			long k2 = 0;
			if (contigPositionInScaffold[id2Index(reverseResult.id)].size() == 2) {
				long k = 0;
				long consensusID1 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(reverseResult.id)][k].id)][0].second.id;
				long tmpOffset = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(reverseResult.id)][k].id)][0].second.offset;
				char h1 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(reverseResult.id)][k].id)][0].first;
				long offset1 = consensusnode[id2Index(consensusID1)].node[h1][tmpOffset].start;
				if (consensusID1 > 0)
					offset1 += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[id2Index(reverseResult.id)][k].offset].start;
				else
					offset1 += node[id2Index(reverseResult.id)].length - node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[id2Index(reverseResult.id)][k].offset].end;
				if (consensusID1 * contigPositionInScaffold[id2Index(reverseResult.id)][k].id > 0)
					offset1 += reverseResult.offset;
				else
					offset1 += contig[id2Index(reverseResult.id)].length - reverseResult.offset - 1;
				k = 1;
				long consensusID2 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(reverseResult.id)][k].id)][0].second.id;
				tmpOffset = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(reverseResult.id)][k].id)][0].second.offset;
				char h2 = nodePositionInConsensus[id2Index(contigPositionInScaffold[id2Index(reverseResult.id)][k].id)][0].first;
				long offset2 = consensusnode[id2Index(consensusID2)].node[h2][tmpOffset].start;
				if (consensusID2 > 0)
					offset2 += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[id2Index(reverseResult.id)][k].offset].start;
				else
					offset2 += node[id2Index(reverseResult.id)].length - node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[id2Index(reverseResult.id)][k].offset].end;
				if (consensusID2 * contigPositionInScaffold[id2Index(reverseResult.id)][k].id > 0)
					offset2 += reverseResult.offset;
				else
					offset2 += contig[id2Index(reverseResult.id)].length - reverseResult.offset - 1;
				if (!(abs(consensusID1) == abs(consensusID2) && h1 != h2 && abs(offset1 - offset2) < tolerence))
					continue;
				if (HiCLinkMode == 2)
					k2 = 1;
			}

			for (long i2 = 0; i2 <= k1; ++i2) {
				for (long j2 = 0; j2 <= k2; ++j2) {


					long i = id2Index(forwardResult.id);
					//if (contigPositionInScaffold[i].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
					forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[i][i2].id : -(contigPositionInScaffold[i][i2].id); //edited by ouchi (contigPositionInScaffold)
					forwardResult.offset = contigPositionInScaffold[i][i2].id > 0 ? forwardResult.offset : contig[i].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
					forwardResult.offset += node[id2Index(forwardResult.id)].contig[contigPositionInScaffold[i][i2].offset].start; //edited by ouchi (contigPositionInScaffold)

					long j = id2Index(reverseResult.id);
					//if (contigPositionInScaffold[j].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
					reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[j][j2].id : -(contigPositionInScaffold[j][j2].id); //edited by ouchi (contigPositionInScaffold)
					reverseResult.offset = contigPositionInScaffold[j][j2].id > 0 ? reverseResult.offset : contig[j].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
					reverseResult.offset += node[id2Index(reverseResult.id)].contig[contigPositionInScaffold[j][j2].offset].start; //edited by ouchi (contigPositionInScaffold)

					HiCLink hicLink1, hicLink2;
					hicLink1.end = reverseResult.id;
					hicLink1.offset1 = forwardResult.offset;
					hicLink1.offset2 = reverseResult.offset;
					hicLink2.end = forwardResult.id;
					hicLink2.offset1 = reverseResult.offset;
					hicLink2.offset2 = forwardResult.offset;
					# pragma omp critical (push)
					{
						node[id2Index(forwardResult.id)].HiCLinks.push_back(hicLink1);
						node[id2Index(reverseResult.id)].HiCLinks.push_back(hicLink2);
					}


				}
			}
		}
	}

}
//


void PairedDBG::makeGraph(const long numThread)
{
	unsigned numPairedEndLibrary = 0;
	if (allLibraryMT != NULL)
		numPairedEndLibrary = (*allLibraryMT).size();

	unsigned numLongReadLibrary = 0;
	if (longReadLibraryMT != NULL)
		numLongReadLibrary = 1;

	if (graphLinkFP.empty())
		graphLinkFP.resize(std::max(numPairedEndLibrary + numLongReadLibrary, 1u), NULL);

    calcLink(this->targetLibraryIndex, this->minLink, numThread);

    cerr << "constructing paired de Bruijn graph" << endl;

    for (int i = 0; i < numNode; ++i) {
		node[i].edge.clear();
		node[i].numEdge = 0;
	}

    rewind(graphLinkFP[this->targetLibraryIndex]);
    long numLink;
    long score;
    while (fread(&numLink, sizeof(long), 1, graphLinkFP[this->targetLibraryIndex])) {
		fread(&score, sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        GraphLinkWithFlag link;
        fread(&link, sizeof(GraphLinkWithFlag), 1, graphLinkFP[this->targetLibraryIndex]);

        long i = id2Index(link.id1);
        long j = id2Index(link.id2);
        long m = node[i].numEdge;
        long n = node[j].numEdge;

		node[i].edge.resize(m + 1);
		node[j].edge.resize(n + 1);


		if (link.overlapFlag) {
			node[i].edge[m].state = DBG_OVERLAP;
			node[j].edge[n].state = DBG_OVERLAP;
			node[i].edge[m].numLink = numLink - 1;
			node[j].edge[n].numLink = numLink - 1;
		}
		else {
			node[i].edge[m].state = 0x0;
			node[j].edge[n].state = 0x0;
			node[i].edge[m].numLink = numLink;
			node[j].edge[n].numLink = numLink;
		}

		node[i].edge[m].score = score;
		node[j].edge[n].score = score;

		node[i].edge[m].breakdown.resize(node[i].numContig, 0);
        for (unsigned int k = 0; k < node[i].numContig; ++k) {
            fread(&(node[i].edge[m].breakdown[k]), sizeof(long), 1, graphLinkFP[this->targetLibraryIndex]);
        }
		node[j].edge[n].breakdown.resize(node[j].numContig, 0);
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

	setOppositeBubbleNodeIDAndStateForEachNode();

	if (getMode() & LENGTH_CUTOFF_MODE) {
		if (getMode() & BUBBLE_AWARE_MODE)
			deleteEdgeFromShortNodeKeepingBubble(this->cutoffLength);
		else
			deleteEdgeFromShortNode(this->cutoffLength);
	}

	if (getMode() & PREVIOUS_DIVISION_AWARE_MODE)
		deleteEdgeFromDifferentPreviousParent(numThread);

//	if (getMode() & SECONDARY_BUBBLE_REMOVAL_MODE)
//		deleteEdgeFromSecondaryBubble();
}

void PairedDBG::makeGraphAllLibraries(const long numThread)
{
	unsigned numPairedEndLibrary = 0;
	if (allLibraryMT != NULL)
		numPairedEndLibrary = this->targetLibraryIndex + 1;
//		numPairedEndLibrary = (*allLibraryMT).size();

	unsigned numLongReadLibrary = 0;
	if (longReadLibraryMT != NULL)
		numLongReadLibrary = 1;

	if (graphLinkFP.empty())
		graphLinkFP.resize(numPairedEndLibrary + numLongReadLibrary, NULL);

	unsigned currentMode = getMode();
	for (unsigned libraryIndex = 0; libraryIndex < numPairedEndLibrary; ++libraryIndex) {
		this->setTargetLibraryIndex(libraryIndex);

//		if (getMode() & LENGTH_CUTOFF_MODE)
//			this->setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * (*allLibraryMT)[libraryIndex][0].getSDInsSize(), 0.1 * (*allLibraryMT)[libraryIndex][0].getAverageInsSize()));
//		else
//			this->setCutoffLength(0);

//		this->setTolerence(this->cutoffLength / 2);

		if (libraryIndex > 0 && (getMode() & OVERLAP_MODE))
			unsetMode(OVERLAP_MODE);

		calcLink(libraryIndex, 1, numThread);
	}

	for (unsigned libraryIndex = numPairedEndLibrary; libraryIndex < numPairedEndLibrary + numLongReadLibrary; ++libraryIndex) {
		this->setTargetLibraryIndex(libraryIndex);
//		this->setTolerence(2 * (this->minOverlap + 1));
//		if (getMode() & LENGTH_CUTOFF_MODE)
//			this->setCutoffLength(this->cutoffLengthFactor * (*longReadLibraryMT)[0].getAverageInsSize());
//		else
//			this->setCutoffLength(0);

		if (libraryIndex > 0 && (getMode() & OVERLAP_MODE))
			unsetMode(OVERLAP_MODE);

		unsetMode(PAIRED_END_LINK_MODE);
		calcLink(libraryIndex, 1, numThread);
	}

	setMode(currentMode);

    cerr << "constructing paired de Bruijn graph using all libraries simultaneously" << endl;

    for (int i = 0; i < numNode; ++i) {
		node[i].edge.clear();
		node[i].numEdge = 0;
	}

	for (unsigned libraryIndex = 0; libraryIndex < numPairedEndLibrary + numLongReadLibrary; ++libraryIndex) {
		rewind(graphLinkFP[libraryIndex]);
		long numLink;
		long score;
		while (fread(&numLink, sizeof(long), 1, graphLinkFP[libraryIndex])) {
			fread(&score, sizeof(long), 1, graphLinkFP[libraryIndex]);
			GraphLinkWithFlag link;
			fread(&link, sizeof(GraphLinkWithFlag), 1, graphLinkFP[libraryIndex]);

			long i = id2Index(link.id1);
			long m = 0 ;
			for (; m < node[i].numEdge; ++m) {
				if (node[i].edge[m].end == sign(link.id1) * link.id2 && node[i].edge[m].direction == sign(link.id1))
					break;
			}
			if (m == node[i].numEdge) {
				node[i].edge.resize(m + 1);
				++(node[i].numEdge);
			}

			long j = id2Index(link.id2);
			long n = 0;
			for (; n < node[j].numEdge; ++n) {
				if (node[j].edge[n].end == sign(link.id2) * link.id1 && node[j].edge[n].direction == -sign(link.id2))
					break;
			}
			if (n == node[j].numEdge) {
				node[j].edge.resize(n + 1);
				++(node[j].numEdge);
			}


			if (link.overlapFlag) {
				node[i].edge[m].state = DBG_OVERLAP;
				node[j].edge[n].state = DBG_OVERLAP;
				node[i].edge[m].numLink += (numLink - 1);
				node[j].edge[n].numLink += (numLink - 1);
			}
			else {
				node[i].edge[m].numLink += numLink;
				node[j].edge[n].numLink += numLink;
			}

			node[i].edge[m].score = score;
			node[j].edge[n].score = score;

			long breakdownNum = 0;
			node[i].edge[m].breakdown.resize(node[i].numContig, 0);
			for (unsigned int k = 0; k < node[i].numContig; ++k) {
				fread(&breakdownNum, sizeof(long), 1, graphLinkFP[libraryIndex]);
				node[i].edge[m].breakdown[k] += breakdownNum;
			}
			node[j].edge[n].breakdown.resize(node[j].numContig, 0);
			for (unsigned int k = 0; k < node[j].numContig; ++k) {
				fread(&breakdownNum, sizeof(long), 1, graphLinkFP[libraryIndex]);
				node[j].edge[n].breakdown[k] += breakdownNum;
			}

			if (!(link.overlapFlag)) {
				node[i].edge[m].length += link.gap * numLink;
				node[j].edge[n].length += link.gap * numLink;
			}
			else {
				node[i].edge[m].length = link.gap;
				node[j].edge[n].length = link.gap;
			}

			node[i].edge[m].direction = sign(link.id1);
			node[j].edge[n].direction = -sign(link.id2);

			if (link.id1 * link.id2 > 0) {
				node[i].edge[m].end = j + 1;
				node[j].edge[n].end = i + 1;
			}
			else {
				node[i].edge[m].end = -(j + 1);
				node[j].edge[n].end = -(i + 1);
			}
		}
		fclose(graphLinkFP[libraryIndex]);
		graphLinkFP[libraryIndex] = NULL;
	}

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &tmp = node[nodeIndex];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);

		for (long edgeIndex = 0; edgeIndex < node[nodeIndex].numEdge; ++edgeIndex) {
			if (!(node[nodeIndex].edge[edgeIndex].state & DBG_OVERLAP))
				node[nodeIndex].edge[edgeIndex].length /= node[nodeIndex].edge[edgeIndex].numLink;
		}
    }

	deleteThinEdgeCostantKeepingOverlap(this->minLink);

	if (this->mode & PAIRED_END_LINK_MODE && allLibraryMT != NULL)
		this->setTargetLibraryIndex(numPairedEndLibrary - 1);

	setOppositeBubbleNodeIDAndStateForEachNode();

	if (getMode() & LENGTH_CUTOFF_MODE) {
		if (getMode() & BUBBLE_AWARE_MODE)
			deleteEdgeFromShortNodeKeepingBubble(this->cutoffLength);
		else
			deleteEdgeFromShortNode(this->cutoffLength);
	}

//	if (getMode() & SECONDARY_BUBBLE_REMOVAL_MODE)
//		deleteEdgeFromSecondaryBubble();
}

void PairedDBG::resetGraph()
{
    this->numNode = this->numContig;
    node.clear();
    node.resize(this->numNode);

    contigPositionInScaffold.clear();
    contigPositionInScaffold.resize(this->numContig);

    numBubble.clear();
    numBubble.resize(this->numContig, 0);

    for (long i = 0; i < numNode; ++i) {
        node[i].state = 0;
        node[i].numContig = 1;
        node[i].contig.clear();
        node[i].contig.resize(1);
        node[i].contig[0].id = i + 1;
        node[i].length = node[i].contig[0].end = this->contig[i].length;
        contigPositionInScaffold[i].resize(1); //added by ouchi (contigPositionInScaffold)
        contigPositionInScaffold[i][0].id = i + 1; //edited by ouchi (contigPositionInScaffold)
        contigPositionInScaffold[i][0].offset = 0; //edited by ouchi (contigPositionInScaffold)
    }

	for (unsigned long i = 0; i < this->contigState.size(); ++i) {
		contigState[i] = 0;
		contigBubbleInfo[i] = PairedDBG::ContigBubbleInfo();
	}
}

void PairedDBG::getOverlappedBubbleNodeIndex(vector<long> &bubbleNodeIndex)
{
	vector<long> nodeIDBuffer;
	vector<vector<long> > sinkNodeIDBuffer(2);
	vector<char> bubbleNodeFlag(this->numNode, 0);

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (sourceDirection == 1) {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, 1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, -1, sinkNodeIDBuffer[i]);
				}
				else {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, -1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, 1, sinkNodeIDBuffer[i]);
				}

				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0] *= sign(nodeIDBuffer[i]);
			}

			if (i == 2 && sinkNodeIDBuffer[0][0] == sinkNodeIDBuffer[1][0]) {
				bubbleNodeFlag[id2Index(nodeIDBuffer[0])] = 1;
				bubbleNodeFlag[id2Index(nodeIDBuffer[1])] = 1;
			}
		}
    }

	bubbleNodeIndex.clear();
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (bubbleNodeFlag[nodeIndex]) {
			bubbleNodeIndex.push_back(nodeIndex);
		}
	}
}

void PairedDBG::getOverlappedBubbleNodePairID(vector<pair<long, long> > &bubbleNodePairID)
{
//	const double HOMO_COVERAGE_THRESHOLD = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<long> nodeIDBuffer;
	vector<vector<long> > sinkNodeIDBuffer(2);
	vector<char> bubbleNodeFlag(this->numNode, 0);

	bubbleNodePairID.clear();

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
//		if (calcNodeCoverage(node[sourceNodeIndex]) > HOMO_COVERAGE_THRESHOLD)
//			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				getOverlappedNode(id2Index(nodeIDBuffer[i]), sign(nodeIDBuffer[i]) * sourceDirection, sinkNodeIDBuffer[i]);
				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0] *= sign(nodeIDBuffer[i]);
			}

			if (i != 2 || sinkNodeIDBuffer[0][0] != sinkNodeIDBuffer[1][0])
				continue;

			vector<long> nodeIDBufferForCheck;
			getOverlappedNode(id2Index(sinkNodeIDBuffer[0][0]), -(sign(sinkNodeIDBuffer[0][0]) * sourceDirection), nodeIDBufferForCheck);
			if (nodeIDBufferForCheck.size() != 2)
				continue;

			bubbleNodePairID.emplace_back(nodeIDBuffer[0], nodeIDBuffer[1]);
		}
    }
}

void PairedDBG::getOverlappedForkNodePairID(vector<pair<long, long> > &bubbleNodePairID, vector<char> &directionBuffer)
{
	const double HOMO_COVERAGE_THRESHOLD = 2.0 * HETERO_FORK_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<long> nodeIDBuffer;
	vector<long> nodeIDBufferForCheck;
	vector<char> bubbleNodeFlag(this->numNode, 0);

	bubbleNodePairID.clear();
	directionBuffer.clear();

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		if (calcNodeCoverage(node[sourceNodeIndex]) > HOMO_COVERAGE_THRESHOLD)
			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				char checkDirection = -1;
				for (; checkDirection <= 1; checkDirection += 2) {
					getOverlappedNode(id2Index(nodeIDBuffer[i]), sign(nodeIDBuffer[i]) * sourceDirection * checkDirection, nodeIDBufferForCheck);
					if (nodeIDBufferForCheck.size() != 1)
						break;
				}
				if (checkDirection <= 1)
					break;
			}
			if (i != 2)
				continue;

			bubbleNodePairID.emplace_back(nodeIDBuffer[0], nodeIDBuffer[1]);
			directionBuffer.push_back(sourceDirection);
		}
    }
}

void PairedDBG::getGappedBubbleNodePairID(vector<pair<long, long> > &bubbleNodePairID)
{
//	const double HOMO_COVERAGE_THRESHOLD = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<NodeIDWithGap> nodeIDBuffer;
	vector<vector<NodeIDWithGap> > sinkNodeIDBuffer(2);

	bubbleNodePairID.clear();

    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
//		if (calcNodeCoverage(node[sourceNodeIndex]) > HOMO_COVERAGE_THRESHOLD)
//			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getUniqueConflictingNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2 || abs(nodeIDBuffer[0].id) == abs(nodeIDBuffer[1].id))
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				getLinkedNode(id2Index(nodeIDBuffer[i].id), sign(nodeIDBuffer[i].id) * (-sourceDirection), sinkNodeIDBuffer[i]);
				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				getLinkedNode(id2Index(nodeIDBuffer[i].id), sign(nodeIDBuffer[i].id) * sourceDirection, sinkNodeIDBuffer[i]);
				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0].id *= sign(nodeIDBuffer[i].id);
			}

//			if (i == 2 && sinkNodeIDBuffer[0][0].id == sinkNodeIDBuffer[1][0].id && sourceNodeIndex != id2Index(sinkNodeIDBuffer[0][0].id)
//				&& calcNodeCoverage(node[id2Index(sinkNodeIDBuffer[0][0].id)]) <= HOMO_COVERAGE_THRESHOLD) {
			if (i == 2 && sinkNodeIDBuffer[0][0].id == sinkNodeIDBuffer[1][0].id && sourceNodeIndex != id2Index(sinkNodeIDBuffer[0][0].id)) {
				bubbleNodePairID.emplace_back(nodeIDBuffer[0].id, nodeIDBuffer[1].id);
			}
		}
    }
}

void PairedDBG::getGappedConflictingNodePairID(vector<pair<long, long> > &bubbleNodePairID, double heteroNodeCoverage)
{
	const double MAX_HETERO_COVERAGE = HETERO_COVERAGE_THRESHOLD_FACTOR * heteroNodeCoverage;
	const double MAX_HOMO_COVERAGE = 3.0 * heteroNodeCoverage;
	vector<NodeIDWithGap> nodeIDBuffer;

	bubbleNodePairID.clear();
    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		if (calcNodeCoverage(this->node[sourceNodeIndex]) > MAX_HOMO_COVERAGE)
			continue;

		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getUniqueConflictingNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() == 2 && calcNodeCoverage(this->node[id2Index(nodeIDBuffer[0].id)]) <= MAX_HETERO_COVERAGE && calcNodeCoverage(this->node[id2Index(nodeIDBuffer[1].id)]) <= MAX_HETERO_COVERAGE) {
				bubbleNodePairID.emplace_back(nodeIDBuffer[0].id, nodeIDBuffer[1].id);
			}
		}
    }
}

void PairedDBG::getOverlappedNode(const long sourceNodeIndex, const char targetDirection, vector<long> &nodeIDBuffer)
{
	nodeIDBuffer.clear();
	GraphNode &node = this->node[sourceNodeIndex];
	for (long edgeIndex = 0; edgeIndex < node.numEdge; ++edgeIndex) {
		GraphEdge &edge = node.edge[edgeIndex];
		if (edge.direction == targetDirection && ((!(this->mode & NON_DBG_OVERLAP_MODE) && (edge.state & DBG_OVERLAP)) || ((this->mode & NON_DBG_OVERLAP_MODE) && edge.length < 0)))
			nodeIDBuffer.push_back(edge.end);
	}
}

long PairedDBG::getNumEdgeOneDirection(const GraphNode &targetNode, const char targetDirection)
{
	long n = 0;
	for (auto itr = targetNode.edge.begin(); itr != targetNode.edge.end(); ++itr) {
		if (itr->direction == targetDirection)
			++n;
	}
	return n;
}

void PairedDBG::markHeteroNode(const double maxHeteroCoverageFactor)
{
	const double maxHeteroCoverage = maxHeteroCoverageFactor * this->heteroCoverage;

	for (auto it = this->node.begin(); it != this->node.end(); ++it) {
		if (calcNodeCoverage(*it) <= maxHeteroCoverage)
			it->state |= DBG_HETERO;
	}
}

void PairedDBG::markBubbleHeteroNode(const vector<long> &candidateNodeIndex, const double maxHeteroCoverageFactor)
{
	const double maxHeteroCoverage = maxHeteroCoverageFactor * this->heteroCoverage;

	for (auto it = candidateNodeIndex.begin(); it != candidateNodeIndex.end(); ++it) {
		if (calcNodeCoverage(this->node[*it]) <= maxHeteroCoverage)
			this->node[*it].state |= DBG_HETERO;
	}
}

void PairedDBG::calculateHeteroCoverage(const vector<long> &bubbleNodeIndex)
{
	const long MIN_NUM_BUBBLE = 10000;
	const double TRUNCATION_FACTOR = 2.0;

	vector<char> bubbleNodeFlag(this->node.size(), 0);

	long numBubble = 0;
	for (auto itr = bubbleNodeIndex.begin(); itr != bubbleNodeIndex.end(); ++itr) {
		bubbleNodeFlag[*itr] = 1;
		++numBubble;
	}

	vector<std::pair<int, int> > buffer(2);
	for (unsigned long i = 0; i < this->node.size(); ++i) {
		if (this->node[i].length <= this->contigMaxK)
			continue;

		if (bubbleNodeFlag[i])
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i]) + 0.5), static_cast<int>(this->node[i].length)));
		else if (numBubble < MIN_NUM_BUBBLE)
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i])/2.0 + 0.5), static_cast<int>(this->node[i].length)));
	}

	if (!buffer.empty()) {
		long sum = 0;
		long totalLength = 0;
		double weightedMean  = 0.0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			sum +=  itr->first * itr->second;
			totalLength += itr->second;
		}
		weightedMean  = std::round(static_cast<double>(sum) / totalLength);

		sum = totalLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			if (weightedMean / TRUNCATION_FACTOR  <= itr->first && itr->first <= weightedMean * TRUNCATION_FACTOR) {
				sum +=  itr->first * itr->second;
				totalLength += itr->second;
			}
		}
		if (totalLength > 0)
			this->heteroCoverage = std::round(static_cast<double>(sum) / totalLength);
		else
			this->heteroCoverage = weightedMean;
	}
	else {
		this->heteroCoverage = 1.0;
	}
	this->averageCoverage = 2.0 * this->heteroCoverage;

/*
	double weightedMedian = 1.0;
	if (!buffer.empty()) {
		std::sort(buffer.begin(), buffer.end(), platanus::PairFirstLess());

		long totalLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr)
			totalLength += itr->second;

		long sumLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			sumLength += itr->second;
			if (sumLength >= totalLength / 2) {
				weightedMedian = itr->first;
				break;
			}
		}
	}
	this->heteroCoverage = weightedMedian;
	this->averageCoverage = 2.0 * this->heteroCoverage;
*/

	cerr << "ESTIMATED_HETERO_COVERAGE = "<< this->heteroCoverage << endl;
}

void PairedDBG::calculateHeteroAndAverageCoverageUnphase()
{
	const double TRUNCATION_FACTOR = 2.0;

	vector<std::pair<int, int> > buffer(2);
	unsigned long i;
	for (i = 0; i < this->node.size() - this->numInputBubbleContig; ++i) {
		if (this->node[i].length > this->contigMaxK)
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i])/2.0 + 0.5), static_cast<int>(this->node[i].length)));
	}
	for (; i < this->node.size(); ++i) {
		if (this->node[i].length > this->contigMaxK)
			buffer.push_back(std::make_pair(static_cast<int>(calcNodeCoverage(this->node[i]) + 0.5), static_cast<int>(this->node[i].length)));
	}

	if (!buffer.empty()) {
		long sum = 0;
		long totalLength = 0;
		double weightedMean  = 0.0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			sum +=  itr->first * itr->second;
			totalLength += itr->second;
		}
		weightedMean  = std::round(static_cast<double>(sum) / totalLength);

		sum = totalLength = 0;
		for (auto itr = buffer.begin(); itr != buffer.end(); ++itr) {
			if (weightedMean / TRUNCATION_FACTOR  <= itr->first && itr->first <= weightedMean * TRUNCATION_FACTOR) {
				sum +=  itr->first * itr->second;
				totalLength += itr->second;
			}
		}
		if (totalLength > 0)
			this->heteroCoverage = std::round(static_cast<double>(sum) / totalLength);
		else
			this->heteroCoverage = weightedMean;
	}
	else {
		this->heteroCoverage = 1.0;
	}
	this->averageCoverage = 2.0 * this->heteroCoverage;

	cerr << "ESTIMATED_HETERO_COVERAGE = "<< this->heteroCoverage << endl;
}

void PairedDBG::extractDBGBubbleInformation()
{
	const double MAX_HETERO_BUBBLE_COVERAGE_FACTOR = 2.0;
	vector<long> bubbleNodeIndex;

	getOverlappedBubbleNodeIndex(bubbleNodeIndex);
	if (this->heteroCoverage <= 0.0) {
		calculateHeteroCoverage(bubbleNodeIndex);
		this->averageCoverage = 2.0 * this->heteroCoverage;
	}
	markBubbleHeteroNode(bubbleNodeIndex, MAX_HETERO_BUBBLE_COVERAGE_FACTOR);
}

long PairedDBG::crushSimpleDBGBubble()
{
    const double coverageThreshold = heteroCoverage * 3.0;

    vector<long> ids;
	vector<long> nodeIDBuffer;
	vector<vector<long> > sinkNodeIDBuffer(2);

	long numCrush = 0;
    for (long sourceNodeIndex = 0; sourceNodeIndex < numNode; ++sourceNodeIndex) {
		for (char sourceDirection = -1; sourceDirection <= 1; sourceDirection += 2) {
			getOverlappedNode(sourceNodeIndex, sourceDirection, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (sourceDirection == 1) {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, 1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, -1, sinkNodeIDBuffer[i]);
				}
				else {
					if (nodeIDBuffer[i] >= 0)
						getOverlappedNode(nodeIDBuffer[i] - 1, -1, sinkNodeIDBuffer[i]);
					else
						getOverlappedNode(-(nodeIDBuffer[i]) - 1, 1, sinkNodeIDBuffer[i]);
				}

				if (sinkNodeIDBuffer[i].size() != 1)
					break;

				sinkNodeIDBuffer[i][0] *= sign(nodeIDBuffer[i]);
			}

			long index1 = id2Index(nodeIDBuffer[0]);
			long index2 = id2Index(nodeIDBuffer[1]);
			GraphNode *node1 = &(node[index1]);
			GraphNode *node2 = &(node[index2]);
			double coverage1 = this->calcNodeCoverage(*node1);
			double coverage2 = this->calcNodeCoverage(*node2);

			if (i != 2 ||
				sinkNodeIDBuffer[0][0] != sinkNodeIDBuffer[1][0] ||
				(node1->state & SC_DEL) ||
				(node2->state & SC_DEL) ||
				coverage1 + coverage2 > coverageThreshold) {

				continue;
			}

			getOverlappedNode(id2Index(sinkNodeIDBuffer[0][0]), -(sourceDirection * sign(sinkNodeIDBuffer[0][0])), nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			if (coverage1 > coverage2) {
				std::swap(node1, node2);
				std::swap(index1, index2);
			}

			node1->state |= SC_DEL;
			for (long n = 0; n < node1->numEdge; ++n) {
				ids.push_back(index1 + 1);
				ids.push_back(node1->edge[n].end);
			}

			writeNodeSeq(*node1, bubbleFP);
			writeNodeSeq(*node2, bubbleOpositeFP);

			++numCrush;
		}
    }

    this->deleteEdges(ids);

	return numCrush;
}

void PairedDBG::deleteNonOverlapHomoEdge()
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing non-hetero-informative edges..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
			if (node1.edge[edgeID].state & DBG_OVERLAP)
				continue;

            GraphNode &node2 = this->node[id2Index(node1.edge[edgeID].end)];
			if ((node1.state & DBG_HETERO) && (node2.state & DBG_HETERO))
				continue;

			ids.push_back(nodeID + 1);
			ids.push_back(node1.edge[edgeID].end);
			++numDelete;
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::joinUnambiguousNodePair(const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & SC_DEL)
				continue;

			for (long i = 0; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getOverlappedNode(initialNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 1)
					continue;

				long altNodeID = nodeIDBuffer[0];

				getOverlappedNode(id2Index(altNodeID), -(initialDirection * sign(altNodeID)), nodeIDBuffer);
				if (nodeIDBuffer.size() != 1)
					continue;

				//added by ouchi for test
				//if (abs(calcNodeCoverage(this->node[initialNodeIndex]) - calcNodeCoverage(this->node[id2Index(altNodeID)])) > this->heteroCoverage)
				//	continue;
				//

				pathBuffer[t].emplace_back(initialNodeIndex, 2);
				if (i == 0) {
					pathBuffer[t].back().nodeID[0] = altNodeID;
					pathBuffer[t].back().nodeID[1] = initialNodeIndex + 1;
				}
				else {
					pathBuffer[t].back().nodeID[0] = initialNodeIndex + 1;
					pathBuffer[t].back().nodeID[1] = altNodeID;
				}

				break;
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numJoin = remakeGraphAccordingToPath(mergedPathBuffer);

	return numJoin;
}

long PairedDBG::joinUnambiguousNodePath(const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		vector<char> visitFlag(this->numNode, 0);
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & SC_DEL)
				continue;

			visitFlag[initialNodeIndex] = 1;
			std::array<vector<long>, 2> bothSideNodeID;
			for (long i = 0; i < 2; ++i) {
				long currentNodeID = initialNodeIndex + 1;
				char initialDirection = 2*i - 1;
				while (1) {
					getOverlappedNode(id2Index(currentNodeID), initialDirection * sign(currentNodeID), nodeIDBuffer);
					if (nodeIDBuffer.size() != 1)
						break;

					long nextNodeID = nodeIDBuffer[0] * sign(currentNodeID);
					getOverlappedNode(id2Index(nextNodeID), -(initialDirection * sign(nextNodeID)), nodeIDBuffer);
					if (nodeIDBuffer.size() != 1)
						break;

					if (visitFlag[id2Index(nextNodeID)])
						break;

					visitFlag[id2Index(nextNodeID)] = 1;
					bothSideNodeID[i].push_back(nextNodeID);
					currentNodeID = nextNodeID;
				}
			}
			pathBuffer[t].emplace_back(initialNodeIndex, 0);
			GraphPath &path = pathBuffer[t].back();
			for (unsigned long i = 0; i < bothSideNodeID[0].size(); ++i)
				path.nodeID.push_back(bothSideNodeID[0][bothSideNodeID[0].size() - i - 1]);
			path.nodeID.push_back(initialNodeIndex + 1);
			for (unsigned long i = 0; i < bothSideNodeID[1].size(); ++i)
				path.nodeID.push_back(bothSideNodeID[1][i]);
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numJoin = remakeGraphAccordingToPath(mergedPathBuffer);

	return numJoin;
}

long PairedDBG::remakeGraphAccordingToPath(vector<GraphPath> &pathBuffer)
{
	long numNewNode = 0;
	long newContigPoolSize = 0;
	long numJoin = 0;

    FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); ++pathItr) {
		if ((node[id2Index(pathItr->nodeID.front())].state & SC_INC) || (node[id2Index(pathItr->nodeID.back())].state & SC_INC))
			continue;

		numJoin += pathItr->nodeID.size();
		newContigPoolSize += writeAndMarkOverlappedNodes(pathItr->nodeID, newContigFP);
		++numNewNode;
	}

	this->writeSingletonNode(numNewNode, newContigPoolSize, newContigFP);
    this->remake(numNewNode, newContigPoolSize, newContigFP);
    fclose(newContigFP);

	return numJoin;
}

long PairedDBG::remakeGraphAccordingToPathPair(vector<GraphPath> &pathBuffer)
{
	long numNewNode = 0;
	long newContigPoolSize = 0;
	long numJoin = 0;
	vector<long> nodeIDForWriting;
	vector<long> nodeIDForWriting2;

    FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); pathItr += 2) {
		nodeIDForWriting.clear();
		auto nodeIDItr = pathItr->nodeID.begin();
		for (; nodeIDItr != pathItr->nodeID.end(); ++nodeIDItr) {
			if (node[id2Index(*nodeIDItr)].state & SC_INC)
				break;
			else
				nodeIDForWriting.push_back(*nodeIDItr);
		}
		if (nodeIDItr != pathItr->nodeID.end())
			continue;

		nodeIDForWriting2.clear();
		auto pathItr2 = pathItr + 1;
		nodeIDItr = pathItr2->nodeID.begin();
		for (; nodeIDItr != pathItr2->nodeID.end(); ++nodeIDItr) {
			if (node[id2Index(*nodeIDItr)].state & SC_INC)
				break;
			else
				nodeIDForWriting2.push_back(*nodeIDItr);
		}
		if (nodeIDItr != pathItr2->nodeID.end())
			continue;

		numJoin += (pathItr->nodeID.size() + pathItr2->nodeID.size());
		newContigPoolSize += writeAndMarkOverlappedNodes(nodeIDForWriting, newContigFP);
		newContigPoolSize += writeAndMarkOverlappedNodes(nodeIDForWriting2, newContigFP);
		numNewNode += 2;
	}

	this->writeSingletonNode(numNewNode, newContigPoolSize, newContigFP);
    this->remake(numNewNode, newContigPoolSize, newContigFP);
    fclose(newContigFP);

	return numJoin;
}

long PairedDBG::solveSimpleCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;
//	const double COVERAGE_THRESHOLD = 2 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<vector<GraphPath> > pathBuffer(numThread);

	//added by ouchi for test
	long numCross = 0;
	long numSolvedCrossbyMP = 0;
	long numMissedbyHiC =0;
	long numCorrectbyHiC = 0;
	long numSolvedCrossbyHiC = 0;
	vector<std::array<long, 4> > contiglengthlist;
	vector<std::array<double, 4> > HiCScorelist;
	vector<std::array<double, 4> > MPScorelist;
	vector<std::array<int, 2> > CrossResultlist;
	//

    # pragma omp parallel for schedule(static, 1) reduction(+: numCross, numSolvedCrossbyMP, numMissedbyHiC, numSolvedCrossbyHiC) //edited by ouchi for test (reduction~)
	for (long t = 0; t < numThread; ++t) {
		vector<long> nodeIDBuffer;
		vector<long> nodeIDBufferForCheck;
		std::array<std::array<long, 2>, 2 > externalNodeID;

		for (long centerNodeIndex = t; centerNodeIndex < numNode; centerNodeIndex += numThread) {
			if (2 * COVERAGE_THRESHOLD < calcNodeCoverage(this->node[centerNodeIndex]) || this->node[centerNodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getOverlappedNode(centerNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 2)
					break;

				if ((this->mode & BUBBLE_AWARE_MODE) && !(checkNodeConvergenceBothSide(nodeIDBuffer[0], nodeIDBuffer[1])))
					break;

				long j = 0;
				for (; j < 2; ++j) {
					externalNodeID[i][j] = nodeIDBuffer[j];
					getOverlappedNode(id2Index(externalNodeID[i][j]),  -(initialDirection *sign(externalNodeID[i][j])), nodeIDBufferForCheck); //edited by ouchi , sign(externalNodeID[i][j],
					if (nodeIDBufferForCheck.size() >= 2)
						break;
				}
				if (j < 2)
					break;
			}
			if (i < 2)
				continue;

			//added by ouchi
			if (externalNodeID[0][0] == externalNodeID[1][0] || externalNodeID[0][0] == externalNodeID[1][1] ||
				externalNodeID[0][1] == externalNodeID[1][0] || externalNodeID[0][1] == externalNodeID[1][1])
				continue;
			//

//			if (COVERAGE_THRESHOLD < std::max(std::max(calcNodeCoverage(this->node[id2Index(externalNodeID[0][0])]), calcNodeCoverage(this->node[id2Index(externalNodeID[0][1])])), std::max(calcNodeCoverage(this->node[id2Index(externalNodeID[1][0])]), calcNodeCoverage(this->node[id2Index(externalNodeID[1][1])]))))
			if (COVERAGE_THRESHOLD < std::min(std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[0][0])]), calcNodeCoverage(this->node[id2Index(externalNodeID[0][1])])), std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[1][0])]), calcNodeCoverage(this->node[id2Index(externalNodeID[1][1])]))))
				continue;

			++numCross; //added by ouchi for test;
			std::array<double, 2> sumLinkForHaplotypeforHiC; //added by ouchi for test
			sumLinkForHaplotypeforHiC.fill(0); //added by ouchi for test

			std::array<long, 2> sumLinkForHaplotype;
			sumLinkForHaplotype.fill(0);
			if (resolutionMode == TAG) {
				if (node[centerNodeIndex].length > maxFragmentLengthOfTag)
					continue;

				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getCommonTagBetweenNodePair(externalNodeID[0][leftIndex], externalNodeID[1][rightIndex]);
					}
				}
			}
			if (resolutionMode == SCORE) {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getScoreFromIDPair(externalNodeID[0][leftIndex], externalNodeID[1][rightIndex]);
					}
				}
			}
			else {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getNumLinkFromIDPair(externalNodeID[0][leftIndex], externalNodeID[1][rightIndex]);
					}
				}
			}

			//deleted by ouchi
			/*if ((resolutionMode == LINK || resolutionMode == TAG)  && std::max(sumLinkForHaplotype[0], sumLinkForHaplotype[1]) < this->minLink)
				continue;

			int crossFlag;
			if (linkRateThreshold * sumLinkForHaplotype[0] >= sumLinkForHaplotype[1])
				crossFlag = 1;
			else if (linkRateThreshold * sumLinkForHaplotype[1] >= sumLinkForHaplotype[0])
				crossFlag = 0;
			else
				continue;*/
			//

//			++numSolvedCrossbyMP; //added by ouchi for test

			//added by ouchi for test
			std::array<double, 4> HiCScores;
			HiCScores.fill(0);
			std::array<double, 4> MPScores;
			MPScores.fill(0);
			std::array<double, 2> sumLinkForHaplotypeforMP;
			sumLinkForHaplotypeforMP.fill(0);
            std::array<std::array<NodeIDWithGap, 2>, 2> externalNodeIDWithGap;
            externalNodeIDWithGap[0][0].id = externalNodeID[0][0];
            externalNodeIDWithGap[0][0].gap = -this->minOverlap;
            externalNodeIDWithGap[0][1].id = externalNodeID[0][1];
            externalNodeIDWithGap[0][1].gap = -this->minOverlap;
            externalNodeIDWithGap[1][0].id = externalNodeID[1][0];
            externalNodeIDWithGap[1][0].gap = -this->minOverlap;
            externalNodeIDWithGap[1][1].id = externalNodeID[1][1];
            externalNodeIDWithGap[1][1].gap = -this->minOverlap;
			sumLinkForHaplotypeforMP[0] = sumLinkForHaplotype[0];
			sumLinkForHaplotypeforMP[1] = sumLinkForHaplotype[1];
			//getHeteroMPLinkNum(centerNodeIndex, externalNodeIDWithGap, sumLinkForHaplotypeforMP, MPScores); //added by ouchi for test
            int crossFlag = -1;
            if (linkRateThreshold * sumLinkForHaplotypeforMP[0] >= sumLinkForHaplotypeforMP[1])
                crossFlag = 1;
            else if (linkRateThreshold * sumLinkForHaplotypeforMP[1] >= sumLinkForHaplotypeforMP[0])
                crossFlag = 0;
			if ((resolutionMode == LINK || resolutionMode == TAG) && std::max(sumLinkForHaplotypeforMP[0], sumLinkForHaplotypeforMP[1]) < this->minLink)
				crossFlag = -1;
            getHeteroHiCLinkScore(centerNodeIndex, externalNodeIDWithGap, sumLinkForHaplotypeforHiC, HiCScores);
			/*long minleftLength = std::min(node[id2Index(externalNodeID[0][0])].length, node[id2Index(externalNodeID[0][1])].length);
			long minrightLength = std::min(node[id2Index(externalNodeID[1][0])].length, node[id2Index(externalNodeID[1][1])].length);
			for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
				for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
					sumLinkForHaplotypeforHiC[leftIndex == rightIndex] += this->getHiCLinkScoreBetweenNodePair(externalNodeID[0][leftIndex], sign(externalNodeID[0][leftIndex]), externalNodeID[1][rightIndex], minleftLength, minrightLength);
				}
			}*/
			int HiCcrossFlag = -1;
			if (linkRateThreshold * sumLinkForHaplotypeforHiC[0] > sumLinkForHaplotypeforHiC[1]) {
				HiCcrossFlag = 1;
			}
			else if (linkRateThreshold * sumLinkForHaplotypeforHiC[1] > sumLinkForHaplotypeforHiC[0]) {
				HiCcrossFlag = 0;
			}
			if (std::max(sumLinkForHaplotypeforHiC[0], sumLinkForHaplotypeforHiC[1]) <= 3)
				HiCcrossFlag = -1;
			if (crossFlag != -1)
				++numSolvedCrossbyMP;
			if (HiCcrossFlag != -1)
				++numSolvedCrossbyHiC;
			if (crossFlag != -1 && HiCcrossFlag != -1) {
				if (crossFlag != HiCcrossFlag)
					++numMissedbyHiC;
				else
					++numCorrectbyHiC;
			}

			# pragma omp critical (push)
			{
				contiglengthlist.push_back(std::array<long, 4>{node[id2Index(externalNodeID[0][0])].length, node[id2Index(externalNodeID[0][1])].length, node[id2Index(externalNodeID[1][0])].length, node[id2Index(externalNodeID[1][1])].length});
				HiCScorelist.push_back(HiCScores);
				MPScorelist.push_back(MPScores);
				CrossResultlist.push_back(std::array<int, 2>{crossFlag, HiCcrossFlag});
			}
			if (crossFlag == -1) {
				if (HiCcrossFlag == -1)
					continue;
				else
					crossFlag = HiCcrossFlag;
			} else {
				if (HiCcrossFlag != -1) {
					if (crossFlag != HiCcrossFlag)
						continue;
				}
			}
			//

			for (long j = 0; j < 2; ++j) {
				pathBuffer[t].emplace_back(centerNodeIndex, 3);
				pathBuffer[t].back().nodeID[0] = externalNodeID[0][j];
				pathBuffer[t].back().nodeID[1] = centerNodeIndex + 1;
				pathBuffer[t].back().nodeID[2] = externalNodeID[1][(j + crossFlag) % 2];
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numSolvedCross = remakeGraphAccordingToPathPair(mergedPathBuffer);

	//added by ouchi for test
	//std::cerr<< "numCross:" << numCross <<"\tnumSolvedbyMP:" << numSolvedCrossbyMP << "\tnumMissedbyHiC:" << numMissedbyHiC << "\tnumCorrectbyHiC:" << numCorrectbyHiC << "\tnumSolvedbyHiC:" << numSolvedCrossbyHiC << std::endl;
	/*std::ofstream ofs;
	static int iteration;
	ofs.open("solveCross" + std::to_string(iteration++));
	for (long i = 0; i < contiglengthlist.size(); ++i)
		ofs << contiglengthlist[i][0] <<"\t"<< contiglengthlist[i][1] <<"\t"<< contiglengthlist[i][2] <<"\t"<< contiglengthlist[i][3] <<"\t";
	ofs << std::endl;
	for (long i = 0; i < HiCScorelist.size(); ++i)
		ofs << HiCScorelist[i][0] <<"\t"<< HiCScorelist[i][1] <<"\t"<< HiCScorelist[i][2] <<"\t"<< HiCScorelist[i][3] <<"\t";
	ofs << std::endl;
	for (long i = 0; i < MPScorelist.size(); ++i)
		ofs << MPScorelist[i][0] <<"\t"<< MPScorelist[i][1] <<"\t"<< MPScorelist[i][2] <<"\t"<< MPScorelist[i][3] <<"\t";
	ofs << std::endl;
	for (long i = 0; i < CrossResultlist.size(); ++i)
		ofs << CrossResultlist[i][0] <<"\t"<< CrossResultlist[i][1] <<"\t";
	ofs << std::endl;
	ofs.close();*/
	//

	return numSolvedCross;
}

long PairedDBG::solveSimpleGappedCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;
//	const double COVERAGE_THRESHOLD = 2 * HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	vector<vector<GraphPathGapped> > pathBuffer(numThread);

	//added by ouchi for test
	long numCross = 0;
	long numSolvedCrossbyMP = 0;
	long numMissedbyHiC =0;
	long numCorrectbyHiC =0;
	long numSolvedCrossbyHiC = 0;
	vector<std::array<long, 4> > contiglengthlist;
	vector<std::array<double, 4> > HiCScorelist;
	vector<std::array<double, 4> > MPScorelist;
	vector<std::array<int, 2> > CrossResultlist;
	//

    # pragma omp parallel for schedule(static, 1) reduction(+: numCross, numSolvedCrossbyMP, numMissedbyHiC, numSolvedCrossbyHiC) //edited by ouchi for test
	for (long t = 0; t < numThread; ++t) {
		vector<NodeIDWithGap> nodeIDBuffer;
		vector<NodeIDWithGap> nodeIDBufferForCheck;
		std::array<std::array<NodeIDWithGap, 2>, 2 > externalNodeID;
		std::array<std::array<long, 2>, 2 > externalLinkedNodeID; //added by ouchi

		for (long centerNodeIndex = t; centerNodeIndex < numNode; centerNodeIndex += numThread) {
			if (2 * COVERAGE_THRESHOLD < calcNodeCoverage(this->node[centerNodeIndex]) || this->node[centerNodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getUniqueConflictingNode(centerNodeIndex, initialDirection, nodeIDBuffer);
				if (nodeIDBuffer.size() != 2) {
					break;
				}
/*
				if (!(this->node[id2Index(nodeIDBuffer[0].id)].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)) ||
					!(this->node[id2Index(nodeIDBuffer[1].id)].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)) ||
					sign(nodeIDBuffer[0].id) * this->node[id2Index(nodeIDBuffer[0].id)].oppositeBubbleNodeID != nodeIDBuffer[1].id) {
					break;
				}
*/
				long j = 0;
				for (; j < 2; ++j) {
					externalNodeID[i][j] = nodeIDBuffer[j];
					getLinkedNode(id2Index(externalNodeID[i][j].id), -(initialDirection *sign(externalNodeID[i][j].id)), nodeIDBufferForCheck); //edited by ouchi , sign(externalNodeID[i][j].id),
					if (nodeIDBufferForCheck.size() >= 3)
						break;
					//added by ouchi
					externalLinkedNodeID[i][j] = 0;
					for (long k = 0; k < nodeIDBufferForCheck.size(); ++k) {
						if (nodeIDBufferForCheck[k].id != (centerNodeIndex+1) * sign(externalNodeID[i][j].id))
							externalLinkedNodeID[i][j] = nodeIDBufferForCheck[k].id;
					}
					//
				}
				if (j < 2)
					break;
			}
			if (i < 2)
				continue;

			//added by ouchi
			for (i = 0; i < 2; ++i) {
				if (externalLinkedNodeID[i][0] != 0 && abs(externalLinkedNodeID[i][0]) != abs(externalNodeID[(i+1)%2][0].id) && abs(externalLinkedNodeID[i][0]) != abs(externalNodeID[(i+1)%2][1].id))
					break;
				if (externalLinkedNodeID[i][1] != 0 && abs(externalLinkedNodeID[i][1]) != abs(externalNodeID[(i+1)%2][0].id) && abs(externalLinkedNodeID[i][1]) != abs(externalNodeID[(i+1)%2][1].id))
					break;
			}
			if (i < 2)
				continue;
			if (externalNodeID[0][0].id == externalNodeID[1][0].id || externalNodeID[0][0].id == externalNodeID[1][1].id ||
				externalNodeID[0][1].id == externalNodeID[1][0].id || externalNodeID[0][1].id == externalNodeID[1][1].id)
				continue;
			//

//			if (COVERAGE_THRESHOLD < std::max(std::max(calcNodeCoverage(this->node[id2Index(externalNodeID[0][0].id)]), calcNodeCoverage(this->node[id2Index(externalNodeID[0][1].id)])), std::max(calcNodeCoverage(this->node[id2Index(externalNodeID[1][0].id)]), calcNodeCoverage(this->node[id2Index(externalNodeID[1][1].id)]))))
			if (COVERAGE_THRESHOLD < std::min(std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[0][0].id)]), calcNodeCoverage(this->node[id2Index(externalNodeID[0][1].id)])), std::min(calcNodeCoverage(this->node[id2Index(externalNodeID[1][0].id)]), calcNodeCoverage(this->node[id2Index(externalNodeID[1][1].id)]))))
				continue;

			++numCross; //added by ouchi for test;
			std::array<double, 2> sumLinkForHaplotypeforHiC; //added by ouchi for test
			sumLinkForHaplotypeforHiC.fill(0); //added by ouchi for test

			std::array<long, 2> sumLinkForHaplotype;
			sumLinkForHaplotype.fill(0);
			if (resolutionMode == TAG) {
				if (node[centerNodeIndex].length > maxFragmentLengthOfTag)
					continue;

				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getCommonTagBetweenNodePair(externalNodeID[0][leftIndex].id, externalNodeID[1][rightIndex].id);
					}
				}
			}
			if (resolutionMode == SCORE) {
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getScoreFromIDPair(externalNodeID[0][leftIndex].id, externalNodeID[1][rightIndex].id);
					}
				}
			}
			if (resolutionMode == LINK || resolutionMode == TAG) { //edited by ouchi (else)
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getNumLinkFromIDPair(externalNodeID[0][leftIndex].id, externalNodeID[1][rightIndex].id);
					}
				}
			}
			//added by ouchi
			if (resolutionMode == HIC) {
				long minleftLength = std::min(node[id2Index(externalNodeID[0][0].id)].length, node[id2Index(externalNodeID[0][1].id)].length);
				long minrightLength = std::min(node[id2Index(externalNodeID[1][0].id)].length, node[id2Index(externalNodeID[1][1].id)].length);
				for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
					for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
						sumLinkForHaplotype[leftIndex == rightIndex] += this->getHiCLinkScoreBetweenNodePair(externalNodeID[0][leftIndex].id, sign(externalNodeID[0][leftIndex].id), externalNodeID[1][rightIndex].id, minleftLength, minrightLength);
					}
				}
			}
			//

/*
if (tagFlag) {
	cerr << "ga " << sumLinkForHaplotype[0] << ", " << sumLinkForHaplotype[1] << endl;
	if (std::max(sumLinkForHaplotype[0], sumLinkForHaplotype[1]) >= this->minLink)
		cerr << "s\n" << endl;
	else
		cerr << "u\n" << endl;
}
*/

			//deleted by ouchi
			/*if ((resolutionMode == LINK || resolutionMode == TAG) && std::max(sumLinkForHaplotype[0], sumLinkForHaplotype[1]) < this->minLink)
				continue;

			int crossFlag;
			if (linkRateThreshold * sumLinkForHaplotype[0] >= sumLinkForHaplotype[1])
				crossFlag = 1;
			else if (linkRateThreshold * sumLinkForHaplotype[1] >= sumLinkForHaplotype[0])
				crossFlag = 0;
			else
				continue;*/
			//

			//++numSolvedCrossbyMP; //added by ouchi for test

			//added by ouchi for test
			std::array<double, 4> HiCScores;
			HiCScores.fill(0);
			std::array<double, 4> MPScores;
			MPScores.fill(0);
			std::array<double, 2> sumLinkForHaplotypeforMP;
			sumLinkForHaplotypeforMP.fill(0);
			//getHeteroMPLinkNum(centerNodeIndex, externalNodeID, sumLinkForHaplotypeforMP, MPScores); //added by ouchi for test
			sumLinkForHaplotypeforMP[0] = sumLinkForHaplotype[0];
			sumLinkForHaplotypeforMP[1] = sumLinkForHaplotype[1];
            int crossFlag = -1;
            if (linkRateThreshold * sumLinkForHaplotypeforMP[0] >= sumLinkForHaplotypeforMP[1])
                crossFlag = 1;
            else if (linkRateThreshold * sumLinkForHaplotypeforMP[1] >= sumLinkForHaplotypeforMP[0])
                crossFlag = 0;
			if ((resolutionMode == LINK || resolutionMode == TAG) && std::max(sumLinkForHaplotypeforMP[0], sumLinkForHaplotypeforMP[1]) < this->minLink)
				crossFlag = -1;
            getHeteroHiCLinkScore(centerNodeIndex, externalNodeID, sumLinkForHaplotypeforHiC, HiCScores);
			/*long minleftLength = std::min(node[id2Index(externalNodeID[0][0].id)].length, node[id2Index(externalNodeID[0][1].id)].length);
			long minrightLength = std::min(node[id2Index(externalNodeID[1][0].id)].length, node[id2Index(externalNodeID[1][1].id)].length);
			for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
				for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
					sumLinkForHaplotypeforHiC[leftIndex == rightIndex] += this->getHiCLinkScoreBetweenNodePair(externalNodeID[0][leftIndex].id, sign(externalNodeID[0][leftIndex].id), externalNodeID[1][rightIndex].id, minleftLength, minrightLength);
				}
			}*/
			int HiCcrossFlag = -1;
			if (linkRateThreshold * sumLinkForHaplotypeforHiC[0] > sumLinkForHaplotypeforHiC[1]) {
				HiCcrossFlag = 1;
			}
			else if (linkRateThreshold * sumLinkForHaplotypeforHiC[1] > sumLinkForHaplotypeforHiC[0]) {
				HiCcrossFlag = 0;
			}
			if (std::max(sumLinkForHaplotypeforHiC[0], sumLinkForHaplotypeforHiC[1]) <= 3)
				HiCcrossFlag = -1;
			if (crossFlag != -1)
				++numSolvedCrossbyMP;
			if (HiCcrossFlag != -1)
				++numSolvedCrossbyHiC;
			if (crossFlag != -1 && HiCcrossFlag != -1) {
				if (crossFlag != HiCcrossFlag)
					++numMissedbyHiC;
				else
					++numCorrectbyHiC;
			}

			# pragma omp critical (push)
			{
				contiglengthlist.push_back(std::array<long, 4>{node[id2Index(externalNodeID[0][0].id)].length, node[id2Index(externalNodeID[0][1].id)].length, node[id2Index(externalNodeID[1][0].id)].length, node[id2Index(externalNodeID[1][1].id)].length});
				HiCScorelist.push_back(HiCScores);
				MPScorelist.push_back(MPScores);
				CrossResultlist.push_back(std::array<int, 2>{crossFlag, HiCcrossFlag});
			}
			if (crossFlag == -1) {
				if (HiCcrossFlag == -1)
					continue;
				else
					crossFlag = HiCcrossFlag;
			} else {
				if (HiCcrossFlag != -1) {
					if (crossFlag != HiCcrossFlag)
						continue;
				}
			}
			//


			//added by ouchi
			/*if (resolutionMode == HIC) {
				for (i = 0; i < 2; ++i) {
					if (!(checkHiCLinkBetweenNodePair(externalNodeID[0][i].id, sign(externalNodeID[0][i].id), externalNodeID[1][(i+crossFlag)%2].id)))
						break;
				}
				if (i < 2)
					continue;
			}*/
			//

			for (long j = 0; j < 2; ++j) {
				pathBuffer[t].emplace_back(centerNodeIndex, 3);
				pathBuffer[t].back().node[0].id = externalNodeID[0][j].id;
				pathBuffer[t].back().node[0].gap = 0;
				pathBuffer[t].back().node[1].id = centerNodeIndex + 1;
				pathBuffer[t].back().node[1].gap = externalNodeID[0][j].gap;
				pathBuffer[t].back().node[2] = externalNodeID[1][(j + crossFlag) % 2];
			}
		}
	}

	vector<GraphPathGapped> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathGappedSelfIDLess());

	long numSolvedCross = remakeGraphAccordingToGappedPathPair(mergedPathBuffer);

	//added by ouchi for test
	//std::cerr<< "numCross:" << numCross <<"\tnumSolvedbyMP:" << numSolvedCrossbyMP << "\tnumMissedbyHiC:" << numMissedbyHiC << "\tnumCorrectbyHiC:" << numCorrectbyHiC << "\tnumSolvedbyHiC:" << numSolvedCrossbyHiC << std::endl;
	/*std::ofstream ofs;
	static int iteration;
	ofs.open("solveGappedCross" + std::to_string(iteration++));
	for (long i = 0; i < contiglengthlist.size(); ++i)
		ofs << contiglengthlist[i][0] <<"\t"<< contiglengthlist[i][1] <<"\t"<< contiglengthlist[i][2] <<"\t"<< contiglengthlist[i][3] <<"\t";
	ofs << std::endl;
	for (long i = 0; i < HiCScorelist.size(); ++i)
		ofs << HiCScorelist[i][0] <<"\t"<< HiCScorelist[i][1] <<"\t"<< HiCScorelist[i][2] <<"\t"<< HiCScorelist[i][3] <<"\t";
	ofs << std::endl;
	for (long i = 0; i < MPScorelist.size(); ++i)
		ofs << MPScorelist[i][0] <<"\t"<< MPScorelist[i][1] <<"\t"<< MPScorelist[i][2] <<"\t"<< MPScorelist[i][3] <<"\t";
	ofs << std::endl;
	for (long i = 0; i < CrossResultlist.size(); ++i)
		ofs << CrossResultlist[i][0] <<"\t"<< CrossResultlist[i][1] <<"\t";
	ofs << std::endl;
	ofs.close();*/
	//

	return numSolvedCross;
}

long PairedDBG::writeAndMarkOverlappedNodes(const vector<long> &nodeID, FILE *storeFP)
{
	long numContig = 0;

	for (auto it = nodeID.begin(); it != nodeID.end(); ++it)
		numContig += node[id2Index(*it)].numContig;

	fwrite(&numContig, sizeof(long), 1, storeFP);

	long start = 0;
	for (auto it = nodeID.begin(); it != nodeID.end(); ++it) {
		long index;
		if (*it > 0) {
			index = *it - 1;
			node[index].state |= SC_INC;
			for (long i = 0; i < node[index].numContig; ++i) {
				ScaffoldPart scaffoldPartForWriting(node[index].contig[i].id, start + node[index].contig[i].start, start + node[index].contig[i].end);
				fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
			}
		}
		else {
			index = -(*it) - 1;
			node[index].state |= SC_INC;
			for (long i = node[index].numContig - 1; i >= 0; --i) {
				ScaffoldPart scaffoldPartForWriting(-(node[index].contig[i].id), start + node[index].length - node[index].contig[i].end, start + node[index].length - node[index].contig[i].start);
				fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
			}
		}

		start += node[index].length - this->minOverlap;
	}

	return numContig;
}

void PairedDBG::cutAndPrintSeq(const long minSeqLength, const unsigned long long readLength, const string &outFilename, const string &componentFilename)
{
    std::ofstream out(outFilename);
//    std::ofstream com(componentFilename);
    vector<long> maxLeftOverlap(numContig, 0);
    vector<long> maxRightOverlap(numContig, 0);

	vector<string> nodeName;
	long numOutputNode = 0;
    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node = this->node[nodeID];
		vector<char> nodeSeq;
		this->node2seq(node, nodeSeq);

		if (nodeSeq.size() < static_cast<unsigned long>(minSeqLength))
			continue;

		++numOutputNode;

		std::ostringstream oss;
		oss << ">seq" << numOutputNode << "_"
			<< "len" << nodeSeq.size() << "_"
			<< "cov" << calcNodeCoverage(node) << "_"
			<< "read" << readLength << "_"
			<< "maxK" << contigMaxK << endl;

		out << oss.str();

		unsigned long i;
		for (i = 0; i < nodeSeq.size(); ++i) {
			out << platanus::Bin2Char(nodeSeq[i]);
            if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
		}
        if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }

    out.close();
//    com.close();

    this->minOverlap = minOverlap;
}

long PairedDBG::solveUniquePathBetweenLinkedNodePair(const long numThread)
{
	vector<vector<GraphPath> > pathBuffer(numThread);
    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		for (long startNodeIndex = t; startNodeIndex < numNode; startNodeIndex += numThread) {
			GraphNode &startNode = this->node[startNodeIndex];
			for (long edgeIndex = 0; edgeIndex < startNode.numEdge; ++edgeIndex) {
				GraphEdge &edge = startNode.edge[edgeIndex];
				if (edge.state & DBG_OVERLAP)
					continue;

				searchUniqueOverlapPathGuidedByEdge(startNodeIndex + 1, edge, pathBuffer[t]);
			}
		}
	}

	vector<GraphPath> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathSelfIDLess());

	long numSolvedPath = remakeGraphAccordingToPath(mergedPathBuffer);
	return numSolvedPath;
}

void PairedDBG::searchUniqueOverlapPathGuidedByEdge(long startNodeID, const GraphEdge &guideEdge, vector<GraphPath> &pathBuffer)
{
	std::unordered_map<long, NodeInfoForPathSearch> nodeInfo;

	long endNodeID = guideEdge.end * sign(startNodeID);
	std::queue<long> nodeIDQueue;

	nodeIDQueue.push(startNodeID);
	nodeInfo[startNodeID] = NodeInfoForPathSearch(1, -(node[id2Index(startNodeID)].length), 0);
	while (nodeIDQueue.size() > 0) {
		long currentNodeID = nodeIDQueue.front();
		nodeIDQueue.pop();

		if (currentNodeID == endNodeID)
			continue;

		GraphNode &currentNode = node[id2Index(currentNodeID)];
		for (long edgeIndex = 0; edgeIndex < currentNode.numEdge; ++edgeIndex) {
			GraphEdge &edge = currentNode.edge[edgeIndex];
			if (!(edge.state & DBG_OVERLAP) || sign(currentNodeID) * edge.direction != guideEdge.direction)
				continue;

			long nextNodeID = edge.end * sign(currentNodeID);
			auto nodeInfoItr = nodeInfo.find(nextNodeID);
			long nextDistance = nodeInfo[currentNodeID].distance + currentNode.length + edge.length;
			if (nextDistance > guideEdge.length + this->tolerence)
				continue;

			if (nodeInfoItr == nodeInfo.end()) {
				nodeInfo[edge.end * sign(currentNodeID)] = NodeInfoForPathSearch(1, nextDistance, currentNodeID);
				nodeIDQueue.push(nextNodeID);
				continue;
			}

/*
			++(nodeInfoItr->second.numVisit);
			if (nextDistance < nodeInfoItr->second.distance) {
				nodeInfoItr->second.distance =  nextDistance;
				nodeInfoItr->second.preNodeID = currentNodeID;
				nodeIDQueue.push(nextNodeID);
			}
*/
		}
	}

	if (nodeInfo.size() <= 1)
		return;

	auto nodeInfoItr = nodeInfo.find(endNodeID);
	if (nodeInfoItr == nodeInfo.end())
		return;

	GraphPath path;
	while (nodeInfoItr->first != startNodeID) {
		 if (nodeInfoItr->second.numVisit != 1)
		 	return;
		path.nodeID.push_back(nodeInfoItr->first);
		nodeInfoItr = nodeInfo.find(nodeInfoItr->second.preNodeID);
	}
	path.nodeID.push_back(nodeInfoItr->first);
	if (guideEdge.direction > 0)
		std::reverse(path.nodeID.begin(), path.nodeID.end());

	pathBuffer.push_back(path);
	pathBuffer.back().sumLink = guideEdge.numLink;
	pathBuffer.back().selfID = startNodeID;
}

double PairedDBG::calcNodeCoverage(const GraphNode &node)
{
    unsigned long long sum = 0;
    unsigned long long num = 0;
    for (long i = 0; i < node.numContig; ++i) {
		unsigned long long contigIndex = id2Index(node.contig[i].id);
		if (contigPositionInScaffold[contigIndex].size() == 1) { //edited by ouchi (contigPositionInScaffold)
			num += contig[contigIndex].length;
			sum += coverage[contigIndex] * contig[contigIndex].length;
		}
    }

    //added by ouchi
    /*if (num == 0) {
        for (long i = 0; i < node.numContig; ++i) {
            unsigned long long contigIndex = id2Index(node.contig[i].id);
            if (contigPositionInScaffold[contigIndex].size() == 0) continue;
            double tmp_coverage = coverage[contigIndex];
            long num_repeat = 1;
            for (long j = 0; j < contigPositionInScaffold[contigIndex].size(); ++j) {
                if (node != this->node[id2Index(contigPositionInScaffold[contigIndex][j].id)]) {
                    const GraphNode &tmpNode = this->node[id2Index(contigPositionInScaffold[contigIndex][j].id)];
                    unsigned long long sum2 = 0;
                    unsigned long long num2 = 0;
                    for (long k = 0; k < tmpNode.numContig; ++k) {
                        unsigned long long contigIndex2 = id2Index(tmpNode.contig[k].id);
                        if (contigPositionInScaffold[contigIndex2].size() == 1) {
                            num2 += contig[contigIndex2].length;
                            sum2 += coverage[contigIndex2] * contig[contigIndex2].length;
                        }
                    }
                    if (num2 == 0)
                        num_repeat++;
                    else
                        tmp_coverage -= static_cast<double>(sum2) / num2;
                }
            }
            if (tmp_coverage > 0) {
                num += contig[contigIndex].length;
                sum += (tmp_coverage / num_repeat) * contig[contigIndex].length;
            }
        }
    }
    // */

    if (num == 0)
        return 0.0;
    else
        return static_cast<double>(sum) / num;

}

unsigned long long PairedDBG::crushHeteroBubble(const double averageCoverage)
{
    long numCrush = 0;
    const double homoCoverageThreshold = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    const double heteroCoverageThreshold = HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
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
                GraphNode *node1 = &(node[id2Index(edge1.end)]);
                if ((node1->state & SC_DEL) || edge1.length + node1->length <= edge2.length) continue;
                GraphNode *node2 = &(node[id2Index(edge2.end)]);
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

//                if (calcNodeCoverage(node[abs(layout1[leftEdge].id) - 1]) >= homoCoverageThreshold) continue;
                for (layoutID = 1; layoutID <= rightEdge; ++layoutID) {
                    if (layout1[layout1.size() - layoutID].id != layout2[layout2.size() - layoutID].id)
                        break;
                }
                if (layoutID == 1) continue;
                rightEdge = layoutID - 1;

				if (abs(layout1[leftEdge].id) == abs(layout1[layout1.size() - rightEdge].id)) continue;

                if (calcNodeCoverage(node[id2Index(layout1[layout1.size() - rightEdge].id)]) > homoCoverageThreshold) continue;
                double coverage1 = this->layoutAverageCoverage(layout1, leftEdge + 1, layout1.size() - rightEdge - leftEdge - 1);
                double coverage2 = this->layoutAverageCoverage(layout2, leftEdge + 1, layout2.size() - rightEdge - leftEdge - 1);
                const vector<GraphLayout> &layoutRef = coverage1 < coverage2 ? layout1 : layout2;
                if (rightEdge + leftEdge + 1 >= static_cast<long>(layoutRef.size()) || coverage1 > heteroCoverageThreshold || coverage2 > heteroCoverageThreshold) continue;
                long numNodeInBubble = layoutRef.size() - rightEdge - leftEdge - 1;
                if (numNodeInBubble != 1 || coverage1 > heteroCoverageThreshold || coverage2 > heteroCoverageThreshold) continue;

                layoutID = leftEdge + 1;

				node1 = &(node[id2Index(layoutRef[layoutID].id)]);
				for (long n = 0; n < node1->numEdge; ++n) {
					ids.push_back(static_cast<long>(id2Index(layoutRef[layoutID].id)) + 1);
					ids.push_back(node1->edge[n].end);
				}
				for (long n = 0; n < node1->numContig; ++n)
					contigPositionInScaffold[id2Index(node1->contig[n].id)].clear(); //edited by ouchi (contigPositionInScaffold)
				node1->state |= SC_DEL;

				if (coverage1 >= coverage2)
					std::swap(layout1, layout2);

				writeNodeSeq(node[id2Index(layout1[layoutID].id)], bubbleFP);
				writeNodeSeq(node[id2Index(layout2[layoutID].id)], bubbleOpositeFP);

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

long PairedDBG::deleteHeteroEdge(void)
{
    long numDelete = 0;
    const double homoCoverageThreshold = 2.0 * HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    const double heteroCoverageThreshold = HETERO_COVERAGE_THRESHOLD_FACTOR * heteroCoverage;
    vector<long> ids;
    if (bubbleThreshold == 0.0) return 0;

    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = node[nodeID].edge[edgeID1];
                const GraphEdge &edge2 = node[nodeID].edge[edgeID2];
                GraphNode *node1 = &(node[id2Index(edge1.end)]);
                GraphNode *node2 = &(node[id2Index(edge2.end)]);
                if (!((edge1.state & DBG_OVERLAP) && (edge2.state & DBG_OVERLAP) && edge1.direction * edge2.direction > 0) && !this->checkDeleteEdge(edge1, edge2, *node1, *node2)) continue;

                if (calcNodeCoverage(node[nodeID]) > homoCoverageThreshold) continue;

				long numExternalEdge1 = getNumEdgeOneDirection(*node1, edge1.direction * sign(edge1.end));
				long numExternalEdge2 = getNumEdgeOneDirection(*node2, edge2.direction * sign(edge2.end));
				if (numExternalEdge1 > 0 && numExternalEdge2 > 0)
					continue;

                double coverage1 = this->calcNodeCoverage(*node1);
                double coverage2 = this->calcNodeCoverage(*node2);
                long id = std::abs(edge1.end);
                if ((numExternalEdge1 > 0 && numExternalEdge2 == 0) || node1->length > node2->length) {
                    node1 = node2;
                    coverage1 = coverage2;
                    id = std::abs(edge2.end);
                }
                if (std::min(coverage1, coverage2) > heteroCoverageThreshold) continue;

                ++numDelete;
                node1->state |= SC_DEL;
                for (long edgeID = 0; edgeID < node1->numEdge; ++edgeID) {
                    ids.push_back(id);
                    ids.push_back(node1->edge[edgeID].end);
                }
                for (long contigID = 0; contigID < node1->numContig; ++contigID) {
                    contigPositionInScaffold[abs(node1->contig[contigID].id) - 1].clear(); //edited by ouchi (contigPositionInScaffold)
                }
            }
        }
    }

    this->deleteEdges(ids);
    cerr << "NUM_DELETED_HETERO_EDGES = " << numDelete << endl;

    return numDelete;
}

void PairedDBG::loadResultSeq(const long minSeqLength, const unsigned long long readLength, const string prefix)
{
	const long MIN_GAP_LENGTH = 10;
	const long MIN_OVERLAP_TO_JOIN = 32;

//	resultSeq.clear();
	resultSeq.resize(node.size());
	resultSeq.shrink_to_fit(); //added by ouchi

    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

    long maxNumContig = 0;
    long scaffoldLength = 0;
    for (long i = 0; i < numNode; ++i) {
        if (node[i].numContig > maxNumContig)
            maxNumContig = node[i].numContig;
    }
    vector<long> leftCut(maxNumContig, 0);
    vector<long> rightCut(maxNumContig, 0);
    vector<long> gap(maxNumContig, 0);

	long numOutputNode = 0;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (node[nodeIndex].state & SC_DEL)
			continue;

        long contigID = 0;
        for (; contigID < node[nodeIndex].numContig; ++contigID) {
            if (contigPositionInScaffold[id2Index(node[nodeIndex].contig[contigID].id)].size() == 1) //edited by ouchi (contigPositionInScaffold)
                break;
        }
        if (contigID == node[nodeIndex].numContig) continue;

        scaffoldLength = 0;
		leftCut[0] = 0;

        for (contigID = 0; contigID < node[nodeIndex].numContig - 1; ++contigID) {
            if (node[nodeIndex].contig[contigID].end > node[nodeIndex].contig[contigID + 1].start) {
				if (this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id) > minOverlap) {
					rightCut[contigID] = 0;
					leftCut[contigID + 1] = this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id);
					gap[contigID] = 0;
				}
				else {
					rightCut[contigID] = 0;
					leftCut[contigID + 1] = 0;
					gap[contigID] = MIN_GAP_LENGTH;
				}
            }
            else if (node[nodeIndex].contig[contigID + 1].start - node[nodeIndex].contig[contigID].end <= this->tolerence) {
				if (this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id) > minOverlap) {
					rightCut[contigID] = 0;
					leftCut[contigID + 1] = this->getOverlap(node[nodeIndex].contig[contigID].id, node[nodeIndex].contig[contigID + 1].id);
					gap[contigID] = 0;
				}
				else {
					rightCut[contigID] = leftCut[contigID + 1] = 0;
					gap[contigID] = node[nodeIndex].contig[contigID + 1].start - node[nodeIndex].contig[contigID].end;
				}
			}
			else {
                rightCut[contigID] = leftCut[contigID + 1] = 0;
                gap[contigID] = node[nodeIndex].contig[contigID + 1].start - node[nodeIndex].contig[contigID].end;
            }

            scaffoldLength += (node[nodeIndex].contig[contigID].end - node[nodeIndex].contig[contigID].start) - leftCut[contigID] - rightCut[contigID] + gap[contigID];
        }

		rightCut[contigID] = 0;
        scaffoldLength += (node[nodeIndex].contig[contigID].end - node[nodeIndex].contig[contigID].start) - leftCut[contigID] - rightCut[contigID];
        gap[contigID] = 0;

        if (scaffoldLength < minSeqLength) continue;

		resultSeq[nodeIndex].seq.clear();
		resultSeq[nodeIndex].seq.shrink_to_fit(); //added by ouchi
        for (contigID = 0; contigID < node[nodeIndex].numContig; ++contigID) {
            if (node[nodeIndex].contig[contigID].id > 0) {
                for (long seqID = leftCut[contigID]; seqID < contig[node[nodeIndex].contig[contigID].id - 1].length - rightCut[contigID]; ++seqID)
                    resultSeq[nodeIndex].seq.push_back(contig[node[nodeIndex].contig[contigID].id - 1].base[seqID]);
            }
			else {
                for (long seqID = contig[-(node[nodeIndex].contig[contigID].id) - 1].length - leftCut[contigID] - 1; seqID >= (long)rightCut[contigID]; --seqID) {
                    if (contig[-(node[nodeIndex].contig[contigID].id) - 1].base[seqID] != 4)
                        resultSeq[nodeIndex].seq.push_back(0x3 ^ contig[-(node[nodeIndex].contig[contigID].id) - 1].base[seqID]);
                    else
                        resultSeq[nodeIndex].seq.push_back(4);
                }
            }

            for (long k = 0; k < gap[contigID]; ++k)
				resultSeq[nodeIndex].seq.push_back(4);
        }

		++numOutputNode;

		if (resultSeq[nodeIndex].name.empty()) {
			std::ostringstream oss;
			oss << prefix << numOutputNode << "_"
				<< "len" << resultSeq[nodeIndex].seq.size() << "_"
				<< "cov" << calcNodeCoverage(node[nodeIndex]) << "_"
				<< "read" << readLength << "_"
				<< "maxK" << contigMaxK;

			resultSeq[nodeIndex].name = oss.str();
		}


		long currentPosition = 0;
		resultSeq[nodeIndex].component.clear();
        for (contigID = 0; contigID < node[nodeIndex].numContig; ++contigID) {

			long startPosition = currentPosition;
			long endPosition = startPosition + contig[id2Index(node[nodeIndex].contig[contigID].id)].length - leftCut[contigID] - rightCut[contigID];

			if (contigID > 0 && gap[contigID - 1] == 0)
				startPosition -= leftCut[contigID];

			if (contigID < node[nodeIndex].numContig - 1 && gap[contigID] == 0)
				endPosition += rightCut[contigID];

			std::ostringstream oss;
			oss << resultSeq[nodeIndex].name << '\t'
				<< startPosition << '\t'
				<< endPosition << '\t'
				<< contigName[id2Index(node[nodeIndex].contig[contigID].id)] << '\t'
				<< 0 << '\t';

			if (node[nodeIndex].contig[contigID].id > 0)
				oss << "+\n";
			else
				oss << "-\n";

			resultSeq[nodeIndex].component += oss.str();

			currentPosition += contig[id2Index(node[nodeIndex].contig[contigID].id)].length - leftCut[contigID] - rightCut[contigID] + gap[contigID];
		}
    }

    this->minOverlap = defaultMinOverlap;
}


void PairedDBG::loadDividedContigResultSeq(const long minSeqLength, const unsigned long long readLength)
{
	resultSeq.clear();

	long currentSeqIndex = 0;
    for (long contigIndex = 0; contigIndex < contig.size(); ++contigIndex) {
		nodeBreakpointPosition[contigIndex].push_back(0);
		nodeBreakpointPosition[contigIndex].push_back(contig[contigIndex].base.size());
		std::sort(nodeBreakpointPosition[contigIndex].begin(), nodeBreakpointPosition[contigIndex].end());

        for (auto itr = nodeBreakpointPosition[contigIndex].begin(), endItr = nodeBreakpointPosition[contigIndex].end(); itr != endItr; ++itr) {
            auto nextItr = itr;
            ++nextItr;
            if (nextItr == endItr)
				break;


            if (*nextItr - *itr < minSeqLength)
				continue;

			long start = *itr;
			long end = *nextItr;
			for (; start < end; ++start) {
				if (contig[contigIndex].base[start] != 4)
					break;
			}
			for (; end > start; --end) {
				if (contig[contigIndex].base[end - 1] != 4)
					break;
			}
            if (end - start < minSeqLength)
				continue;

			resultSeq.resize(currentSeqIndex + 1);
			resultSeq[currentSeqIndex].seq.resize(end - start);
			for (long i = start; i < end; ++i)
				resultSeq[currentSeqIndex].seq[i - start] = contig[contigIndex].base[i];

			std::ostringstream nameOss;
			nameOss << "seq" << currentSeqIndex + 1 << "_"
				<< "len" << end - start << "_"
				<< "cov" << coverage[contigIndex] << "_"
				<< "read" << readLength << "_"
				<< "maxK" << contigMaxK;
			resultSeq[currentSeqIndex].name = nameOss.str();

			std::ostringstream compOss;
			compOss << resultSeq[currentSeqIndex].name << '\t'
				<< 0 << '\t'
				<< end - start << '\t'
				<< contigName[contigIndex] << ':' << start << '-' << end << '\t'
				<< 0 << '\t'
				<< "+\n";
			resultSeq[currentSeqIndex].component += compOss.str();

			++currentSeqIndex;
		}
	}
}


void PairedDBG::outputResultSeqWithBubble(const string filePrefix, const string &primaryBubbleSuffix, const string &secondaryBubbleSuffix, const string &primaryForkSuffix, const string &secondaryForkSuffix, const string &nestedBubbleSuffix, const string &nonBubbleOtherSuffix, const string &pairSuffix, const long minSeqLength, const long readLength)
{
    cerr << "writing scaffold and bubble files..." << endl;

    string outFilename;

    outFilename = filePrefix;
    outFilename += primaryBubbleSuffix;
    std::ofstream primaryBubbleOut(outFilename);

    outFilename = filePrefix;
    outFilename += secondaryBubbleSuffix;
    std::ofstream secondaryBubbleOut(outFilename);

    outFilename = filePrefix;
    outFilename += primaryForkSuffix;
    std::ofstream primaryForkOut(outFilename);

    outFilename = filePrefix;
    outFilename += secondaryForkSuffix;
    std::ofstream secondaryForkOut(outFilename);

    outFilename = filePrefix;
    outFilename += nestedBubbleSuffix;
    std::ofstream nestedBubbleOut(outFilename);

    outFilename = filePrefix;
    outFilename += nonBubbleOtherSuffix;
    std::ofstream nonBubbleOtherOut(outFilename);

    outFilename = filePrefix;
    outFilename += pairSuffix;
    std::ofstream pairOut(outFilename);

	vector<long> scaffoldIndexOrder(0);
	vector<long> bubbleIndexOrder(0);
	vector<char> pairFlag(numNode, false);
	vector<char> unprintFlag(numNode, false);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (resultSeq[nodeIndex].seq.empty() || (node[nodeIndex].state & SC_DEL)) {
			unprintFlag[nodeIndex] = true;
			continue;
		}

		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0) {
			if (node[id2Index(oppositeNodeID)].oppositeBubbleNodeID != 0 &&
				id2Index(node[id2Index(oppositeNodeID)].oppositeBubbleNodeID) == nodeIndex &&
				checkNodeConvergenceBothSide(nodeIndex + 1, oppositeNodeID)) {

				pairFlag[nodeIndex] = true;
			}
		}
	}

	long numSeq = 0;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (!(pairFlag[nodeIndex]) || (node[nodeIndex].state & DBG_SECONDARY_BUBBLE)) {
			continue;
		}

		long altNodeIndex = id2Index(node[nodeIndex].oppositeBubbleNodeID);
		if (resultSeq[nodeIndex].redundantFlag && resultSeq[altNodeIndex].redundantFlag)
			continue;

		if (node[nodeIndex].oppositeBubbleNodeID < 0)
			reverseComplement(resultSeq[altNodeIndex].seq);

        ++numSeq;

		std::ostringstream primaryBubbleSS;
		primaryBubbleSS << "primary_bubble" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK;
		primaryBubbleOut << '>' << primaryBubbleSS.str() << '\n';
		resultSeq[nodeIndex].name = primaryBubbleSS.str();
		printBinSeqConstLineLength(resultSeq[nodeIndex].seq, primaryBubbleOut);

		std::ostringstream secondaryBubbleSS;
		secondaryBubbleSS << "secondary_bubble" << numSeq << "_len" << resultSeq[altNodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[altNodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK;
		secondaryBubbleOut << '>' << secondaryBubbleSS.str() << '\n';
		resultSeq[altNodeIndex].name = secondaryBubbleSS.str();
		printBinSeqConstLineLength(resultSeq[altNodeIndex].seq, secondaryBubbleOut);

		pairOut << primaryBubbleSS.str() << '\t' << secondaryBubbleSS.str() << '\n';
	}

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (unprintFlag[nodeIndex] || pairFlag[nodeIndex] || resultSeq[nodeIndex].redundantFlag)
			continue;

        ++numSeq;

		bool bubbleFlag = checkSingleNodeConvergenceBothSide(nodeIndex + 1);

		if (bubbleFlag && (node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))) {
			std::ostringstream oss;
			oss << "non_bubble_hetero" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK;
			nestedBubbleOut << '>' << oss.str() << '\n';
			resultSeq[nodeIndex].name = oss.str();
			printBinSeqConstLineLength(resultSeq[nodeIndex].seq, nestedBubbleOut);
		}
		else if (!bubbleFlag && (node[nodeIndex].state & DBG_PRIMARY_BUBBLE)) {
			std::ostringstream oss;
			oss << "primary_fork" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK;
			primaryForkOut << '>' << oss.str() << '\n';
			resultSeq[nodeIndex].name = oss.str();
			printBinSeqConstLineLength(resultSeq[nodeIndex].seq, primaryForkOut);
		}
		else if (!bubbleFlag && (node[nodeIndex].state & DBG_SECONDARY_BUBBLE)) {
			std::ostringstream oss;
			oss << "secondary_fork" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK;
			secondaryForkOut << '>' << oss.str() << '\n';
			resultSeq[nodeIndex].name = oss.str();
			printBinSeqConstLineLength(resultSeq[nodeIndex].seq, secondaryForkOut);
		}
		else {
			std::ostringstream oss;
			oss << "non_bubble_other" << numSeq << "_len" << resultSeq[nodeIndex].seq.size() << "_cov" << static_cast<long>(calcNodeCoverage(node[nodeIndex]) + 0.5) << "_read" << readLength << "_maxK" << contigMaxK;
			nonBubbleOtherOut << '>' << oss.str() << '\n';
			resultSeq[nodeIndex].name = oss.str();
			printBinSeqConstLineLength(resultSeq[nodeIndex].seq, nonBubbleOtherOut);
		}
	}

	primaryBubbleOut.close();
	secondaryBubbleOut.close();
	primaryForkOut.close();
	secondaryForkOut.close();
	nestedBubbleOut.close();
	nonBubbleOtherOut.close();
	pairOut.close();
}

void PairedDBG::outputResultSeqComponent(const string filePrefix, const string &fileSuffix)
{
    string outFilename;

    outFilename = filePrefix;
    outFilename += fileSuffix;
    std::ofstream out(outFilename);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (!(resultSeq[nodeIndex].seq.empty() || (node[nodeIndex].state & SC_DEL)))
			out << resultSeq[nodeIndex].component;
	}

	out.close();
}

bool PairedDBG::checkNodeConvergenceBothSide(const long nodeID1, const long nodeID2)
{
	GraphNode &node1 = this->node[id2Index(nodeID1)];
	GraphNode &node2 = this->node[id2Index(nodeID2)];

	long leftOppositeContigID1;
	if (node1.contig.front().id > 0)
		leftOppositeContigID1 = this->contigBubbleInfo[id2Index(node1.contig.front().id)].oppositeContigID[0];
	else
		leftOppositeContigID1 = -1 * this->contigBubbleInfo[id2Index(node1.contig.front().id)].oppositeContigID[1];

	long rightOppositeContigID1;
	if (node1.contig.back().id > 0)
		rightOppositeContigID1 = this->contigBubbleInfo[id2Index(node1.contig.back().id)].oppositeContigID[1];
	else
		rightOppositeContigID1 = -1 * this->contigBubbleInfo[id2Index(node1.contig.back().id)].oppositeContigID[0];

	long leftOppositeContigID2;
	if (node2.contig.front().id > 0)
		leftOppositeContigID2 = this->contigBubbleInfo[id2Index(node2.contig.front().id)].oppositeContigID[0];
	else
		leftOppositeContigID2 = -1 * this->contigBubbleInfo[id2Index(node2.contig.front().id)].oppositeContigID[1];

	long rightOppositeContigID2;
	if (node2.contig.back().id > 0)
		rightOppositeContigID2 = this->contigBubbleInfo[id2Index(node2.contig.back().id)].oppositeContigID[1];
	else
		rightOppositeContigID2 = -1 * this->contigBubbleInfo[id2Index(node2.contig.back().id)].oppositeContigID[0];


	if (leftOppositeContigID1 * rightOppositeContigID1 * leftOppositeContigID2 * rightOppositeContigID2 == 0)
		return false;


	long leftOppositeNodeID1 = 0; //edited by ouchi (contigPositionInScaffold)
	if (this->contigPositionInScaffold[id2Index(leftOppositeContigID1)].size() == 1) //edited by ouchi (contigPositionInScaffold)
		leftOppositeNodeID1 = abs(this->contigPositionInScaffold[id2Index(leftOppositeContigID1)][0].id); //edited by ouchi (contigPositionInScaffold)
	long rightOppositeNodeID1 = 0; //edited by ouchi (contigPositionInScaffold)
	if (this->contigPositionInScaffold[id2Index(rightOppositeContigID1)].size() == 1) //edited by ouchi (contigPositionInScaffold)
		rightOppositeNodeID1 = abs(this->contigPositionInScaffold[id2Index(rightOppositeContigID1)][0].id); //edited by ouchi (contigPositionInScaffold)

	long leftOppositeNodeID2 = 0; //edited by ouchi (contigPositionInScaffold)
	if (this->contigPositionInScaffold[id2Index(leftOppositeContigID2)].size() == 1) //edited by ouchi (contigPositionInScaffold)
		leftOppositeNodeID2 = abs(this->contigPositionInScaffold[id2Index(leftOppositeContigID2)][0].id); //edited by ouchi (contigPositionInScaffold)
	long rightOppositeNodeID2 = 0; //edited by ouchi (contigPositionInScaffold)
	if (this->contigPositionInScaffold[id2Index(rightOppositeContigID2)].size() == 1) //edited by ouchi (contigPositionInScaffold)
		rightOppositeNodeID2 = abs(this->contigPositionInScaffold[id2Index(rightOppositeContigID2)][0].id); //edited by ouchi (contigPositionInScaffold)

	if (leftOppositeNodeID1 == abs(nodeID2) && rightOppositeNodeID1 == abs(nodeID2) &&
		leftOppositeNodeID2 == abs(nodeID1) && rightOppositeNodeID2 == abs(nodeID1)) {
		return true;
	}
	else {
		return false;
	}
}

bool PairedDBG::checkSingleNodeConvergenceBothSide(const long nodeID)
{
	GraphNode &node1 = this->node[id2Index(nodeID)];

	long leftOppositeContigID1;
	if (node1.contig.front().id > 0)
		leftOppositeContigID1 = this->contigBubbleInfo[id2Index(node1.contig.front().id)].oppositeContigID[0];
	else
		leftOppositeContigID1 = -1 * this->contigBubbleInfo[id2Index(node1.contig.front().id)].oppositeContigID[1];

	long rightOppositeContigID1;
	if (node1.contig.back().id > 0)
		rightOppositeContigID1 = this->contigBubbleInfo[id2Index(node1.contig.back().id)].oppositeContigID[1];
	else
		rightOppositeContigID1 = -1 * this->contigBubbleInfo[id2Index(node1.contig.back().id)].oppositeContigID[0];

	if (leftOppositeContigID1 * rightOppositeContigID1 == 0)
		return false;

	long leftOppositeNodeID1 = 0; //edited by ouchi (contigPositionInScaffold)
	if (this->contigPositionInScaffold[id2Index(leftOppositeContigID1)].size() == 1) //edited by ouchi (contigPositionInScaffold)
		leftOppositeNodeID1 = abs(this->contigPositionInScaffold[id2Index(leftOppositeContigID1)][0].id); //edited by ouchi (contigPositionInScaffold)
	long rightOppositeNodeID1 = 0; //edited by ouchi (contigPositionInScaffold)
	if (this->contigPositionInScaffold[id2Index(rightOppositeContigID1)].size() == 1) //edited by ouchi (contigPositionInScaffold)
		rightOppositeNodeID1 = abs(this->contigPositionInScaffold[id2Index(rightOppositeContigID1)][0].id); //edited by ouchi (contigPositionInScaffold)

	if (leftOppositeNodeID1 * rightOppositeNodeID1 == 0 || leftOppositeNodeID1 != rightOppositeNodeID1)
		return false;

	return true;
}

void PairedDBG::solveSimpleCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
	const unsigned long currentMinLink = this->minLink;

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a de Bruijn graph.." << endl;
	do {
		setMinLink(1);
		makeGraph(numThread);

		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);

		makeHiCLink(numThread); //added by ouchi for test
		//makeMPLink(numThread, this->targetLibraryIndex, this->targetLibraryIndex + 1); //added by ouchi for test

		if (resolutionMode == SCORE)
			num = solveSimpleCrossStructure(CROSS_SCORE_RATE_THRESHOLD, resolutionMode, numThread);
		else
			num = solveSimpleCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::solveSimpleCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
	const unsigned long currentMinLink = this->minLink;

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a de Bruijn graph.." << endl;
	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);

		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);

		makeHiCLink(numThread); //added by ouchi for test
		//makeMPLink(numThread, 0, this->targetLibraryIndex + 1); //added by ouchi for test

		num = solveSimpleCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::solveSimpleGappedCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
	const unsigned long currentMinLink = this->minLink;

	long total = 0;
	long num;
	long iteration = 0;
	cerr << endl << "solving simple cross-structurs in a scaffold graph using all libraries..." << endl;
	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);
		if (!((*allLibraryMT).empty()))
			deleteEdgeFromShortNode(2 * this->tolerenceFactor * (*allLibraryMT)[targetLibraryIndex][0].getSDInsSize());
		deleteLongEdge((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize());

		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);
		if (resolutionMode == HIC) //added by ouchi
			makeHiCLink(numThread); //added by ouchi

		makeHiCLink(numThread); //added by ouchi for test
		//makeMPLink(numThread, 0, this->targetLibraryIndex + 1); //added by ouchi for test

		if (resolutionMode == SCORE)
			num = solveSimpleGappedCrossStructure(CROSS_SCORE_RATE_THRESHOLD, resolutionMode, numThread);
		else
			num = solveSimpleGappedCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::solveSimpleGappedCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread)
{
	long total = 0;
	long num;
	long iteration = 0;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "solving simple cross-structurs in a scaffold graph..." << endl;

	do {
		setMinLink(1);
		makeGraph(numThread);
		setOppositeBubbleNodeIDForEachNode(numThread);
		setMinLink(currentMinLink);
		if (resolutionMode == TAG)
			countMappedTagForEachScaffold(numThread);
		if (resolutionMode == HIC) //added by ouchi
			makeHiCLink(numThread); //added by ouchi

		makeHiCLink(numThread); //added by ouchi for test
		//makeMPLink(numThread, this->targetLibraryIndex, (*allLibraryMT).size()); //this->targetLibraryIndex + 1); //added by ouchi for test

		if (resolutionMode == SCORE)
			num = solveSimpleGappedCrossStructure(CROSS_SCORE_RATE_THRESHOLD, resolutionMode, numThread);
		else
			num = solveSimpleGappedCrossStructure(CROSS_LINK_RATE_THRESHOLD, resolutionMode, numThread);

		total += num;
		++iteration;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0 && iteration < MAX_ITERATION_OF_CROSS_SOLUTION);

	cerr << "TOTAL_NUM_NODES_IN_SOLVED_CROSSES=" << total << endl << endl;
}

void PairedDBG::joinUnambiguousNodePairIterative(const long numThread)
{
//	const double MAX_HETERO_COVERAGE_FACTOR = 1.5;

	unsigned currentMode = getMode();
	setMode(OVERLAP_MODE);

	long total = 0;
	long num;

	cerr << endl << "joining unambiguous pair of nodes in a de Bruijn graph.." << endl;
	do {
		makeGraph(numThread);
		extractDBGBubbleInformation();
//		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
//		deleteNonOverlapHomoEdge();
		num = joinUnambiguousNodePair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
	setMode(currentMode);
}

void PairedDBG::joinUnambiguousNodePairGappedIterativeAllLibraries(const long numThread)
{
	long total = 0;
	long num;
	long iteration = 0;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "joining unambiguous pair of nodes in a scaffold graph using all libraries.." << endl;

	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);
		setMinLink(currentMinLink);
//		if (!(libraryMT.empty()))
//			pairedDBG.deleteEdgeFromShortNode(2 * this->tolerenceFactor * (*allLibraryMT)[targetLibraryIndex][0].getSDInsSize());
//		pairedDBG.deleteLongEdge((*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize());

		num = joinUnambiguousNodePairGapped(numThread);
		total += num;
		++iteration;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}

void PairedDBG::solveUniquePathBetweenLinkedNodePairIterative(const long numThread)
{
	const double MAX_HETERO_COVERAGE_FACTOR = HETERO_COVERAGE_THRESHOLD_FACTOR;
	const double MAX_ITERATION = 2;

	long total = 0;
	long num;
	long numIteration = 0;
	cerr << endl << "solving unique paths guided by linked path in a de Bruijn graph.." << endl;
	do {
		makeGraph(numThread);
		extractDBGBubbleInformation();
		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
		deleteNonOverlapHomoEdge();
		num = solveUniquePathBetweenLinkedNodePair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;

		++numIteration;
	} while (num > 0 && numIteration < MAX_ITERATION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_PATHS =" << total << endl << endl;
}

void PairedDBG::solveUniquePathBetweenLinkedNodePairAllLibrariesIterative(const long numThread)
{
	const double MAX_HETERO_COVERAGE_FACTOR = HETERO_COVERAGE_THRESHOLD_FACTOR;
	const double MAX_ITERATION = 2;

	long total = 0;
	long num;
	long numIteration = 0;
	cerr << endl << "solving unique paths guided by linked path in a de Bruijn graph using all libraries.." << endl;
	do {
		makeGraphAllLibraries(numThread);
		extractDBGBubbleInformation();
		markHeteroNode(MAX_HETERO_COVERAGE_FACTOR);
		deleteNonOverlapHomoEdge();
		num = solveUniquePathBetweenLinkedNodePair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;

		++numIteration;
	} while (num > 0 && numIteration < MAX_ITERATION);
	cerr << "TOTAL_NUM_NODES_IN_SOLVED_PATHS =" << total << endl << endl;
}

void PairedDBG::setOppositeBubbleContigIDOverlapped(const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigBubbleInfo.empty())
		this->contigBubbleInfo.resize(this->numContig);

	vector<pair<long, long> > bubbleNodePairID;
	getOverlappedBubbleNodePairID(bubbleNodePairID);

	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		GraphNode &node1 = this->node[id2Index(itr->first)];
		GraphNode &node2 = this->node[id2Index(itr->second)];

		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD || calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
			continue;

		long contigID2 = sign(itr->second) * node2.contig.front().id;
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				contigID2 = sign(itr->second) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			long contigID1 = sign(itr->first) * contigItr->id;
			long contigIndex1 = id2Index(contigID1);
			if (this->contigPositionInScaffold[contigIndex1].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				this->contigBubbleInfo[contigIndex1].oppositeContigID.fill(sign(contigID1) * contigID2);
			}
		}

		long contigID1 = sign(itr->first) * node1.contig.front().id;
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				contigID1 = sign(itr->first) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			long contigID2 = sign(itr->second) * contigItr->id;
			long contigIndex2 = id2Index(contigID2);
			if (this->contigPositionInScaffold[contigIndex2].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				this->contigBubbleInfo[contigIndex2].oppositeContigID.fill(sign(contigID2) * contigID1);
			}
		}
	}
}

void PairedDBG::setOppositeBubbleContigIDGapped(const long numThread)
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigBubbleInfo.empty())
		this->contigBubbleInfo.resize(this->numContig);

	vector<pair<long, long> > bubbleNodePairID;
	getGappedBubbleNodePairID(bubbleNodePairID);

	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		GraphNode &node1 = this->node[id2Index(itr->first)];
		GraphNode &node2 = this->node[id2Index(itr->second)];

		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD || calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
			continue;

		long contigID2 = sign(itr->second) * node2.contig[0].id;
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				contigID2 = sign(itr->second) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			long contigID1 = sign(itr->first) * contigItr->id;
			if (this->contigPositionInScaffold[id2Index(contigID1)].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				this->contigBubbleInfo[id2Index(contigID1)].oppositeContigID.fill(sign(contigID1) * contigID2);
			}
		}

		long contigID1 = sign(itr->first) * node1.contig[0].id;
		for (auto contigItr = node1.contig.begin(); contigItr != node1.contig.end(); ++contigItr) {
			if (this->contigPositionInScaffold[id2Index(contigItr->id)].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				contigID1 = sign(itr->first) * contigItr->id;
				break;
			}
		}
		for (auto contigItr = node2.contig.begin(); contigItr != node2.contig.end(); ++contigItr) {
			long contigID2 = sign(itr->second) * contigItr->id;
			if (this->contigPositionInScaffold[id2Index(contigID2)].size() == 1) { //edited by ouchi (contigPositionInScaffold)
				this->contigBubbleInfo[id2Index(contigID2)].oppositeContigID.fill(sign(contigID2) * contigID1);
			}
		}
	}
}

void PairedDBG::setOppositeForkContigIDOverlapped(const long numThread)
{
//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;
	const double COVERAGE_THRESHOLD = HETERO_FORK_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigBubbleInfo.empty())
		this->contigBubbleInfo.resize(this->numContig);

	vector<pair<long, long> > bubbleNodePairID;
	vector<char> directionBuffer;
	getOverlappedForkNodePairID(bubbleNodePairID, directionBuffer);

	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		GraphNode &node1 = this->node[id2Index(itr->first)];
		GraphNode &node2 = this->node[id2Index(itr->second)];

		if (calcNodeCoverage(node1) > COVERAGE_THRESHOLD || calcNodeCoverage(node2) > COVERAGE_THRESHOLD)
			continue;

		long contigID1;
		if (directionBuffer[itr - bubbleNodePairID.begin()] * sign(itr->first) > 0)
			contigID1 = sign(itr->first) * node1.contig.front().id;
		else
			contigID1 = sign(itr->first) * node1.contig.back().id;

		long contigID2;
		if (directionBuffer[itr - bubbleNodePairID.begin()] * sign(itr->second) > 0)
			contigID2 = sign(itr->second) * node2.contig.front().id;
		else
			contigID2 = sign(itr->second) * node2.contig.back().id;


		if (directionBuffer[itr - bubbleNodePairID.begin()] * sign(itr->first) > 0) {
			if (node1.contig.front().id > 0)
				this->contigBubbleInfo[id2Index(contigID1)].oppositeContigID[0] = sign(contigID1) * contigID2;
			else
				this->contigBubbleInfo[id2Index(contigID1)].oppositeContigID[1] = sign(contigID1) * contigID2;
		}
		else {
			if (node1.contig.back().id > 0)
				this->contigBubbleInfo[id2Index(contigID1)].oppositeContigID[1] = sign(contigID1) * contigID2;
			else
				this->contigBubbleInfo[id2Index(contigID1)].oppositeContigID[0] = sign(contigID1) * contigID2;
		}

		if (directionBuffer[itr - bubbleNodePairID.begin()] * sign(itr->second) > 0) {
			if (node2.contig.front().id > 0)
				this->contigBubbleInfo[id2Index(contigID2)].oppositeContigID[0] = sign(contigID2) * contigID1;
			else
				this->contigBubbleInfo[id2Index(contigID2)].oppositeContigID[1] = sign(contigID2) * contigID1;
		}
		else {
			if (node2.contig.back().id > 0)
				this->contigBubbleInfo[id2Index(contigID2)].oppositeContigID[1] = sign(contigID2) * contigID1;
			else
				this->contigBubbleInfo[id2Index(contigID2)].oppositeContigID[0] = sign(contigID2) * contigID1;
		}
	}
}

long PairedDBG::divideNodeUsingBubbleContigPair(const long numThread)
{
//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    cerr << "dividing scaffolds to adjust bubble-breakpoints..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
	long totalNumDivision = 0;

	vector<long> numContigUsed(this->numContig, 0);

    omp_set_num_threads(numThread);
//    pragma omp parallel for schedule(dynamic) reduction(+: numNewNode, newContigPoolSize, totalNumDivision)
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		vector<std::array<long, 2> > oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, contigRef);

		pair<long, long> ends(0, oppositeBubbleNodeID.size());
		vector<pair<long, long> > endsStack(1, ends);

		while (endsStack.size() > 0) {
			ends = endsStack.back();
			endsStack.pop_back();

			pair<long, long> newEnds = fillMajorityIDRunConvergenceAware(oppositeBubbleNodeID, this->node[nodeIndex], ends, 1.0);
			if (newEnds != ends) {
				endsStack.emplace_back(ends.first, newEnds.first);
				endsStack.emplace_back(newEnds.second, ends.second);
			}
		}

		vector<char> breakpointFlag(oppositeBubbleNodeID.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		for (unsigned long i = 1; i < oppositeBubbleNodeID.size(); ++i) {
			if (oppositeBubbleNodeID[i - 1] != oppositeBubbleNodeID[i]) {
				breakpointFlag[i] = 1;
				++totalNumDivision;
			}
		}

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
	//			pragma omp critical
	//			{
					fwrite(&k, sizeof(long), 1, scaffoldFP);
					for (k = j; k < i; ++k) {
						fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
						++numContigUsed[id2Index(contigRef[k].id)];
					}
	//			}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return totalNumDivision;
}

long PairedDBG::divideNodeUsingBubbleContigPairStrandAware(const long numThread)
{
    cerr << "dividing scaffolds to adjust bubble-breakpoints distiguishing strands ..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
	long totalNumDivision = 0;

	vector<long> numContigUsed(this->numContig, 0);

    omp_set_num_threads(numThread);
//    pragma omp parallel for schedule(dynamic) reduction(+: numNewNode, newContigPoolSize, totalNumDivision)
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		vector<std::array<long, 2> > oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, contigRef);

		pair<long, long> ends(0, oppositeBubbleNodeID.size());
		vector<pair<long, long> > endsStack(1, ends);

		fillMajorityIDRunConvergenceAware(oppositeBubbleNodeID, this->node[nodeIndex], ends, 1.0);
		flipOppositeBubbleNodeID(oppositeBubbleNodeID, contigRef);

		while (endsStack.size() > 0) {
			ends = endsStack.back();
			endsStack.pop_back();

			pair<long, long> newEnds = fillMajorityIDRunConvergenceAware(oppositeBubbleNodeID, this->node[nodeIndex], ends, 1.0);
			if (newEnds != ends) {
				endsStack.emplace_back(ends.first, newEnds.first);
				endsStack.emplace_back(newEnds.second, ends.second);
			}
		}

		vector<char> breakpointFlag(oppositeBubbleNodeID.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		for (unsigned long i = 0; i < oppositeBubbleNodeID.size(); ++i) {
			if (oppositeBubbleNodeID[i] == oppositeBubbleNodeID.back()) {
				breakpointFlag[i] = 1;
				if (i != 0)
					++totalNumDivision;
				break;
			}
		}

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
	//			pragma omp critical
	//			{
					fwrite(&k, sizeof(long), 1, scaffoldFP);
					for (k = j; k < i; ++k) {
						fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
						++numContigUsed[id2Index(contigRef[k].id)];
					}
	//			}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return totalNumDivision;
}

void PairedDBG::setOppositeBubbleNodeID(vector<std::array<long, 2> > &nodeIDVector, const vector<ScaffoldPart> &partVector)
{
	nodeIDVector.resize(partVector.size());
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		for (long i = 0; i < 2; ++i) {
			long j = partItr->id > 0 ? i : (i + 1)%2;

			long tmpNodeID = 0; //edited by ouchi (contigPositionInScaffold)
			if (this->contigPositionInScaffold[id2Index(partItr->id)].size() == 1) //edited by ouchi (contigPositionInScaffold)
				tmpNodeID = this->contigPositionInScaffold[id2Index(partItr->id)][0].id; //edited by ouchi (contigPositionInScaffold)
			long tmpOppositeNodeID = 0; //edited by ouchi (contigPositionInScaffold)
			if (this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i] != 0) { //added by ouchi
				if (this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])].size() == 1) //edited by ouchi (contigPositionInScaffold)
					tmpOppositeNodeID = this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])][0].id; //edited by ouchi (contigPositionInScaffold)
			} //added by ouchi

			if (this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i] == 0 || abs(tmpNodeID) == abs(tmpOppositeNodeID)) //edited by ouchi (contigPositionInScaffold)
				nodeIDVector[partItr - partVector.begin()][j] = 0;
			else
				nodeIDVector[partItr - partVector.begin()][j] = abs(tmpOppositeNodeID); //edited by ouchi (contigPositionInScaffold)
		}
	}
}

void PairedDBG::setOppositeBubbleNodeIDStrandAware(vector<std::array<long, 2> > &nodeIDVector, const vector<ScaffoldPart> &partVector)
{
	nodeIDVector.resize(partVector.size());
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		for (unsigned i = 0; i < 2; ++i) {
			long j = partItr->id > 0 ? i : (i + 1)%2;

			long tmpNodeID = 0; //edited by ouchi (contigPositionInScaffold)
			if (this->contigPositionInScaffold[id2Index(partItr->id)].size() == 1) //edited by ouchi (contigPositionInScaffold)
				tmpNodeID = this->contigPositionInScaffold[id2Index(partItr->id)][0].id; //edited by ouchi (contigPositionInScaffold)
			long tmpOppositeNodeID = 0; //edited by ouchi (contigPositionInScaffold)
			if (this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i] != 0) { //added by ouchi
				if (this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])].size() == 1) //edited by ouchi (contigPositionInScaffold)
					tmpOppositeNodeID = this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])][0].id; //edited by ouchi (contigPositionInScaffold)
			} //added by ouchi

			if (this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i] == 0 || abs(tmpNodeID) == abs(tmpOppositeNodeID)) //edited by ouchi (contigPositionInScaffold)
				nodeIDVector[partItr - partVector.begin()][j] = 0;
			else
				nodeIDVector[partItr - partVector.begin()][j] = tmpOppositeNodeID; //edited by ouchi (contigPositionInScaffold)
		}
	}
}

void PairedDBG::flipOppositeBubbleNodeID(vector<std::array<long, 2> > &nodeIDVector, const vector<ScaffoldPart> &partVector)
{
	for (auto partItr = partVector.begin(); partItr != partVector.end(); ++partItr) {
		for (unsigned i = 0; i < 2; ++i) {
			long oppositeNodeID;
			long j = partItr->id > 0 ? i : (i + 1)%2;

			if (this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i] != 0)
				if (this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])].size() == 1) //edited by ouchi (contigPositionInScaffold)
					oppositeNodeID = this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])][0].id; //edited by ouchi (contigPositionInScaffold)
				else //edited by ouchi (contigPositionInScaffold)
					oppositeNodeID = 0; //edited by ouchi (contigPositionInScaffold)
			else
				oppositeNodeID = 0;

			if (oppositeNodeID == 0 || abs(oppositeNodeID) == abs(nodeIDVector[partItr - partVector.begin()][j]))
				nodeIDVector[partItr - partVector.begin()][j] = sign(partItr->id) * oppositeNodeID;
		}
	}
}

double PairedDBG::calcNodeCoveragePartial(const GraphNode &node, const long start, const long end)
{
    unsigned long long sum = 0;
    unsigned long long num = 0;
    for (long i = start; i < end; ++i) {
		unsigned long long contigIndex = id2Index(node.contig[i].id);
		if (contigPositionInScaffold[contigIndex].size() == 1) { //edited by ouchi (contigPositionInScaffold)
			num += contig[contigIndex].length;
			sum += coverage[contigIndex] * contig[contigIndex].length;
		}
    }

    if (num == 0)
        return 0.0;
    else
        return static_cast<double>(sum) / num;
}


long PairedDBG::maxLengthContigID(vector<std::array<long, 2 > > &IDs, const long start, const long end)
{
	std::unordered_map<long, long> countMap;

	for (long i = start; i < end; ++i) {
		if (IDs[i][0] * IDs[i][1] == 0 || IDs[i][0] != IDs[i][1])
			continue;

		auto countItr = countMap.find(IDs[i][0]);
		if (countItr != countMap.end())
			countItr->second += 2 * this->node[id2Index(IDs[i][0])].length;
		else
			countMap[IDs[i][0]] = 2 * this->node[id2Index(IDs[i][0])].length;
	}

	if (countMap.empty()) {
		for (long i = start; i < end; ++i) {
			for (long j = 0; j < 2; ++j) {
				if (IDs[i][j] == 0)
					continue;

				auto countItr = countMap.find(IDs[i][j]);
				if (countItr != countMap.end())
					countItr->second += this->node[id2Index(IDs[i][j])].length;
				else
					countMap[IDs[i][j]] = this->node[id2Index(IDs[i][j])].length;
			}
		}
	}

	long maxCount = 0;
	long maxID = 0;
	for (auto countItr = countMap.begin(); countItr != countMap.end(); ++countItr) {
		if (maxCount < countItr->second) {
			maxID = countItr->first;
			maxCount = countItr->second;
		}
	}

	return maxID;
}

pair<long, long> PairedDBG::fillMajorityIDRunConvergenceAware(vector<std::array<long, 2> > &IDs, const GraphNode &targetNode , const pair<long, long> &ends, const double scoreFactor)
{
	long maxID = maxLengthContigID(IDs, ends.first, ends.second);

	pair<long, long> newEnds(ends);
	if (maxID == 0)
		return newEnds;

	for (long i = newEnds.first; i < newEnds.second; ++i) {
		if (IDs[i][0] == maxID) {
			newEnds.first = i;
			break;
		}
	}
	for (long i = newEnds.second - 1; i >= newEnds.first; --i) {
		if (IDs[i][1] == maxID) {
			newEnds.second = i + 1;
			break;
		}
	}
	if (newEnds.first >= newEnds.second) {
		newEnds.first = ends.first;
		newEnds.second = ends.first + 1;
		return newEnds;
	}


	long score = 0;
	for (long i = newEnds.first; i < newEnds.second; ++i) {
		if (score > 0) {
			score = 0;
			newEnds.first = i;
		}

		if (IDs[i][0] * IDs[i][1] != 0 && IDs[i][0] == IDs[i][1]) {
			if (IDs[i][0] == maxID)
				score -= 2 * this->node[id2Index(IDs[i][0])].length;
			else
				score += 2 * this->node[id2Index(IDs[i][0])].length * scoreFactor;
		}
		else {
			for (long j = 0; j < 2; ++j) {
				if (IDs[i][j] == maxID)
					score -= this->node[id2Index(IDs[i][j])].length;
			}
		}
	}

	score = 0;
	for (long i = newEnds.second - 1; i >= newEnds.first; --i) {
		if (score > 0) {
			score = 0;
			newEnds.second = i + 1;
		}

		if (IDs[i][0] * IDs[i][1] != 0 && IDs[i][0] == IDs[i][1]) {
			if (IDs[i][0] == maxID)
				score -= 2 * this->node[id2Index(IDs[i][0])].length;
			else
				score += 2 * this->node[id2Index(IDs[i][0])].length * scoreFactor;
		}
		else {
			for (long j = 0; j < 2; ++j) {
				if (IDs[i][j] == maxID)
					score -= this->node[id2Index(IDs[i][j])].length;
			}
		}
	}


	for (long i = newEnds.first; i < newEnds.second; ++i) {
		if (IDs[i][0] == maxID) {
			newEnds.first = i;
			break;
		}
	}
	for (long i = newEnds.second - 1; i >= newEnds.first; --i) {
		if (IDs[i][1] == maxID) {
			newEnds.second = i + 1;
			break;
		}
	}
	if (newEnds.first >= newEnds.second) {
		newEnds.first = ends.first;
		newEnds.second = ends.first + 1;
		return newEnds;
	}

	for (auto itr = IDs.begin()+ newEnds.first; itr !=  IDs.begin() + newEnds.second - 1; ++itr)
		itr->fill(maxID);

	return newEnds;
}

void PairedDBG::setOppositeBubbleNodeIDForEachNode(const long numThread)
{
	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		this->node[nodeIndex].oppositeBubbleNodeID = 0;
		this->node[nodeIndex].state &= ~(DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE);
	}

	if (this->contigBubbleInfo.empty())
		return;

//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    cerr << "pairing bubble nodes..." << endl;

    omp_set_num_threads(numThread);
	# pragma omp parallel for schedule(static, 1)
	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		vector<std::array<long, 2> > oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, this->node[nodeIndex].contig);

		long oppositeNodeID = maxLengthContigID(oppositeBubbleNodeID, 0, oppositeBubbleNodeID.size());
//		if (oppositeNodeID == 0 || calcNodeCoverage(this->node[nodeIndex]) > COVERAGE_THRESHOLD || calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > COVERAGE_THRESHOLD)
		if (oppositeNodeID == 0 ||
			calcNodeCoverage(this->node[nodeIndex]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[nodeIndex].length)) ||
			calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[id2Index(oppositeNodeID)].length)))
			continue;

		if (id2Index(oppositeNodeID) == nodeIndex)
			continue;

		this->node[nodeIndex].oppositeBubbleNodeID = oppositeNodeID;
	}
}

void PairedDBG::setOppositeBubbleNodeIDAndStateForEachNode()
{
	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		this->node[nodeIndex].oppositeBubbleNodeID = 0;
		this->node[nodeIndex].state &= ~(DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE);
	}

	for (unsigned long i = 0; i < this->contigState.size(); ++i)
		contigState[i] &= ~(DBG_CONTIG_PRIMARY_BUBBLE | DBG_CONTIG_SECONDARY_BUBBLE);

	if (this->contigBubbleInfo.empty())
		return;

//	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    cerr << "pairing bubble nodes..." << endl;

	for (long nodeIndex = 0; nodeIndex < static_cast<long>(this->node.size()); ++nodeIndex) {
		vector<std::array<long, 2> > oppositeBubbleNodeID;
		setOppositeBubbleNodeID(oppositeBubbleNodeID, this->node[nodeIndex].contig);

		long oppositeNodeID = maxLengthContigID(oppositeBubbleNodeID, 0, oppositeBubbleNodeID.size());
//		if (oppositeNodeID == 0 || calcNodeCoverage(this->node[nodeIndex]) > COVERAGE_THRESHOLD || calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > COVERAGE_THRESHOLD)
		if (oppositeNodeID == 0 ||
			calcNodeCoverage(this->node[nodeIndex]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[nodeIndex].length)) ||
			calcNodeCoverage(this->node[id2Index(oppositeNodeID)]) > this->heteroCoverage * std::max(1.25, HETERO_COVERAGE_THRESHOLD_FACTOR - (0.25 * 0.00001 * this->node[id2Index(oppositeNodeID)].length)))
			continue;

		long oppositeNodeIndex = id2Index(oppositeNodeID);
		if (oppositeNodeIndex == nodeIndex)
			continue;

		this->node[nodeIndex].oppositeBubbleNodeID = oppositeNodeID;

		long numEdgeDirection = getNumEdgeDirectionOfNode(this->node[nodeIndex]);
		long oppositeNumEdgeDirection = getNumEdgeDirectionOfNode(this->node[oppositeNodeIndex]);
		if (numEdgeDirection > oppositeNumEdgeDirection) {
			this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
		}
		else if (numEdgeDirection < oppositeNumEdgeDirection) {
			this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
		}
		else {
			long nonGapLength = getNonGapContigLengthOfNode(this->node[nodeIndex]);
			long oppositeNonGapLength = getNonGapContigLengthOfNode(this->node[oppositeNodeIndex]);
			if (nonGapLength > oppositeNonGapLength) {
				this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
			}
			else if (nonGapLength < oppositeNonGapLength) {
				this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
			}
			else {
				double nodeCoverage = calcNodeCoverage(this->node[nodeIndex]);
				double oppositeNodeCoverage = calcNodeCoverage(this->node[oppositeNodeIndex]);
				if (nodeCoverage > oppositeNodeCoverage) {
					this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
				}
				else if (nodeCoverage < oppositeNodeCoverage) {
					this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
				}
				else {
					if (nodeIndex < oppositeNodeIndex) {
						this->node[oppositeNodeIndex].state |= DBG_SECONDARY_BUBBLE;
					}
					else {
						this->node[nodeIndex].state |= DBG_SECONDARY_BUBBLE;
					}
				}
			}
		}

		if (this->node[nodeIndex].state & DBG_SECONDARY_BUBBLE)
			this->node[oppositeNodeIndex].state |= DBG_PRIMARY_BUBBLE;
		else
			this->node[nodeIndex].state |= DBG_PRIMARY_BUBBLE;
	}


	for (unsigned long i = 0; i < this->contigState.size(); ++i) {
		if (contigPositionInScaffold[i].size() != 1) //edited by ouchi (contigPositionInScaffold)
			continue;
		else if (node[id2Index(contigPositionInScaffold[i][0].id)].state & DBG_PRIMARY_BUBBLE) //edited by ouchi (contigPositionInScaffold)
			contigState[i] |= DBG_CONTIG_PRIMARY_BUBBLE;
		else if (node[id2Index(contigPositionInScaffold[i][0].id)].state & DBG_SECONDARY_BUBBLE) //edited by ouchi (contigPositionInScaffold)
			contigState[i] |= DBG_CONTIG_SECONDARY_BUBBLE;
	}
}

long PairedDBG::getNonGapContigLengthOfNode(const GraphNode &targetNode)
{
	if (targetNode.contig.empty())
		return 0;

	long gapLength = 0;
	for (unsigned i = 0; i < targetNode.contig.size() - 1; ++i)
		gapLength += targetNode.contig[i + 1].start - targetNode.contig[i].end;

	return targetNode.contig.back().end - gapLength;
}

long PairedDBG::getNumEdgeDirectionOfNode(const GraphNode &targetNode)
{
	long numLeft = 0;
	long numRight = 0;

	for (auto itr = targetNode.edge.begin(); itr != targetNode.edge.end(); ++itr) {
		if (itr->direction > 0)
			numLeft = 1;
		else
			numRight = 1;
	}

	return (numLeft + numRight);
}

long PairedDBG::deleteDifferentBubbleEdge(const long numThread)
{
    vector<long> ids;
    long numDelete = 0;

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		const GraphNode &sourceNode = node[nodeIndex];
        for (long edgeID1 = 0; edgeID1 < sourceNode.numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeIndex].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (sourceNode.edge[edgeID1]);
                const GraphEdge &edge2 = (sourceNode.edge[edgeID2]);
                const GraphNode &node1 = (node[id2Index(edge1.end)]);
                const GraphNode &node2 = (node[id2Index(edge2.end)]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2) || sourceNode.oppositeBubbleNodeID == 0) continue;

                if (sourceNode.oppositeBubbleNodeID != sign(edge1.end)*node1.oppositeBubbleNodeID && sourceNode.oppositeBubbleNodeID == sign(edge2.end)*node2.oppositeBubbleNodeID) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeIndex + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (sourceNode.oppositeBubbleNodeID != sign(edge2.end)*node2.oppositeBubbleNodeID && sourceNode.oppositeBubbleNodeID == sign(edge2.end)*node1.oppositeBubbleNodeID) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeIndex + 1);
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

void PairedDBG::deleteDifferentBubbleEdgeIterative(const long numThread)
{
    cerr << "removing edges between nodes with different contig-bubble assignments..." << endl << endl;;

	this->setOppositeBubbleNodeIDForEachNode(numThread);

    long totalDelete = 0;
    long numDelete;
    do {
        numDelete = this->deleteDifferentBubbleEdge(numThread);
        totalDelete += numDelete;
        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
    } while (numDelete > 0);

    cerr << "TOTAL_REMOVED_EDGES =" << totalDelete << endl << endl;;
}

void PairedDBG::deleteThinEdgeCostantKeepingOverlap(const long linkThreshold)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing thin edges (NUM_LINK < " <<  linkThreshold << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (!(node1.edge[edgeID].state & DBG_OVERLAP) && node1.edge[edgeID].numLink < linkThreshold) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::deleteConflictingBubbleEdge(const long numThread)
{
	const double CROSS_LINK_RATE_THRESHOLD = 0.25;

    omp_set_num_threads(numThread);
	this->setOppositeBubbleNodeIDForEachNode(numThread);

    vector<long> ids;
    long numDelete = 0;

    cerr << "removing conflicting edges between bubbles..." << endl;

    #pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
		if (node[nodeID].oppositeBubbleNodeID == 0)
			continue;

        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {

                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);
                const GraphNode &node1 = (node[abs(edge1.end) - 1]);
                const GraphNode &node2 = (node[abs(edge2.end) - 1]);
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2))
					continue;

				if (node1.oppositeBubbleNodeID != 0 && node1.oppositeBubbleNodeID != sign(edge1.end) * edge2.end)
					continue;

                if (edge1.numLink < CROSS_LINK_RATE_THRESHOLD * edge2.numLink) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (edge2.numLink < CROSS_LINK_RATE_THRESHOLD * edge1.numLink) {
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
    cerr << "TOTAL_NUM_DELETE=" << numDelete << endl << endl;
    this->deleteEdges(ids);
    return numDelete;
}

long PairedDBG::deleteSecondaryBubbleNodeAndEdge(const long numThread)
{
    omp_set_num_threads(numThread);
	this->setOppositeBubbleNodeIDAndStateForEachNode();

    cerr << "removing secondary bubbles from scaffold graph..." << endl;

    vector<long> ids;
    unsigned long long numDelete = 0;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
		if (!(this->node[nodeIndex].state & DBG_SECONDARY_BUBBLE))
			continue;

        GraphNode &bubbleNode = this->node[nodeIndex];

		++numDelete;
		bubbleNode.state |= SC_DEL;
        for (long edgeIndex = 0; edgeIndex < bubbleNode.numEdge; ++edgeIndex) {
			ids.push_back(nodeIndex + 1);
			ids.push_back(bubbleNode.edge[edgeIndex].end);
        }
    }
    std::cerr << "TOTAL_NUM_DELETED_NODES =" << numDelete << std::endl;
    this->deleteEdges(ids);

	return numDelete;
}

long PairedDBG::deleteShortAndLowCoverageBranch(const long lengthThreshold, const double coverageThreshold, const long numThread)
{
    omp_set_num_threads(numThread);

    vector<long> ids;
    unsigned long long numDelete = 0;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &targetNode = this->node[nodeIndex];
		if (targetNode.length >= lengthThreshold || calcNodeCoverage(targetNode) >= coverageThreshold || getNumEdgeDirectionOfNode(this->node[nodeIndex]) >= 2)
			continue;

		++numDelete;
		targetNode.state |= SC_DEL;
        for (long edgeIndex = 0; edgeIndex < targetNode.numEdge; ++edgeIndex) {
			ids.push_back(nodeIndex + 1);
			ids.push_back(targetNode.edge[edgeIndex].end);
        }
    }
    this->deleteEdges(ids);

	return numDelete;
}

void PairedDBG::deleteShortAndLowCoverageBranchIterative(const long lengthThreshold, const double coverageThreshold, const long numThread)
{
	long total = 0;
	long num;
    cerr << "removing short and low-coverage branches..." << endl;
	do {
		joinUnambiguousNodePairIterative(numThread);

		makeGraph(numThread);
		num = deleteShortAndLowCoverageBranch(lengthThreshold, coverageThreshold, numThread);

		total += num;
		cerr << "NUM_SOLVED_CROSSES = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_REMOVED_NODES =" << total << endl << endl;
}

void PairedDBG::setBubbleJunctionContigIDOverlapped()
{
	if (this->contigBubbleInfo.empty()) {
		this->contigBubbleInfo.resize(this->numContig);
		for (auto itr = this->contigBubbleInfo.begin(); itr != this->contigBubbleInfo.end(); ++itr)
			itr->joinedContigID[0] = itr->joinedContigID[1] = 0;
	}

	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigState.empty())
		this->contigState.resize(this->numContig, 0);

	vector<pair<long, long> > bubbleNodePairID;
	getOverlappedBubbleNodePairID(bubbleNodePairID);

	vector<char> bubbleFlag(this->node.size(), false);
	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		if (calcNodeCoverage(this->node[id2Index(itr->first)]) <= COVERAGE_THRESHOLD && calcNodeCoverage(this->node[id2Index(itr->second)]) <= COVERAGE_THRESHOLD) {
			bubbleFlag[id2Index(itr->first)] = true;
			bubbleFlag[id2Index(itr->second)] = true;
		}
	}

	vector<long> nodeIDBuffer;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		for (char direction = -1; direction <= 1; direction += 2) {
			getOverlappedNode(nodeIndex, direction, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (!bubbleFlag[id2Index(nodeIDBuffer[i])])
					break;
			}
			if (i != 2)
				continue;

			if (direction > 0) {
				this->contigState[id2Index(this->node[nodeIndex].contig.back().id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
				if (this->node[nodeIndex].contig.back().id > 0) {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[1] = this->node[id2Index(nodeIDBuffer[0])].contig.front().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[1] = -(this->node[id2Index(nodeIDBuffer[0])].contig.back().id);
				}
				else {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[0] = this->node[id2Index(nodeIDBuffer[0])].contig.front().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[0] = -(this->node[id2Index(nodeIDBuffer[0])].contig.back().id);
				}
			}
			else {
				this->contigState[id2Index(this->node[nodeIndex].contig.front().id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
				if (this->node[nodeIndex].contig.front().id > 0) {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[0] = this->node[id2Index(nodeIDBuffer[0])].contig.back().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[0] = -(this->node[id2Index(nodeIDBuffer[0])].contig.front().id);
				}
				else {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[1] = this->node[id2Index(nodeIDBuffer[0])].contig.back().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[1] = -(this->node[id2Index(nodeIDBuffer[0])].contig.front().id);
				}
			}
		}
	}
}

void PairedDBG::setForkJunctionContigIDOverlapped()
{
	if (this->contigBubbleInfo.empty()) {
		this->contigBubbleInfo.resize(this->numContig);
		for (auto itr = this->contigBubbleInfo.begin(); itr != this->contigBubbleInfo.end(); ++itr)
			itr->joinedContigID[0] = itr->joinedContigID[1] = 0;
	}

	const double COVERAGE_THRESHOLD = HETERO_FORK_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigState.empty())
		this->contigState.resize(this->numContig, 0);

	vector<pair<long, long> > bubbleNodePairID;
	vector<char> directionBuffer;
	getOverlappedForkNodePairID(bubbleNodePairID, directionBuffer);

	vector<char> bubbleFlag(this->node.size(), false);
	for (auto itr = bubbleNodePairID.begin(); itr != bubbleNodePairID.end(); ++itr) {
		if (calcNodeCoverage(this->node[id2Index(itr->first)]) <= COVERAGE_THRESHOLD && calcNodeCoverage(this->node[id2Index(itr->second)]) <= COVERAGE_THRESHOLD) {
			bubbleFlag[id2Index(itr->first)] = true;
			bubbleFlag[id2Index(itr->second)] = true;
		}
	}

	vector<long> nodeIDBuffer;

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		for (char direction = -1; direction <= 1; direction += 2) {
			getOverlappedNode(nodeIndex, direction, nodeIDBuffer);
			if (nodeIDBuffer.size() != 2)
				continue;

			long i = 0;
			for (; i < 2; ++i) {
				if (!bubbleFlag[id2Index(nodeIDBuffer[i])])
					break;
			}
			if (i != 2)
				continue;

			if (direction > 0) {
				this->contigState[id2Index(this->node[nodeIndex].contig.back().id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
				if (this->node[nodeIndex].contig.back().id > 0) {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[1] = this->node[id2Index(nodeIDBuffer[0])].contig.front().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[1] = -(this->node[id2Index(nodeIDBuffer[0])].contig.back().id);
				}
				else {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[0] = this->node[id2Index(nodeIDBuffer[0])].contig.front().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.back().id)].joinedContigID[0] = -(this->node[id2Index(nodeIDBuffer[0])].contig.back().id);
				}
			}
			else {
				this->contigState[id2Index(this->node[nodeIndex].contig.front().id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
				if (this->node[nodeIndex].contig.front().id > 0) {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[0] = this->node[id2Index(nodeIDBuffer[0])].contig.back().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[0] = -(this->node[id2Index(nodeIDBuffer[0])].contig.front().id);
				}
				else {
					if (nodeIDBuffer[0] > 0)
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[1] = this->node[id2Index(nodeIDBuffer[0])].contig.back().id;
					else
						this->contigBubbleInfo[id2Index(this->node[nodeIndex].contig.front().id)].joinedContigID[1] = -(this->node[id2Index(nodeIDBuffer[0])].contig.front().id);
				}
			}
		}
	}
}

void PairedDBG::markJunctionContigJoinedToBubble()
{
	vector<char> bubbleEdgeFlag(this->numContig, 0);

	for (auto itr = this->contigState.begin(); itr != this->contigState.end(); ++itr)
		*itr &= ~DBG_CONTIG_BUBBLE_JUNCTION;

	for (auto nodeIt = node.begin(); nodeIt != node.end(); ++nodeIt) {
		if (nodeIt->state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)) {
			bubbleEdgeFlag[id2Index(nodeIt->contig.front().id)] = 1;
			bubbleEdgeFlag[id2Index(nodeIt->contig.back().id)] = 1;
		}
	}

	for (auto nodeIt = node.begin(); nodeIt != node.end(); ++nodeIt) {
		for (auto contigIt = nodeIt->contig.begin(); contigIt != nodeIt->contig.end(); ++contigIt) {
			for (long i = 0; i < 2; ++i) {
				if (contigBubbleInfo[id2Index(contigIt->id)].joinedContigID[i] != 0 && bubbleEdgeFlag[id2Index(contigBubbleInfo[id2Index(contigIt->id)].joinedContigID[i])])
					this->contigState[id2Index(contigIt->id)] |= DBG_CONTIG_BUBBLE_JUNCTION;
			}
		}
	}
}


long PairedDBG::divideBubbleJunctionNode(const bool gapDivideFlag)
{
    cerr << "dividing scaffolds at bubble-junctions..." << endl;

	this->setOppositeBubbleNodeIDAndStateForEachNode();

	this->markJunctionContigJoinedToBubble();

	const long MIN_OVERLAP_TO_JOIN = 20;
    long defaultMinOverlap = this->minOverlap;
	if (gapDivideFlag)
		this->minOverlap = MIN_OVERLAP_TO_JOIN;
	else
		this->minOverlap = this->contigMaxK - 1;

	long numDivision = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		if (gapDivideFlag) {
			if (this->node[nodeIndex].oppositeBubbleNodeID == 0) {
				for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
					if ((((this->contigState[id2Index(contigRef[i - 1].id)] | this->contigState[id2Index(contigRef[i].id)]) & DBG_CONTIG_BUBBLE_JUNCTION) && this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < this->contigMaxK) ||
						(contigRef[i].start > contigRef[i - 1].end || this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < minOverlap)) {

						breakpointFlag[i] = 1;
						++numDivision;
					}
				}
			}
		}
		else {
			for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
				if (((this->contigState[id2Index(contigRef[i - 1].id)] | this->contigState[id2Index(contigRef[i].id)]) & DBG_CONTIG_BUBBLE_JUNCTION)
					&& this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < this->contigMaxK - 1) {

					breakpointFlag[i] = 1;
					++numDivision;
				}
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->minOverlap = defaultMinOverlap;
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return numDivision;
}

long PairedDBG::divideBubbleContigInNonHeteroNode()
{
    cerr << "dividing scaffolds at bubble-contig in non-hetero scaffolds ..." << endl;

	this->setOppositeBubbleNodeIDAndStateForEachNode();

	long numDivision = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		if (this->node[nodeIndex].oppositeBubbleNodeID == 0 && contigRef.size() > 1) {
			for (i = 0; i < static_cast<long>(contigRef.size()); ++i) {
				if (this->contigBubbleInfo[id2Index(contigRef[i].id)].oppositeContigID[0] != 0 || this->contigBubbleInfo[id2Index(contigRef[i].id)].oppositeContigID[1] != 0) {
					breakpointFlag[i] = 1;
					breakpointFlag[i + 1] = 1;
					++numDivision;
				}
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return numDivision;
}

void PairedDBG::divideGappedNode(const long minGapSize)
{
    cerr << "dividing scaffolds at gaps..." << endl;

	const long MIN_OVERLAP_TO_JOIN = 20;
    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;

		for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
//			if (contigRef[i].start > contigRef[i - 1].end || this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < minOverlap)
			if (contigRef[i].start - contigRef[i - 1].end > minGapSize)
				breakpointFlag[i] = 1;
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->minOverlap = defaultMinOverlap;
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::trimSparseEnd()
{
    cerr << "trimming sparse ends of scaffolds..." << endl;

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (this->node[nodeIndex].state & SC_DEL)
			continue;

		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;


		if (contigRef.size() > 1) {
			if (contigRef[1].start - contigRef[0].end > contigRef[0].end - contigRef[0].start)
				breakpointFlag[1] = 1;

			if (contigRef[contigRef.size() - 1].start - contigRef[contigRef.size() - 2].end > contigRef[contigRef.size() - 1].end - contigRef[contigRef.size() - 1].start)
				breakpointFlag[contigRef.size() - 2] = 1;
		}


		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::copyAllNodes(PairedDBG &targetGraph)
{
    targetGraph.numNode = this->numNode;
    targetGraph.numContig = this->numContig;
    targetGraph.node= this->node;
    targetGraph.contigPositionInScaffold = this->contigPositionInScaffold;
	targetGraph.contigBubbleInfo = this->contigBubbleInfo;
	targetGraph.contigState = this->contigState;
}

void PairedDBG::smoothNodeIDVector(vector<std::array<long, 2> > &nodeIDVector, const GraphNode &targetNode, const double scoreFactor)
{
	pair<long, long> ends(0, nodeIDVector.size());
	vector<pair<long, long> > endsStack(1, ends);

	while (endsStack.size() > 0) {
		ends = endsStack.back();
		endsStack.pop_back();

		pair<long, long> newEnds = fillMajorityIDRunConvergenceAware(nodeIDVector, targetNode, ends, scoreFactor);
		if (newEnds != ends) {
			endsStack.emplace_back(ends.first, newEnds.first);
			endsStack.emplace_back(newEnds.second, ends.second);
		}
	}
}

long PairedDBG::getQuartileLengthOfBubble(const unsigned long quartileNumber)
{
	vector<long> lengthBuffer;

    for (long nodeIndex = 0; nodeIndex < this->numNode; ++nodeIndex) {
		if (this->node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))
			lengthBuffer.push_back(this->node[nodeIndex].length);
	}

    sort(lengthBuffer.begin(), lengthBuffer.end());

	for (unsigned i = 0; i < lengthBuffer.size(); ++i) {
		if (i/lengthBuffer.size() >= quartileNumber/4)
			return lengthBuffer[i];
	}

	return 1;
}

void PairedDBG::detectRepeat(const double averageCoverage)
{
    const double coverageThreshold = averageCoverage * 1.75;

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
                if ((abs(edge1.length + node1.length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
                    || abs(edge2.length + node2.length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1))
					&& abs(edge1.end) != abs(node2.oppositeBubbleNodeID)) continue;

                node[nodeID].state |= SC_REP;
                edgeID1 = node[nodeID].numEdge;
                break;
            }
        }
    }
}

void PairedDBG::deleteRepeatEdge()
{
//    const double coverageThreshold = averageCoverage * 1.75;

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

				bool collisionFlag = true;
                if (!this->checkDeleteEdge(edge1, edge2, node1, node2) && abs(edge1.end) != abs(node2.oppositeBubbleNodeID))
					collisionFlag = false;

                for (long m = 0; m < node[nodeID].numContig; ++m) {
//                    if (coverage[id2Index(node[nodeID].contig[m].id)] <= coverageThreshold && !collisionFlag)
                    if (!collisionFlag)
						continue;

					long numLinkToSubtract = std::min(edge1.breakdown[m], edge2.breakdown[m]);

                    for (long n = 0; n < node[nodeID].numEdge; ++n) {
                        node[nodeID].edge[n].numLink -= numLinkToSubtract;
                        node[nodeID].edge[n].breakdown[m] -= numLinkToSubtract;
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

void PairedDBG::deleteEdgeFromDifferentPreviousParent(const long numThread)
{
	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		return;

    vector<long> ids;

    cerr << "deleting edges from previously divided points ..." << endl;

    omp_set_num_threads(numThread);
    #pragma omp parallel for schedule(dynamic)
    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		const GraphNode &node1 = node[nodeIndex];
        for (long edgeIndex = 0; edgeIndex < node[nodeIndex].numEdge - 1; ++edgeIndex) {
			const GraphNode &node2 = node[id2Index(node[nodeIndex].edge[edgeIndex].end)];
			for (long contigIndex1 = 0; contigIndex1 < node1.numContig; ++contigIndex1) {
				if (contigPreviousParentNodeID[id2Index(node1.contig[contigIndex1].id)] == 0)
					break;

				for (long contigIndex2 = 0; contigIndex2 < node2.numContig; ++contigIndex2) {
					if (contigPreviousParentNodeID[id2Index(node2.contig[contigIndex2].id)] == 0)
						break;

					if (contigPreviousParentNodeID[id2Index(node1.contig[contigIndex1].id)] == contigPreviousParentNodeID[id2Index(node2.contig[contigIndex2].id)]) {
                        node[nodeIndex].edge[edgeIndex].numLink -= node1.edge[edgeIndex].breakdown[contigIndex1];
						node[nodeIndex].edge[edgeIndex].breakdown[contigIndex1] = 0;
					}
				}
			}
		}
	}

    long numDelete = 0;
    for (long i = 0; i < numNode; ++i) {
        for (long j = 0; j < node[i].numEdge; ++j) {
            if (node[i].edge[j].numLink < minLink) {
                ids.push_back(i + 1);
                ids.push_back(node[i].edge[j].end);
				++numDelete;
            }
        }
    }
    this->deleteEdges(ids);

	cerr << "NUM_REMOVED_EDGES =" << numDelete / 2 << endl;
}

void PairedDBG::getUniqueConflictingNode(const long sourceNodeIndex, const char targetDirection, vector<NodeIDWithGap> &nodeIDBuffer)
{
	nodeIDBuffer.clear();
	GraphNode &sourceNode = this->node[sourceNodeIndex];

	for (long edgeIndex1 = 0; edgeIndex1 < sourceNode.numEdge - 1; ++edgeIndex1) {
		GraphEdge &edge1 = sourceNode.edge[edgeIndex1];
		if (edge1.direction != targetDirection || edge1.numLink < minLink)
			continue;

		for (long edgeIndex2 = edgeIndex1 + 1; edgeIndex2 < sourceNode.numEdge; ++edgeIndex2) {
			GraphEdge &edge2 = sourceNode.edge[edgeIndex2];
			if (edge2.direction != targetDirection || edge2.numLink < minLink)
				continue;

			if (this->checkDeleteEdge(edge1, edge2, this->node[id2Index(edge1.end)], this->node[id2Index(edge2.end)]) || abs(edge1.end) == abs(this->node[id2Index(edge2.end)].oppositeBubbleNodeID)) {
				if (!nodeIDBuffer.empty()) {
					nodeIDBuffer.clear();
					return;
				}
				nodeIDBuffer.emplace_back(edge1.end, edge1.length);
				nodeIDBuffer.emplace_back(edge2.end, edge2.length);
			}
		}
	}
}

void PairedDBG::divideNestedBubbleNode(const long numThread)
{
    cerr << "dividing nested bubble scaffolds..." << endl;

	setOppositeBubbleNodeIDForEachNode(numThread);

	const long MIN_OVERLAP_TO_JOIN = 20;
    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

	vector<char> nestedBubbleFlag(this->node.size(), false);

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (node[nodeIndex].state & SC_DEL)
			continue;

		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0) {
			long oppositeNodeIndex = id2Index(oppositeNodeID);
			if (node[oppositeNodeIndex].oppositeBubbleNodeID != 0 && id2Index(node[oppositeNodeIndex].oppositeBubbleNodeID) != nodeIndex)
				nestedBubbleFlag[nodeIndex] = nestedBubbleFlag[oppositeNodeIndex] = true;
		}
	}


    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		if (node[nodeIndex].state & SC_DEL)
			continue;

		long oppositeNodeID = node[nodeIndex].oppositeBubbleNodeID;
		if (oppositeNodeID != 0 && nestedBubbleFlag[id2Index(oppositeNodeID)])
			nestedBubbleFlag[nodeIndex] =  true;
	}


    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		long i;
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.back() = 1;

		if (nestedBubbleFlag[nodeIndex]) {
			for (i = 1; i < static_cast<long>(contigRef.size()); ++i) {
				if ((this->contigState[id2Index(contigRef[i].id)] & DBG_CONTIG_BUBBLE_JUNCTION) || contigRef[i].start > contigRef[i - 1].end || this->getOverlap(contigRef[i - 1].id, contigRef[i].id) < minOverlap)
					breakpointFlag[i] = 1;
			}
		}

		i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->minOverlap = defaultMinOverlap;
    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::deleteLongEdge(const long maxEdgeLength)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing long-gap edges (GAP_LENGTH > " <<  maxEdgeLength << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
        for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
            if (node1.edge[edgeID].length > maxEdgeLength) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
                ++numDelete;
            }
        }
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::deleteErroneousEdgeNumTagRate(const long numThread)
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

				long numTag1 = this->getCommonTagBetweenNodePair(nodeID + 1, edge1.end);
				long numTag2 = this->getCommonTagBetweenNodePair(nodeID + 1, edge2.end);

                if (numTag1 < RATE_THRESHOLD * numTag2) {
                    #pragma omp critical (delete_edge)
                    {
                        ids.push_back(nodeID + 1);
                        ids.push_back(edge1.end);
                    }
                    ++numDelete;
                } else if (numTag2 < RATE_THRESHOLD * numTag1) {
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

void PairedDBG::deleteErroneousEdgeNumTagRateIterative(const long numThread)
{
	if (this->tagLibraryMT == NULL)
		return;

	countMappedTagForEachScaffold(numThread);

    cerr << "removing erroneous edges using linked-reads (barcodes)..." << endl << endl;;
    long totalDelete = 0;
    long numDelete;
    do {
        numDelete = this->deleteErroneousEdgeNumTagRate(numThread);
        totalDelete += numDelete;
        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
    } while (numDelete > 0);

    cerr << "TOTAL_NUM_REMOVED_EDGES_BY_BARCODE =" << totalDelete << endl << endl;
}


//added by ouchi
long PairedDBG::deleteErroneousEdgebyHiC(const long numThread)
{
    vector<long> ids;
    long numDelete = 0;

    omp_set_num_threads(numThread);

    //#pragma omp parallel for schedule(dynamic) reduction(+: numDelete)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        for (long edgeID = 0; edgeID < node[nodeID].numEdge; ++edgeID) {
            if (node[nodeID].length < 200000 || node[id2Index(node[nodeID].edge[edgeID].end)].length < 200000) continue;
            if (this->checkHiCLinkBetweenNodePair(nodeID+1, node[nodeID].edge[edgeID].direction, node[nodeID].edge[edgeID].end)) continue;
            //#pragma omp critical (delete_edge)
            {
//				std::cerr << "nodeID:" << nodeID << "\tnodeID2:" << node[nodeID].edge[edgeID].end << std::endl;
                ids.push_back(nodeID + 1);
                ids.push_back(node[nodeID].edge[edgeID].end);
            }
            ++numDelete;
        }
    }
    this->deleteEdges(ids);
    return numDelete;
}
//

//added by ouchi
void PairedDBG::deleteErroneousEdgebyHiCIterative(const long numThread)
{
	if (this->HiCLibraryMT == NULL)
		return;

    cerr << "removing erroneous edges using HiC reads ..." << endl << endl;;

    makeHiCLink(numThread);

    long totalDelete = 0;
    long numDelete;
//    do {
        numDelete = this->deleteErroneousEdgebyHiC(numThread);
        totalDelete += numDelete;
//        cerr << "NUM_REMOVED_EDGES =" << numDelete << endl;
//    } while (numDelete > 0);

    cerr << "TOTAL_NUM_REMOVED_EDGES_BY_HIC =" << totalDelete << endl << endl;
}
//

void PairedDBG::deleteEdgeFromShortNodeKeepingBubble(const long lengthThreshold)
{
    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from short nodes (LENGTH < " <<  lengthThreshold << ") ..." << std::endl;

    for (long nodeID = 0 ; nodeID < numNode; ++nodeID) {
        GraphNode &node1 = this->node[nodeID];
		if (node1.length < lengthThreshold && node1.oppositeBubbleNodeID == 0) {
			for (long edgeID = 0; edgeID < node1.numEdge; ++edgeID) {
				ids.push_back(nodeID + 1);
				ids.push_back(node1.edge[edgeID].end);
				++numDelete;
			}
		}
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

long PairedDBG::divideInconsistentBubbleEnd()
{
    cerr << "dividing inconsistent bubble-ends..." << endl;

	this->setOppositeBubbleNodeIDAndStateForEachNode();

    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

	vector<vector<char> > breakpointFlag(this->node.size());
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		breakpointFlag[nodeIndex].assign(this->node[nodeIndex].contig.size() + 1, 0);
		breakpointFlag[nodeIndex].front() = 1;
		breakpointFlag[nodeIndex].back() = 1;
	}

	vector<char> leftEndFlag(this->numContig, 0);
	vector<char> rightEndFlag(this->numContig, 0);

	long totalDivision = -1;
	long numDivision = 1;
	while (numDivision > 0) {
		totalDivision += numDivision;
		numDivision = 0;
		for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
			if (!(node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
				continue;

			vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
			for (unsigned long i = 0; i < contigRef.size(); ++i) {
				if (breakpointFlag[nodeIndex][i]) {
					if (contigRef[i].id > 0)
						leftEndFlag[id2Index(contigRef[i].id)] = 1;
					else
						rightEndFlag[id2Index(contigRef[i].id)] = 1;
				}

				if (breakpointFlag[nodeIndex][i + 1]) {
					if (contigRef[i].id > 0)
						rightEndFlag[id2Index(contigRef[i].id)] = 1;
					else
						leftEndFlag[id2Index(contigRef[i].id)] = 1;
				}
			}
		}

		for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
			if (!(node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
				continue;

			vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
			for (unsigned long i = 0; i < contigRef.size(); ++i) {
				long j = contigRef[i].id > 0 ? 0 : 1;
				long oppositeID = sign(contigRef[i].id) *  this->contigBubbleInfo[id2Index(contigRef[i].id)].oppositeContigID[j];
				long oppositeNodeID = 0; //edited by ouchi (contigPositionInScaffold)
				if (oppositeID != 0) { //added by ouchi
					if (contigPositionInScaffold[id2Index(oppositeID)].size() == 1) //edited by ouchi (contigPositionInScaffold)
						oppositeNodeID = contigPositionInScaffold[id2Index(oppositeID)][0].id; //edited by ouchi (contigPositionInScaffold)
				} //added by ouchi
				if (!(oppositeID == 0 || (abs(oppositeNodeID) != (nodeIndex + 1) && abs(oppositeNodeID) != abs(node[nodeIndex].oppositeBubbleNodeID)))) { //edited by ouchi (contigPositionInScaffold)
					if (oppositeID > 0) {
						if (leftEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i]) {
							breakpointFlag[nodeIndex][i] = 1;
							++numDivision;
							break;
						}
					}
					else {
						if (rightEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i]) {
							breakpointFlag[nodeIndex][i] = 1;
							++numDivision;
							break;
						}
					}
				}

				j = contigRef[i].id > 0 ? 1 : 0;
				oppositeID = sign(contigRef[i].id) *  this->contigBubbleInfo[id2Index(contigRef[i].id)].oppositeContigID[j];
				oppositeNodeID = 0; //edited by ouchi (contigPositionInScaffold)
				if (oppositeID != 0) { //added by ouchi
					if (contigPositionInScaffold[id2Index(oppositeID)].size() == 1) //edited by ouchi (contigPositionInScaffold)
						oppositeNodeID = contigPositionInScaffold[id2Index(oppositeID)][0].id; //edited by ouchi (contigPositionInScaffold)
				} //added by ouchi
				if (!(oppositeID == 0 || (abs(oppositeNodeID) != (nodeIndex + 1) && abs(oppositeNodeID) != abs(node[nodeIndex].oppositeBubbleNodeID)))) { //edited by ouchi (contigPositionInScaffold)
					if (oppositeID > 0) {
						if (rightEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i + 1]) {
							breakpointFlag[nodeIndex][i + 1] = 1;
							++numDivision;
							break;
						}
					}
					else {
						if (leftEndFlag[id2Index(oppositeID)] && !breakpointFlag[nodeIndex][i + 1]) {
							breakpointFlag[nodeIndex][i + 1] = 1;
							++numDivision;
							break;
						}
					}
				}
			}
		}
	}


	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[nodeIndex][i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;


			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

	return totalDivision;
}

void PairedDBG::adjustOppositeBubbleNodeIDDirection()
{
	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		if (!(node[nodeIndex].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))
			continue;

		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;
		for (long j = 0; j < 2; ++j) {
			long contigID;
			long oppositeContigID;
			if (j == 0) {
				contigID = contigRef.front().id;
				if (contigID > 0)
					oppositeContigID = this->contigBubbleInfo[id2Index(contigID)].oppositeContigID[0];
				else
					oppositeContigID = this->contigBubbleInfo[id2Index(contigID)].oppositeContigID[1];
			}
			else {
				contigID = contigRef.back().id;
				if (contigID > 0)
					oppositeContigID = this->contigBubbleInfo[id2Index(contigID)].oppositeContigID[1];
				else
					oppositeContigID = this->contigBubbleInfo[id2Index(contigID)].oppositeContigID[0];
			}

			if (oppositeContigID != 0) {
				long oppositeNodeID = 0;
				if (contigPositionInScaffold[id2Index(oppositeContigID)].size() == 1)
					oppositeNodeID = this->contigPositionInScaffold[id2Index(oppositeContigID)][0].id; //edited by ouchi (contigPositionInScaffold)
				node[nodeIndex].oppositeBubbleNodeID = sign(contigID) * sign(oppositeContigID) * sign(oppositeNodeID) * abs(node[nodeIndex].oppositeBubbleNodeID);
				break;
			}
		}
	}
}

void PairedDBG::deleteEdgeFromSecondaryBubble()
{
	vector<char> secondaryBubbleFlag(numNode, 0);
    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &node1 = this->node[nodeIndex];

		for (long contigIndex = 0; contigIndex < node1.numContig; ++contigIndex) {
			if (contigState[id2Index(node1.contig[contigIndex].id)] & DBG_CONTIG_SECONDARY_BUBBLE) {
				secondaryBubbleFlag[nodeIndex] = 1;
				break;
			}
		}
	}

    vector<long> ids;
    unsigned long long numDelete = 0;

    std::cerr << "removing edges from bubble nodes ..." << std::endl;

    for (long nodeIndex = 0 ; nodeIndex < numNode; ++nodeIndex) {
        GraphNode &node1 = this->node[nodeIndex];

		for (long edgeIndex = 0; edgeIndex < node1.numEdge; ++edgeIndex) {
//			if (node1.oppositeBubbleNodeID == 0 && node[id2Index(node1.edge[edgeIndex].end)].oppositeBubbleNodeID == 0)
			if (!(secondaryBubbleFlag[nodeIndex]) && !(secondaryBubbleFlag[id2Index(node1.edge[edgeIndex].end)]))
				continue;

			ids.push_back(nodeIndex + 1);
			ids.push_back(node1.edge[edgeIndex].end);
			++numDelete;
		}
    }
    std::cerr << "TOTAL_NUM_DELETE=" << numDelete << std::endl;
    this->deleteEdges(ids);
}

void PairedDBG::divideNodeBasedOnBubblesIterative(const bool strandFlag, const long numThread)
{
	unsigned currentMode = getMode();
	setMode(0);

	long total = 0;
	long num;

	cerr << endl << "dividing nodes based on bubbles ..." << endl;
	do {
		num = divideNodeUsingBubbleContigPair(numThread);
		num += divideInconsistentBubbleEnd();
		if (strandFlag)
			num += divideNodeUsingBubbleContigPairStrandAware(numThread);

		total += num;
		cerr << "NUM_DIVISION = " << num << endl;
	} while (num > 0);

	cerr << "TOTAL_NUM_DIVISIONS =" << total << endl << endl;
	setMode(currentMode);
}

void PairedDBG::divideNodeByPrimaryBubbleBoundary(const PairedDBG &bubbleGraph)
{
    cerr << "dividing scaffolds by bubble-boundaries ..." << endl;

	std::vector<std::vector<long> > bubbleContigID;
    for (long nodeIndex = 0; nodeIndex < bubbleGraph.numNode; ++nodeIndex) {
		if (!(bubbleGraph.node[nodeIndex].state & DBG_PRIMARY_BUBBLE))
			continue;

		bubbleContigID.resize(bubbleContigID.size() + 1);
		bubbleContigID.back().resize(bubbleGraph.node[nodeIndex].contig.size());
		for (unsigned i = 0; i < bubbleGraph.node[nodeIndex].contig.size(); ++i)
			bubbleContigID.back()[i] = bubbleGraph.node[nodeIndex].contig[i].id;
	}

	vector<vector<long> > bubbleHeadIndex(this->numContig);
	for (unsigned long i = 0; i < bubbleContigID.size(); ++i)
		bubbleHeadIndex[id2Index(bubbleContigID[i][0])].push_back(i);


    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
	vector<long> numContigUsed(this->numContig, 0);

	for (unsigned long nodeIndex = 0; nodeIndex < this->node.size(); ++nodeIndex) {
		vector<ScaffoldPart> &contigRef = this->node[nodeIndex].contig;

		vector<char> breakpointFlag(contigRef.size() + 1, 0);
		breakpointFlag.front() = 1;
		breakpointFlag.back() = 1;

		for (unsigned long i = 0; i < contigRef.size(); ++i) {
			long contigIndex = id2Index(contigRef[i].id);
			if (bubbleHeadIndex[contigIndex].empty())
				continue;

			for (unsigned long j = 0; j < bubbleHeadIndex[contigIndex].size(); ++j) {
				std::vector<long> &bubbleContigIDRef = bubbleContigID[bubbleHeadIndex[contigIndex][j]];
				unsigned long k;
				if (contigRef[i].id == bubbleContigIDRef[0]) {
					for (k = 1; k < bubbleContigIDRef.size(); ++k) {
						if (contigRef[i + k].id != bubbleContigIDRef[k])
							break;
					}
					if (k == bubbleContigIDRef.size()) {
						breakpointFlag[i] = 1;
						breakpointFlag[i + k] = 1;
					}
				}
				else {
					for (k = 1; k < bubbleContigIDRef.size(); ++k) {
						if (contigRef[i - k].id != -(bubbleContigIDRef[k]))
							break;
					}
					if (k == bubbleContigIDRef.size()) {
						breakpointFlag[i + 1] = 1;
						breakpointFlag[i - k + 1] = 1;
					}
				}
			}
		}

		long i = 0;
		while (i < static_cast<long>(contigRef.size())) {
			long start = contigRef[i].start;
			long j = i;
			while (breakpointFlag[i + 1] == 0) {
				contigRef[i].start -= start;
				contigRef[i].end -= start;
				++i;
			}
			contigRef[i].start -= start;
			contigRef[i].end -= start;
			++i;

			bool uniqueFlag = false;
			for (long k = j; k < i; ++k) {
				if (numContigUsed[id2Index(contigRef[k].id)] == 0)
					uniqueFlag = true;
			}

			if (uniqueFlag) {
				long k = i - j;
				++numNewNode;
				newContigPoolSize += k;
				fwrite(&k, sizeof(long), 1, scaffoldFP);
				for (k = j; k < i; ++k) {
					fwrite(&(contigRef[k]), sizeof(ScaffoldPart), 1, scaffoldFP);
					++numContigUsed[id2Index(contigRef[k].id)];
				}
			}
		}
	}

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::remakeGraphRecoveringSecondaryBubble(PairedDBG &bubbleGraph)
{
    long numNewNode = 0;
    long newContigPoolSize = 0;
    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long i = 0; i < numNode; ++i) {
        if (node[i].state & SC_DEL)
			continue;

		fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
		for (long j = 0; j < node[i].numContig; ++j)
			fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
		newContigPoolSize += node[i].numContig;
		++numNewNode;
    }

    for (long i = 0; i < bubbleGraph.numNode; ++i) {
        if ((bubbleGraph.node[i].state & SC_DEL) || !(bubbleGraph.node[i].state & DBG_SECONDARY_BUBBLE))
			continue;

		fwrite(&(bubbleGraph.node[i].numContig), sizeof(long), 1, scaffoldFP);
		for (long j = 0; j < bubbleGraph.node[i].numContig; ++j)
			fwrite(&(bubbleGraph.node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
		newContigPoolSize += bubbleGraph.node[i].numContig;
		++numNewNode;
    }

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);
}

void PairedDBG::divideErroneousLink(const vector<vector<unsigned> >& numErroneousPair, const vector<vector<unsigned> >& numSpanningPair, const vector<vector<double> >& sumExpectedLink, std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink, const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize)
{
	const long MIN_OVERLAP_TO_JOIN = 32;

    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		this->contigPreviousParentNodeID.resize(this->contig.size(), 0);

    unsigned long numDivided = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    bool isSplit = false;
    FILE *scaffoldFP = platanus::makeTemporaryFile();
    vector<char> breakpoint;
    for (long i = 0; i < numNode; ++i) {
        if (node[i].numContig == 1) {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            fwrite(&(node[i].contig[0]), sizeof(ScaffoldPart), 1, scaffoldFP);
            ++numNewNode;
            ++newContigPoolSize;
            continue;
        }

        isSplit = false;
        breakpoint.resize(node[i].numContig + 1, 0);
        std::fill(breakpoint.begin(), breakpoint.end(), 0);
        breakpoint.back() = 1;
        for (unsigned long j = 0; j < numSpanningPair[i].size(); ++j) {
			if (divisionMode == SWITCH) {
//				if (numErroneousPair[i][j] <= numSpanningPair[i][j])
				if (numErroneousPair[i][j] <= numSpanningPair[i][j] || numErroneousPair[i][j] < minLink)
					continue;
			}
			else if (divisionMode == GAP) {
				if (node[i].contig[j].end - node[i].contig[j + 1].start < maxGapSize || numErroneousPair[i][j] < minLink)
					continue;
			}
			else {
//				long overlapLen = this->getOverlap(node[i].contig[j].id, node[i].contig[j + 1].id);
//				if (!(!(overlapLen >= minOverlap && node[i].contig[j].end - node[i].contig[j + 1].start <= 0) && numErroneousPair[i][j] >= minLink &&  numSpanningPair[i][j] < minLink && sumExpectedLink[i][j] > 1.0 && (double)numSpanningPair[i][j] < sumExpectedLink[i][j] * CHECK_USING_LONGER_LIB_TH))
				if (numErroneousPair[i][j] <= numSpanningPair[i][j])
					continue;
			}

			++numDivided;
			breakpoint[j + 1] = 1;
			isSplit = true;
        }
        if (isSplit) {
            long j = 0;
            while (j < node[i].numContig) {
                long start = node[i].contig[j].start;
                long k = j;
                while (breakpoint[j + 1] == 0) {
                    node[i].contig[j].start -= start;
                    node[i].contig[j].end -= start;
                    ++j;
                }
                node[i].contig[j].start -= start;
                node[i].contig[j].end -= start;
                ++j;
                long tmp = j - k;
                fwrite(&tmp, sizeof(long), 1, scaffoldFP);
                for (tmp = k; tmp < j; ++tmp)
                    fwrite(&(node[i].contig[tmp]), sizeof(ScaffoldPart), 1, scaffoldFP);
                ++numNewNode;
                newContigPoolSize += tmp;

				for (auto itr = node[i].contig.begin(); itr != node[i].contig.end(); ++itr) {
					if (this->contigPositionInScaffold[id2Index(itr->id)].size() == 1) //edited by ouchi (contigPositionInScaffold)
						this->contigPreviousParentNodeID[id2Index(itr->id)] = i + 1;
					else
						this->contigPreviousParentNodeID[id2Index(itr->id)] = 0;
				}
            }
        }
        else {
            fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
            auto contigIterator = node[i].contig.begin();
            auto contigEnd = node[i].contig.end();
            for (; contigIterator != contigEnd; ++contigIterator)
                fwrite(&(*contigIterator), sizeof(ScaffoldPart), 1, scaffoldFP);
            ++numNewNode;
            newContigPoolSize += node[i].numContig;
        }
    }

    std::cerr << "NUM_DIVIDED_SWITCH_ERROR_CANDIDATES = " << numDivided << endl;

    this->remake(numNewNode, newContigPoolSize, scaffoldFP);
    fclose(scaffoldFP);

    this->minOverlap = defaultMinOverlap;
}

void PairedDBG::countPairsSpanningGap(vector<vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, long numThread)
{
    const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numSpanningPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numSpanningPairThread[threadID].resize(numSpanningPair.size());
		for (unsigned i = 0; i < numSpanningPair.size(); ++i)
			numSpanningPairThread[threadID][i].resize(numSpanningPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		unsigned numPairLink;
		unsigned numReadLink;

		rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				unsigned long forwardIndex = id2Index(forwardResult.id);
				unsigned long reverseIndex = id2Index(reverseResult.id);

				if ((this->contigBubbleInfo[forwardIndex].oppositeContigID[0] != this->contigBubbleInfo[forwardIndex].oppositeContigID[1]) ||
					(this->contigBubbleInfo[reverseIndex].oppositeContigID[0] != this->contigBubbleInfo[reverseIndex].oppositeContigID[1])) {
					continue;
				}

				if (divisionMode == SWITCH && (this->contigBubbleInfo[forwardIndex].oppositeContigID[0] == 0 || this->contigBubbleInfo[reverseIndex].oppositeContigID[0] == 0)) {
					continue;
				}

				if (contigPositionInScaffold[forwardIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex][0].id : -(contigPositionInScaffold[forwardIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset = contigPositionInScaffold[forwardIndex][0].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				if (contigPositionInScaffold[reverseIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseIndex][0].id : -(contigPositionInScaffold[reverseIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset = contigPositionInScaffold[reverseIndex][0].id > 0 ? reverseResult.offset : contig[reverseIndex].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset += node[abs(reverseResult.id) - 1].contig[contigPositionInScaffold[reverseIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				if (forwardResult.id != -reverseResult.id || abs(forwardResult.offset - reverseResult.offset) - averageInsSize > tolerence)
					continue;

				unsigned long leftMostContigNum = std::min(contigPositionInScaffold[forwardIndex][0].offset, contigPositionInScaffold[reverseIndex][0].offset); //edited by ouchi (contigPositionInScaffold)
				unsigned long rightMostContigNum = std::max(contigPositionInScaffold[forwardIndex][0].offset, contigPositionInScaffold[reverseIndex][0].offset); //edited by ouchi (contigPositionInScaffold)

				for (unsigned long i = leftMostContigNum; i < rightMostContigNum; ++i) {
					std::pair<int, int> redundancyCheckKey = std::make_pair(abs(forwardResult.id) - 1, i);
					if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
						redundancyCheckSet.insert(redundancyCheckKey);
						++(numSpanningPairThread[threadID][abs(forwardResult.id) - 1][i]);
					}
				}
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
			}
        }
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numSpanningPair[i].size(); ++j) {
				numSpanningPair[i][j] += numSpanningPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLongReadSpanningGap(vector<vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, long numThread)
{
    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numSpanningPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numSpanningPairThread[threadID].resize(numSpanningPair.size());
		for (unsigned i = 0; i < numSpanningPair.size(); ++i)
			numSpanningPairThread[threadID][i].resize(numSpanningPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numResult;
		long score;
		vector<platanus::Position> positionBuffer;

        rewind((*longReadLibraryMT)[threadID].mappedReadFP);
		while (fread(&numResult, sizeof(unsigned), 1, (*longReadLibraryMT)[threadID].mappedReadFP)) {
			positionBuffer.resize(numResult);
			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&(positionBuffer[i]), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fseek((*longReadLibraryMT)[threadID].mappedReadFP, 2 * sizeof(int), SEEK_CUR);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
			}

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (auto positionItr1 = positionBuffer.begin(); positionItr1 != positionBuffer.end() - 1; ++positionItr1) {
				for (auto positionItr2 = positionItr1 + 1; positionItr2 != positionBuffer.end(); ++positionItr2) {
					unsigned long forwardIndex = id2Index(positionItr1->id);
					unsigned long reverseIndex = id2Index(positionItr2->id);
					if ((this->contigBubbleInfo[forwardIndex].oppositeContigID[0] != this->contigBubbleInfo[forwardIndex].oppositeContigID[1]) ||
						(this->contigBubbleInfo[reverseIndex].oppositeContigID[0] != this->contigBubbleInfo[reverseIndex].oppositeContigID[1])) {
						continue;
					}

					if (divisionMode == SWITCH && (this->contigBubbleInfo[forwardIndex].oppositeContigID[0] == 0 || this->contigBubbleInfo[reverseIndex].oppositeContigID[0] == 0)) {
						continue;
					}

					if (contigPositionInScaffold[forwardIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
					positionItr1->id = positionItr1->id > 0 ? contigPositionInScaffold[forwardIndex][0].id : -(contigPositionInScaffold[forwardIndex][0].id); //edited by ouchi (contigPositionInScaffold)
					positionItr1->offset = contigPositionInScaffold[forwardIndex][0].id > 0 ? positionItr1->offset : contig[forwardIndex].length - positionItr1->offset - 1; //edited by ouchi (contigPositionInScaffold)
					positionItr1->offset += node[abs(positionItr1->id) - 1].contig[contigPositionInScaffold[forwardIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

					if (contigPositionInScaffold[reverseIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
					positionItr2->id = positionItr2->id > 0 ? contigPositionInScaffold[reverseIndex][0].id : -(contigPositionInScaffold[reverseIndex][0].id); //edited by ouchi (contigPositionInScaffold)
					positionItr2->offset = contigPositionInScaffold[reverseIndex][0].id > 0 ? positionItr2->offset : contig[reverseIndex].length - positionItr2->offset - 1; //edited by ouchi (contigPositionInScaffold)
					positionItr2->offset += node[abs(positionItr2->id) - 1].contig[contigPositionInScaffold[reverseIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

					if (positionItr1->id != positionItr2->id) continue;

					unsigned long leftMostContigNum = std::min(contigPositionInScaffold[forwardIndex][0].offset, contigPositionInScaffold[reverseIndex][0].offset); //edited by ouchi (contigPositionInScaffold)
					unsigned long rightMostContigNum = std::max(contigPositionInScaffold[forwardIndex][0].offset, contigPositionInScaffold[reverseIndex][0].offset); //edited by ouchi (contigPositionInScaffold)

					for (unsigned long i = leftMostContigNum; i < rightMostContigNum; ++i) {
						std::pair<int, int> redundancyCheckKey = std::make_pair(abs(positionItr1->id) - 1, i);
						if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
							redundancyCheckSet.insert(redundancyCheckKey);
							++(numSpanningPairThread[threadID][abs(positionItr1->id) - 1][i]);
						}
					}
				}
			}
		}
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numSpanningPair[i].size(); ++j) {
				numSpanningPair[i][j] += numSpanningPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLinksInsideContigs(vector<vector<unsigned> >& numPair, const long numThread)
{
	const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position contigResultF;
		platanus::Position contigResultR;
		platanus::Position scaffoldResultF;
		platanus::Position scaffoldResultR;
		long scaffoldOverhangF;
		long scaffoldOverhangR;
		long contigOverhangF;
		long contigOverhangR;
		unsigned long contigIndexF;
		unsigned long contigIndexR;
		unsigned numPairLink;
		unsigned numReadLink;

		FILE *mappedFP = (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP;
        rewind(mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, mappedFP);

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&contigResultF, sizeof(platanus::Position), 1, mappedFP);
				fread(&contigResultR, sizeof(platanus::Position), 1, mappedFP);

				contigIndexF = abs(contigResultF.id) - 1;
				if (contigPositionInScaffold[contigIndexF].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				scaffoldResultF.id = contigResultF.id > 0 ? contigPositionInScaffold[contigIndexF][0].id : -(contigPositionInScaffold[contigIndexF][0].id); //edited by ouchi (contigPositionInScaffold)
				scaffoldResultF.offset = contigPositionInScaffold[contigIndexF][0].id > 0 ? contigResultF.offset : contig[contigIndexF].length - contigResultF.offset - 1; //edited by ouchi (contigPositionInScaffold)
				scaffoldResultF.offset += node[abs(scaffoldResultF.id) - 1].contig[contigPositionInScaffold[contigIndexF][0].offset].start; //edited by ouchi (contigPositionInScaffold)
				scaffoldOverhangF = scaffoldResultF.id > 0  ? node[scaffoldResultF.id - 1].length - scaffoldResultF.offset : scaffoldResultF.offset;

				contigIndexR = abs(contigResultR.id) - 1;
				if (contigPositionInScaffold[contigIndexR].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				scaffoldResultR.id = contigResultR.id > 0 ? contigPositionInScaffold[contigIndexR][0].id : -(contigPositionInScaffold[contigIndexR][0].id); //edited by ouchi (contigPositionInScaffold)
				scaffoldResultR.offset = contigPositionInScaffold[contigIndexR][0].id > 0 ? contigResultR.offset : contig[contigIndexR].length - contigResultR.offset - 1; //edited by ouchi (contigPositionInScaffold)
				scaffoldResultR.offset += node[abs(scaffoldResultR.id) - 1].contig[contigPositionInScaffold[contigIndexR][0].offset].start; //edited by ouchi (contigPositionInScaffold)
				scaffoldOverhangR = scaffoldResultR.id > 0  ? node[scaffoldResultR.id - 1].length - scaffoldResultR.offset : scaffoldResultR.offset;

				if (scaffoldResultF.id == -scaffoldResultR.id || scaffoldOverhangF + scaffoldOverhangR <= averageInsSize + tolerence)
					continue;

				contigOverhangF = contigResultF.id > 0  ? contig[contigIndexF].length - contigResultF.offset : contigResultF.offset;
				if (contigOverhangF > averageInsSize + tolerence) {
					if (scaffoldResultF.id > 0) {
						for (long i = contigPositionInScaffold[contigIndexF][0].offset; i < node[scaffoldResultF.id - 1].numContig - 1; ++i) { //edited by ouchi (contigPositionInScaffold)
							if (node[scaffoldResultF.id - 1].contig[i].end - scaffoldResultF.offset <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultF.id - 1, i);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][scaffoldResultF.id - 1][i]);
								}
							}
							else {
								break;
							}
						}
					}
					else {
						for (long i = contigPositionInScaffold[contigIndexF][0].offset; i > 0; --i) { //edited by ouchi (contigPositionInScaffold)
							if (scaffoldResultF.offset - node[-(scaffoldResultF.id) - 1].contig[i].start <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultF.id) - 1, i - 1);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][-(scaffoldResultF.id) - 1][i - 1]);
								}
							}
							else {
								break;
							}
						}
					}
				}

				contigOverhangR = contigResultR.id > 0  ? contig[contigIndexR].length - contigResultR.offset : contigResultR.offset;
				if (contigOverhangR > averageInsSize + tolerence) {
					if (scaffoldResultR.id > 0) {
						for (long i = contigPositionInScaffold[contigIndexR][0].offset; i < node[scaffoldResultR.id - 1].numContig - 1; ++i) { //edited by ouchi (contigPositionInScaffold)
							if (node[scaffoldResultR.id - 1].contig[i].end - scaffoldResultR.offset <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultR.id - 1, i);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][scaffoldResultR.id - 1][i]);
								}
							}
							else {
								break;
							}
						}
					}
					else {
						for (long i = contigPositionInScaffold[contigIndexR][0].offset; i > 0; --i) { //edited by ouchi (contigPositionInScaffold)
							if (scaffoldResultR.offset - node[-(scaffoldResultR.id) - 1].contig[i].start <= averageInsSize) {
								std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultR.id) - 1, i - 1);
								if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
									redundancyCheckSet.insert(redundancyCheckKey);
									++(numPairThread[threadID][-(scaffoldResultR.id) - 1][i - 1]);
								}
							}
							else {
								break;
							}
						}
					}
				}
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&contigResultF, sizeof(platanus::Position), 1, mappedFP);
				fread(&contigResultR, sizeof(platanus::Position), 1, mappedFP);
			}
		}
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLongReadLinksInsideContigs(vector<vector<unsigned> >& numPair, const long numThread)
{
	const long averageInsSize = (*longReadLibraryMT)[0].getAverageInsSize();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position scaffoldResultF;
		platanus::Position scaffoldResultR;
		long scaffoldOverhangF;
		long scaffoldOverhangR;
		long contigOverhangF;
		long contigOverhangR;
		unsigned long contigIndexF;
		unsigned long contigIndexR;

		FILE *mappedFP = (*longReadLibraryMT)[threadID].mappedReadFP;
		unsigned numResult;
		long score;
		vector<platanus::Position> positionBuffer;

        rewind(mappedFP);
		while (fread(&numResult, sizeof(unsigned), 1, mappedFP)) {
			positionBuffer.resize(numResult);
			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&(positionBuffer[i]), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fseek((*longReadLibraryMT)[threadID].mappedReadFP, 2 * sizeof(int), SEEK_CUR);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
			}

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (auto positionItr1 = positionBuffer.begin(); positionItr1 != positionBuffer.end() - 1; ++positionItr1) {
				for (auto positionItr2 = positionItr1 + 1; positionItr2 != positionBuffer.end(); ++positionItr2) {
					contigIndexF = abs(positionItr1->id) - 1;
					if (contigPositionInScaffold[contigIndexF].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
					scaffoldResultF.id = positionItr1->id > 0 ? contigPositionInScaffold[contigIndexF][0].id : -(contigPositionInScaffold[contigIndexF][0].id); //edited by ouchi (contigPositionInScaffold)
					scaffoldResultF.offset = contigPositionInScaffold[contigIndexF][0].id > 0 ? positionItr1->offset : contig[contigIndexF].length - positionItr1->offset - 1; //edited by ouchi (contigPositionInScaffold)
					scaffoldResultF.offset += node[abs(scaffoldResultF.id) - 1].contig[contigPositionInScaffold[contigIndexF][0].offset].start; //edited by ouchi (contigPositionInScaffold)
					scaffoldOverhangF = scaffoldResultF.id > 0  ? node[scaffoldResultF.id - 1].length - scaffoldResultF.offset : scaffoldResultF.offset;

					contigIndexR = abs(positionItr2->id) - 1;
					if (contigPositionInScaffold[contigIndexR].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
					scaffoldResultR.id = positionItr2->id > 0 ? contigPositionInScaffold[contigIndexR][0].id : -(contigPositionInScaffold[contigIndexR][0].id); //edited by ouchi (contigPositionInScaffold)
					scaffoldResultR.offset = contigPositionInScaffold[contigIndexR][0].id > 0 ? positionItr2->offset : contig[contigIndexR].length - positionItr2->offset - 1; //edited by ouchi (contigPositionInScaffold)
					scaffoldResultR.offset += node[abs(scaffoldResultR.id) - 1].contig[contigPositionInScaffold[contigIndexR][0].offset].start; //edited by ouchi (contigPositionInScaffold)
					scaffoldOverhangR = scaffoldResultR.id > 0  ? node[scaffoldResultR.id - 1].length - scaffoldResultR.offset : scaffoldResultR.offset;

					if (scaffoldResultF.id == scaffoldResultR.id || scaffoldOverhangF + scaffoldOverhangR <= averageInsSize + tolerence) continue;

					contigOverhangF = positionItr1->id > 0  ? contig[contigIndexF].length - positionItr1->offset : positionItr1->offset;
					if (contigOverhangF > averageInsSize + tolerence) {
						if (scaffoldResultF.id > 0) {
							for (long i = contigPositionInScaffold[contigIndexF][0].offset; i < node[scaffoldResultF.id - 1].numContig - 1; ++i) { //edited by ouchi (contigPositionInScaffold)
								if (node[scaffoldResultF.id - 1].contig[i].end - scaffoldResultF.offset <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultF.id - 1, i);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										redundancyCheckSet.insert(redundancyCheckKey);
										++(numPairThread[threadID][scaffoldResultF.id - 1][i]);
									}
								}
								else {
									break;
								}
							}
						}
						else {
							for (long i = contigPositionInScaffold[contigIndexF][0].offset; i > 0; --i) { //edited by ouchi (contigPositionInScaffold)
								if (scaffoldResultF.offset - node[-(scaffoldResultF.id) - 1].contig[i].start <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultF.id) - 1, i - 1);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										redundancyCheckSet.insert(redundancyCheckKey);
										++(numPairThread[threadID][-(scaffoldResultF.id) - 1][i - 1]);
									}
								}
								else {
									break;
								}
							}
						}
					}

					contigOverhangR = positionItr2->id > 0  ? contig[contigIndexR].length - positionItr2->offset : positionItr2->offset;
					if (contigOverhangR > averageInsSize + tolerence) {
						if (scaffoldResultR.id > 0) {
							for (long i = contigPositionInScaffold[contigIndexR][0].offset; i < node[scaffoldResultR.id - 1].numContig - 1; ++i) { //edited by ouchi (contigPositionInScaffold)
								if (node[scaffoldResultR.id - 1].contig[i].end - scaffoldResultR.offset <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldResultR.id - 1, i - 1);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										redundancyCheckSet.insert(redundancyCheckKey);
										++(numPairThread[threadID][scaffoldResultR.id - 1][i]);
									}
								}
								else {
									break;
								}
							}
						}
						else {
							for (long i = contigPositionInScaffold[contigIndexR][0].offset; i > 0; --i) { //edited by ouchi (contigPositionInScaffold)
								if (scaffoldResultR.offset - node[-(scaffoldResultR.id) - 1].contig[i].start <= averageInsSize) {
									std::pair<int, int> redundancyCheckKey = std::make_pair(-(scaffoldResultR.id) - 1, i - 1);
									if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
										++(numPairThread[threadID][-(scaffoldResultR.id) - 1][i - 1]);
									}
								}
								else {
									break;
								}
							}
						}
					}
				}
			}
        }
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countSwitchErrorLinks(vector<vector<unsigned> >& numPair, const long numThread)
{
    const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
	setOppositeBubbleNodeIDAndStateForEachNode();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		std::array<std::array<platanus::Position, 2>, 2> contigResult;
		std::array<std::array<platanus::Position, 2>, 2> scaffoldResult;
		std::array<std::array<unsigned long, 2>, 2> contigIndex;
		std::array<std::array<unsigned long, 2>, 2> scaffoldIndex;
		unsigned numPairLink;
		unsigned numReadLink;
		unsigned i, j;

		FILE *mappedFP;
		mappedFP = (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP;

        rewind(mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, mappedFP);

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&(contigResult[0][0]), sizeof(platanus::Position), 1, mappedFP);
				fread(&(contigResult[0][1]), sizeof(platanus::Position), 1, mappedFP);

				if (this->contigBubbleInfo[id2Index(contigResult[0][0].id)].oppositeContigID[0] != this->contigBubbleInfo[id2Index(contigResult[0][0].id)].oppositeContigID[1] ||
					this->contigBubbleInfo[id2Index(contigResult[0][1].id)].oppositeContigID[0] != this->contigBubbleInfo[id2Index(contigResult[0][1].id)].oppositeContigID[1]) {
					continue;
				}

				if (abs(contigResult[0][0].id) == abs(this->contigBubbleInfo[id2Index(contigResult[0][1].id)].oppositeContigID[0]))
					continue;

				for (i = 0; i < 2; ++i) {
					contigIndex[0][i] = abs(contigResult[0][i].id) - 1;
					if (contigPositionInScaffold[contigIndex[0][i]].size() != 1 || contigBubbleInfo[contigIndex[0][i]].oppositeContigID[0] == 0) //edited by ouchi (contigPositionInScaffold)
						break;

					contigResult[1][i].id = sign(contigResult[0][i].id) * contigBubbleInfo[contigIndex[0][i]].oppositeContigID[0];
					contigResult[1][i].offset = contigResult[0][i].offset;
					contigIndex[1][i] = abs(contigResult[1][i].id) - 1;
					if (contigPositionInScaffold[contigIndex[1][i]].size() != 1) //edited by ouchi (contigPositionInScaffold)
						break;
				}
				if (i != 2)
					continue;

				for (i = 0; i < 2; ++i) {
					for (j = 0; j < 2; ++j) {
						scaffoldResult[j][i].id = contigResult[j][i].id > 0 ? contigPositionInScaffold[contigIndex[j][i]][0].id : -(contigPositionInScaffold[contigIndex[j][i]][0].id); //edited by ouchi (contigPositionInScaffold)
						scaffoldResult[j][i].offset = contigPositionInScaffold[contigIndex[j][i]][0].id > 0 ? contigResult[j][i].offset : contig[contigIndex[j][i]].length - contigResult[j][i].offset - 1; //edited by ouchi (contigPositionInScaffold)
						scaffoldResult[j][i].offset += node[abs(scaffoldResult[j][i].id) - 1].contig[contigPositionInScaffold[contigIndex[j][i]][0].offset].start; //edited by ouchi (contigPositionInScaffold)

						scaffoldIndex[j][i] = abs(scaffoldResult[j][i].id) - 1;
					}
				}
				if (scaffoldIndex[0][0] != scaffoldIndex[1][1] || scaffoldIndex[1][0] != scaffoldIndex[0][1])
					continue;

				for (j = 0; j < 2; ++j) {
					if (scaffoldIndex[j][0] == scaffoldIndex[j][1] || (int)scaffoldIndex[j][0] != abs(node[id2Index(scaffoldResult[j][1].id)].oppositeBubbleNodeID) - 1)
						break;
				}
				if (j != 2)
					continue;


				for (j = 0; j < 2; ++j) {
					long insertLength = abs(scaffoldResult[j][0].offset - scaffoldResult[j^1][1].offset);
					if (insertLength - averageInsSize > tolerence)
						continue;

					unsigned long leftMostContigNum = std::min(contigPositionInScaffold[contigIndex[j][0]][0].offset, contigPositionInScaffold[contigIndex[j^1][1]][0].offset); //edited by ouchi (contigPositionInScaffold)
					unsigned long rightMostContigNum = std::max(contigPositionInScaffold[contigIndex[j][0]][0].offset, contigPositionInScaffold[contigIndex[j^1][1]][0].offset); //edited by ouchi (contigPositionInScaffold)

					for (unsigned long k = leftMostContigNum; k < rightMostContigNum; ++k) {
						std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldIndex[j][0], k);
						if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
							redundancyCheckSet.insert(redundancyCheckKey);
							++(numPairThread[threadID][scaffoldIndex[j][0]][k]);
						}
					}
				}
			}

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&(contigResult[0][0]), sizeof(platanus::Position), 1, mappedFP);
				fread(&(contigResult[0][1]), sizeof(platanus::Position), 1, mappedFP);
			}
		}
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::countLongReadSwitchErrorLinks(vector<vector<unsigned> >& numPair, const long numThread)
{
	setOppositeBubbleNodeIDAndStateForEachNode();

    omp_set_num_threads(numThread);

    vector<vector<vector<unsigned> > > numPairThread(numThread);
	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		numPairThread[threadID].resize(numPair.size());
		for (unsigned i = 0; i < numPair.size(); ++i)
			numPairThread[threadID][i].resize(numPair[i].size(), 0);
	}

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		std::array<std::array<platanus::Position, 2>, 2> contigResult;
		std::array<std::array<platanus::Position, 2>, 2> scaffoldResult;
		std::array<std::array<unsigned long, 2>, 2> contigIndex;
		std::array<std::array<unsigned long, 2>, 2> scaffoldIndex;
		unsigned i, j;

		long score;
		unsigned numResult;
		FILE *mappedFP;
		vector<platanus::Position> positionBuffer;
		mappedFP = (*longReadLibraryMT)[threadID].mappedReadFP;

        rewind(mappedFP);
		while (fread(&numResult, sizeof(unsigned), 1, mappedFP)) {
			positionBuffer.resize(numResult);
			for (i = 0; i < numResult; ++i) {
				fread(&(positionBuffer[i]), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fseek((*longReadLibraryMT)[threadID].mappedReadFP, 2 * sizeof(int), SEEK_CUR);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
			}

			std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> redundancyCheckSet;

			for (auto positionItr1 = positionBuffer.begin(); positionItr1 != positionBuffer.end() - 1; ++positionItr1) {
				for (auto positionItr2 = positionItr1 + 1; positionItr2 != positionBuffer.end(); ++positionItr2) {
					contigResult[0][0] = *positionItr1;
					contigResult[0][1] = *positionItr2;

					if (this->contigBubbleInfo[id2Index(contigResult[0][0].id)].oppositeContigID[0] != this->contigBubbleInfo[id2Index(contigResult[0][0].id)].oppositeContigID[1] ||
						this->contigBubbleInfo[id2Index(contigResult[0][1].id)].oppositeContigID[0] != this->contigBubbleInfo[id2Index(contigResult[0][1].id)].oppositeContigID[1]) {
						continue;
					}

					if (abs(contigResult[0][0].id) == abs(this->contigBubbleInfo[id2Index(contigResult[0][1].id)].oppositeContigID[0]))
						continue;

					for (i = 0; i < 2; ++i) {
						contigIndex[0][i] = abs(contigResult[0][i].id) - 1;
						if (contigPositionInScaffold[contigIndex[0][i]].size() != 1 || contigBubbleInfo[contigIndex[0][i]].oppositeContigID[0] == 0) //edited by ouchi (contigPositionInScaffold)
							break;

						contigResult[1][i].id = sign(contigResult[0][i].id) * contigBubbleInfo[contigIndex[0][i]].oppositeContigID[0];
						contigResult[1][i].offset = contigResult[0][i].offset;
						contigIndex[1][i] = abs(contigResult[1][i].id) - 1;
						if (contigPositionInScaffold[contigIndex[1][i]].size() != 1) //edited by ouchi (contigPositionInScaffold)
							break;
					}
					if (i != 2)
						continue;

					for (i = 0; i < 2; ++i) {
						for (j = 0; j < 2; ++j) {
							scaffoldResult[j][i].id = contigResult[j][i].id > 0 ? contigPositionInScaffold[contigIndex[j][i]][0].id : -(contigPositionInScaffold[contigIndex[j][i]][0].id); //edited by ouchi (contigPositionInScaffold)
							scaffoldResult[j][i].offset = contigPositionInScaffold[contigIndex[j][i]][0].id > 0 ? contigResult[j][i].offset : contig[contigIndex[j][i]].length - contigResult[j][i].offset - 1; //edited by ouchi (contigPositionInScaffold)
							scaffoldResult[j][i].offset += node[abs(scaffoldResult[j][i].id) - 1].contig[contigPositionInScaffold[contigIndex[j][i]][0].offset].start; //edited by ouchi (contigPositionInScaffold)

							scaffoldIndex[j][i] = abs(scaffoldResult[j][i].id) - 1;
						}
					}
					if (scaffoldIndex[0][0] != scaffoldIndex[1][1] || scaffoldIndex[1][0] != scaffoldIndex[0][1])
						continue;

					for (j = 0; j < 2; ++j) {
						if (scaffoldIndex[j][0] == scaffoldIndex[j][1] || (int)scaffoldIndex[j][0] != abs(node[id2Index(scaffoldResult[j][1].id)].oppositeBubbleNodeID) - 1)
							break;
					}
					if (j != 2)
						continue;


					for (j = 0; j < 2; ++j) {
						unsigned long leftMostContigNum = std::min(contigPositionInScaffold[contigIndex[j][0]][0].offset, contigPositionInScaffold[contigIndex[j^1][1]][0].offset); //edited by ouchi (contigPositionInScaffold)
						unsigned long rightMostContigNum = std::max(contigPositionInScaffold[contigIndex[j][0]][0].offset, contigPositionInScaffold[contigIndex[j^1][1]][0].offset); //edited by ouchi (contigPositionInScaffold)

						for (unsigned long k = leftMostContigNum; k < rightMostContigNum; ++k) {
							std::pair<int, int> redundancyCheckKey = std::make_pair(scaffoldIndex[j][0], k);
							if (redundancyCheckSet.find(redundancyCheckKey) == redundancyCheckSet.end()) {
								redundancyCheckSet.insert(redundancyCheckKey);
								++(numPairThread[threadID][scaffoldIndex[j][0]][k]);
							}
						}
					}
				}
			}
        }
    }

    for (long threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned i = 0; i < numNode; ++i) {
			for (unsigned j = 0; j < numPair[i].size(); ++j) {
				numPair[i][j] += numPairThread[threadID][i][j];
			}
		}
	}
}

void PairedDBG::divideErroneousNode(const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize)
{
	long currentLibraryIndex = this->targetLibraryIndex;

	cerr << "dividing erroneous scaffolds..." << endl;

    vector<vector<unsigned> > numSpanningPair(numNode);
    vector<vector<unsigned> > numErroneousPair(numNode);
    for (unsigned i = 0; i < numNode; ++i) {
        numSpanningPair[i].resize(node[i].numContig - 1, 0);
        numErroneousPair[i].resize(node[i].numContig - 1, 0);
	}

    std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> errorLink;
    vector<vector<double> > sumExpectedLink(numNode);
    omp_set_num_threads(numThread);

	if (allLibraryMT != NULL) {
		for (unsigned libraryID = 0; libraryID < (*allLibraryMT).size(); ++libraryID) {
			this->setTargetLibraryIndex(libraryID);
			this->setTolerence(minTolerenceFactor * (*allLibraryMT)[libraryID][0].getSDInsSize());

			this->countPairsSpanningGap(numSpanningPair, divisionMode, numThread);
			if (divisionMode == SWITCH)
				this->countSwitchErrorLinks(numErroneousPair, numThread);
			else
				this->countLinksInsideContigs(numErroneousPair, numThread);

			if (divisionMode != SWITCH) {
				for (unsigned i = 0; i < numNode; ++i) {
					sumExpectedLink[i].resize(node[i].numContig - 1, 0.0);
					for (unsigned j = 0; j < node[i].numContig - 1; ++j) {
						sumExpectedLink[i][j] += this->calcExpectedLink((*allLibraryMT)[libraryID][0].getAverageCoverage(), node[i].contig[j].end, node[i].length - node[i].contig[j + 1].start, node[i].contig[j + 1].start - node[i].contig[j].end);
					}
				}
			}
		}
	}

	if (longReadLibraryMT != NULL) {
		this->setTolerence((*longReadLibraryMT)[0].getAverageInsSize());

		this->countLongReadSpanningGap(numSpanningPair, divisionMode, numThread);
		if (divisionMode == SWITCH)
			this->countLongReadSwitchErrorLinks(numErroneousPair, numThread);
		else
			this->countLongReadLinksInsideContigs(numErroneousPair, numThread);

		if (divisionMode != SWITCH) {
			for (unsigned i = 0; i < numNode; ++i) {
				sumExpectedLink[i].resize(node[i].numContig - 1, 0.0);
				for (unsigned j = 0; j < node[i].numContig - 1; ++j) {
					sumExpectedLink[i][j] += this->calcExpectedLink((*longReadLibraryMT)[0].getAverageCoverage(), node[i].contig[j].end, node[i].length - node[i].contig[j + 1].start, node[i].contig[j + 1].start - node[i].contig[j].end);
				}
			}
		}
	}


    cerr << "dividing low coverage links..." << endl;
	this->divideErroneousLink(numErroneousPair, numSpanningPair, sumExpectedLink, errorLink, minLink, divisionMode, numThread, maxGapSize);

	this->setTargetLibraryIndex(currentLibraryIndex);
}

void PairedDBG::makeScaffold(void)
{
    long numNewContig = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;
    vector<GraphLayout> include, candidate;
    include.reserve(10000);
    candidate.reserve(10000);
    vector<long> numOccurrence(numNode, 0); //added by ouchi

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
            if ((node[tmpNodeID].state & SC_INC) && (!(node[tmpNodeID].state & SC_REP) || (node[tmpNodeID].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))))
				continue;

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


            long tmpNodeID = abs(candidate[minCandidateID].id) - 1;

            if ((node[tmpNodeID].state & SC_INC) && (!(node[tmpNodeID].state & SC_REP) || (node[tmpNodeID].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))) {
                auto candidateIterator = candidate.begin();
                candidateIterator = candidate.erase(candidateIterator + minCandidateID);
                continue;
            }
            unsigned j = 0;
            for (; j < include.size(); ++j) {
				long overlapTolerence = std::min(tolerence, std::min(candidate[minCandidateID].end - candidate[minCandidateID].start, include[j].end - include[j].start) / 2);

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
						if ((node[tmpNodeID].state & SC_INC) && (!(node[tmpNodeID].state & SC_REP) || (node[tmpNodeID].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE))))
							continue;

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
        long includeSize = include.size();
        long j = 0;
        for (; node[abs(include[j].id) - 1].state & SC_REP; ++j)
            numNewContig -= node[abs(include[j].id) - 1].numContig;
        for (; node[abs(include[includeSize - 1].id) - 1].state & SC_REP; --includeSize)
            numNewContig -= node[abs(include[includeSize - 1].id) - 1].numContig;
        fwrite(&numNewContig, sizeof(long), 1, scaffoldFP);

        long minStart = include[j].start;
/*
		for (long k = j + 1; k < includeSize; ++k) {
			if (minStart > include[k].start)
				minStart = include[k].start;
		}
*/

        for (; j < includeSize; ++j) {
            long tmpNodeID = abs(include[j].id) - 1;

            node[tmpNodeID].state |= SC_INC;
            ++numOccurrence[tmpNodeID]; //added by ouchi

            include[j].start -= minStart;
            include[j].end -= minStart;

            if (include[j].start != 0) {
                long overlapLength = this->getScaffoldOverlap(include[j-1].id, include[j].id);
                if (overlapLength > 0 && overlapLength + include[j].start - include[j-1].end <= tolerence) {
                    overlapLength = include[j-1].end - include[j].start - overlapLength;
                    for (long k = j; k < includeSize; ++k) {
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
            //edited by ouchi
            /*vector<long> leftEdge;
            vector<long> rightEdge;
            std::sort(node[i].edge.begin(), node[i].edge.end(), [](GraphEdge const& left, GraphEdge const& right){
                return left.length < right.length;
            });
            for (long edgeID = 0; edgeID < node[i].numEdge; ++edgeID) {
                if (node[i].edge[edgeID].direction < 0) {
                    long j = 0;
                    for (;j < leftEdge.size(); ++j) {
                        if (!(this->checkNodeConflict(node[i].edge[edgeID], node[i].edge[leftEdge[j]], node[id2Index(node[i].edge[edgeID].end)], node[id2Index(node[i].edge[leftEdge[j]].end)])))
                            break;
                    }
                    if (j < leftEdge.size())
                        continue;
                    leftEdge.push_back(edgeID);
                    if (leftEdge.size() > numOccurrence[i])
                        break;
                } else {
                    long j = 0;
                    for (;j < rightEdge.size(); ++j) {
                        if (!(this->checkNodeConflict(node[i].edge[edgeID], node[i].edge[rightEdge[j]], node[id2Index(node[i].edge[edgeID].end)], node[id2Index(node[i].edge[rightEdge[j]].end)])))
                            break;
                    }
                    if (j < rightEdge.size())
                        continue;
                    rightEdge.push_back(edgeID);
                    if (rightEdge.size() > numOccurrence[i])
                        break;
                }
            }
            if (leftEdge.size() > numOccurrence[i] || rightEdge.size() > numOccurrence[i]) {
                fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP);
                for (long j = 0; j < numContig; ++j)
                    fwrite(&(node[i].contig[j]), sizeof(ScaffoldPart), 1, scaffoldFP);
                newContigPoolSize += node[i].numContig;
                ++numNewNode;
            }
            // */

            for (long j = 0; j < numContig; ++j) //deleted by ouchi
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

void PairedDBG::clearContigPreviousParentNodeID()
{
	this->contigPreviousParentNodeID.clear();
}

void PairedDBG::setOppositeBubbleContigIDByEndMatch()
{
	const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigBubbleInfo.empty())
		this->contigBubbleInfo.resize(this->numContig);

	std::unordered_map<string, std::array<long, 3> > endSeqMap;

	for (unsigned long i = this->node.size() - this->numInputBubbleContig; i < this->node.size(); ++i) {
		if (this->node[i].length >= this->contigMaxK) {
			string endSeq = contig[i].base.substr(0, contigMaxK - 1) + contig[i].base.substr(this->contig[i].base.size() - (contigMaxK - 1), contigMaxK - 1);
			auto mapItr = endSeqMap.find(endSeq);

			if (mapItr == endSeqMap.end()) {
				std::array<long, 3> endSeqInfo;
				endSeqInfo[0] = i + 1;
				endSeqInfo[1] = 0;
				endSeqInfo[2] = 1;
				endSeqMap[endSeq] = endSeqInfo;
			}
			else {
				if ((mapItr->second)[2] == 1)
					(mapItr->second)[1] = i + 1;
				++(mapItr->second)[2];
			}
		}
	}

	for (auto mapItr = endSeqMap.begin(); mapItr != endSeqMap.end(); ++mapItr) {
		if ((mapItr->second)[2] != 2)
			continue;

		if (calcNodeCoverage(node[(mapItr->second)[0] - 1]) > COVERAGE_THRESHOLD || calcNodeCoverage(node[(mapItr->second)[1] - 1]) > COVERAGE_THRESHOLD)
			continue;

		this->contigBubbleInfo[(mapItr->second)[0] - 1].oppositeContigID.fill((mapItr->second)[1]);
		this->contigBubbleInfo[(mapItr->second)[1] - 1].oppositeContigID.fill((mapItr->second)[0]);
	}
}

void PairedDBG::setOppositeBubbleContigIDByOneEndMatch()
{
	const double COVERAGE_THRESHOLD = HETERO_FORK_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

	if (this->contigBubbleInfo.empty())
		this->contigBubbleInfo.resize(this->numContig);

	std::unordered_map<string, std::array<long, 3> > endSeqMap;
	std::array<string, 2> endSeq;
	endSeq[0].resize(contigMaxK - 1);
	endSeq[1].resize(contigMaxK - 1);

	std::array<bool, 2> gapFlag;

	for (unsigned long i = this->node.size() - this->numInputBubbleContig; i < this->node.size(); ++i) {
		if (this->node[i].length < this->contigMaxK)
			continue;

		gapFlag.fill(false);

		for (long j = 0; j < contigMaxK - 1; ++j) {
			if (contig[i].base[j] < 4)
				endSeq[0][j] = contig[i].base[j];
			else
				gapFlag[0] = true;

			if (contig[i].base[contig[i].base.size() - (contigMaxK - 1) + j] < 4)
				endSeq[1][j] = 3^(contig[i].base[contig[i].base.size() - 1 - j]);
			else
				gapFlag[1] = true;
		}

		for (long j = 0; j < 2; ++j) {
			if (gapFlag[j])
				continue;

			auto mapItr = endSeqMap.find(endSeq[j]);

			if (mapItr == endSeqMap.end()) {
				std::array<long, 3> endSeqInfo;
				endSeqInfo[0] = (j == 0 ? i + 1 : -(i + 1));
				endSeqInfo[1] = 0;
				endSeqInfo[2] = 1;
				endSeqMap[endSeq[j]] = endSeqInfo;
			}
			else {
				if ((mapItr->second)[2] == 1)
					(mapItr->second)[1] = (j == 0 ? i + 1 : -(i + 1));
				++(mapItr->second)[2];
			}
		}
	}

	for (auto mapItr = endSeqMap.begin(); mapItr != endSeqMap.end(); ++mapItr) {
		if ((mapItr->second)[2] != 2)
			continue;

		if (calcNodeCoverage(node[id2Index((mapItr->second)[0])]) > COVERAGE_THRESHOLD || calcNodeCoverage(node[id2Index((mapItr->second)[1])]) > COVERAGE_THRESHOLD)
			continue;

		for (long j = 0; j < 2; ++j) {
			if ((mapItr->second)[j] > 0 && contigBubbleInfo[(mapItr->second)[j] - 1].oppositeContigID[0] == 0)
				contigBubbleInfo[(mapItr->second)[j] - 1].oppositeContigID[0] = (mapItr->second)[(j + 1)%2];
			else if ((mapItr->second)[j] < 0 && contigBubbleInfo[-(mapItr->second)[j] - 1].oppositeContigID[1] == 0)
				contigBubbleInfo[-(mapItr->second)[j] - 1].oppositeContigID[1] = -(mapItr->second)[(j + 1)%2];
		}
	}

/*
	for (unsigned long i = 0; i < this->node.size(); ++i) {
		std::cout << contigName[i] <<  '\t' << contigBubbleInfo[i].oppositeContigID[0] <<  '\t' << contigBubbleInfo[i].oppositeContigID[1] << endl;
	}
	exit(0);
*/
}

long PairedDBG::getScoreFromIDPair(long leftNodeID, long rightNodeID)
{
	GraphNode &leftNode = this->node[abs(leftNodeID) - 1];
	char directionFromLeft = sign(leftNodeID);
	char rightStrand = directionFromLeft * sign(rightNodeID);

	for (long i = 0; i < leftNode.numEdge; ++i) {
		if (abs(leftNode.edge[i].end) == abs(rightNodeID) &&
			leftNode.edge[i].direction == directionFromLeft &&
			sign(leftNode.edge[i].end) == rightStrand) {

			return leftNode.edge[i].score;
		}
	}

	return 0;
}

void PairedDBG::insertSizeDistribution(vector<SeqLib>& library, vector<long>& distribution, const long numThread)
{
    long insertLength;
    platanus::Position forwardResult;
    platanus::Position reverseResult;

    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numPairLink;
		unsigned numReadLink;

        rewind(library[threadID].mappedFP);
        while (fread(&numPairLink, sizeof(unsigned), 1, library[threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, library[threadID].mappedFP);

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);

				if (pairIndex > 0)
					continue;

				unsigned long forwardIndex = abs(forwardResult.id) - 1;
				if (contigPositionInScaffold[forwardIndex].size() != 1) continue; //edited by ouchi (contigPositionInScaffold)
				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardIndex][0].id : -(contigPositionInScaffold[forwardIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset = contigPositionInScaffold[forwardIndex][0].id > 0 ? forwardResult.offset : contig[forwardIndex].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset += node[abs(forwardResult.id) - 1].contig[contigPositionInScaffold[forwardIndex][0].offset].start;

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

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, library[threadID].mappedFP);
			}
        }
    }
}

void PairedDBG::updateInsertLengthFP(vector<SeqLib>& lib, const long numThread)
{
    platanus::Position forwardResult;
    platanus::Position reverseResult;

	fseek(lib[0].insertLengthFP, 0, SEEK_END);

    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numPairLink;
		unsigned numReadLink;

        rewind(lib[threadID].mappedFP);
        while (fread(&numPairLink, sizeof(unsigned), 1, lib[threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, lib[threadID].mappedFP);

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);

				if (pairIndex > 0)
					continue;

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

			for (unsigned pairIndex = 0; pairIndex < numReadLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, lib[threadID].mappedFP);
			}
        }
    }
}

void PairedDBG::markRedundantResultSeq(const long numThread)
{
	unsigned long prefixLength = this->contigMaxK;

	for (unsigned long seqIndex = 0; seqIndex < resultSeq.size(); ++seqIndex) {
		if (resultSeq[seqIndex].seq.empty() || (node[seqIndex].state & SC_DEL))
			continue;

		if (prefixLength > resultSeq[seqIndex].seq.size())
			prefixLength = resultSeq[seqIndex].seq.size();
	}

	std::unordered_map<string, vector<long> > prefixToSeqIndex;

	for (unsigned long seqIndex = 0; seqIndex < resultSeq.size(); ++seqIndex) {
		auto mapItr = prefixToSeqIndex.find(resultSeq[seqIndex].seq.substr(0, prefixLength));
		if (mapItr == prefixToSeqIndex.end())
			prefixToSeqIndex[resultSeq[seqIndex].seq.substr(0, prefixLength)] = std::vector<long>(1, seqIndex);
		else
		 	mapItr->second.push_back(seqIndex);
	}

    omp_set_num_threads(numThread);
	vector<vector<char> > redundantFlag(numThread);

	# pragma omp parallel for schedule(static, 1)
	for (unsigned threadID = 0; threadID < numThread; ++threadID) {
		redundantFlag[threadID].assign(resultSeq.size(), false);

		for (unsigned long seqIndex = threadID; seqIndex < resultSeq.size(); seqIndex += numThread) {
			if (resultSeq[seqIndex].seq.empty() || (node[seqIndex].state & SC_DEL))
				continue;

			for (unsigned strandIndex = 0; strandIndex < 2; ++strandIndex) {
				string targetSeq;
				if (strandIndex == 0) {
					targetSeq = resultSeq[seqIndex].seq;
				}
				else {
					targetSeq.resize(resultSeq[seqIndex].seq.size());
					for (unsigned long baseIndex = 0; baseIndex < resultSeq[seqIndex].seq.size(); ++baseIndex) {
						if (resultSeq[seqIndex].seq[baseIndex] < 4)
							targetSeq[targetSeq.size() - baseIndex - 1] = 3^(resultSeq[seqIndex].seq[baseIndex]);
						else
							targetSeq[targetSeq.size() - baseIndex - 1] = 4;
					}
				}

				for (unsigned long baseIndex = 0; baseIndex < targetSeq.size() - prefixLength + 1; ++baseIndex) {
					auto mapItr = prefixToSeqIndex.find(targetSeq.substr(baseIndex, prefixLength));
					if (mapItr == prefixToSeqIndex.end())
						continue;

					for (auto vecItr = mapItr->second.begin(); vecItr != mapItr->second.end(); ++vecItr) {
						if (*vecItr != seqIndex &&
							!(node[*vecItr].state & SC_DEL) &&
							(targetSeq.size() > resultSeq[*vecItr].seq.size() || (targetSeq.size() == resultSeq[*vecItr].seq.size() && seqIndex < *vecItr)) &&
							targetSeq.size() - baseIndex >= resultSeq[*vecItr].seq.size() &&
							std::search(resultSeq[*vecItr].seq.begin()+prefixLength, resultSeq[*vecItr].seq.end(), targetSeq.begin()+baseIndex+prefixLength, targetSeq.begin()+baseIndex+resultSeq[*vecItr].seq.size()) == resultSeq[*vecItr].seq.begin()+prefixLength) {

							redundantFlag[threadID][*vecItr] = true;
						}
					}
				}
			}
		}
	}

	for (unsigned threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned long seqIndex = 0; seqIndex < resultSeq.size(); ++seqIndex)
			if (resultSeq[seqIndex].redundantFlag == false && redundantFlag[threadID][seqIndex])
				resultSeq[seqIndex].redundantFlag = true;
	}
}


void PairedDBG::divideErroneousNodeBaseLevel(const long pairLibraryBegin, const long pairLibraryEnd, const long readLength, const bool bubbleFlag, const bool longLibraryFlag, const bool storeOnlyFlag, const long numThread)
{
	long currentLibraryIndex = this->targetLibraryIndex;

	cerr << "dividing erroneous scaffolds based on base-level coverages ..." << endl;

    vector<vector<unsigned> > physicalCoverage(numNode);
    vector<vector<unsigned> > diffCoverage(numNode);
    for (unsigned i = 0; i < numNode; ++i) {
        physicalCoverage[i].resize(node[i].length, 0);
        diffCoverage[i].resize(node[i].length, 0);
	}

    std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> errorLink;
    omp_set_num_threads(numThread);

	if (allLibraryMT != NULL) {
		for (unsigned libraryID = pairLibraryBegin; libraryID < pairLibraryEnd; ++libraryID) {
			this->setTargetLibraryIndex(libraryID);
			this->setTolerence(minTolerenceFactor * (*allLibraryMT)[libraryID][0].getSDInsSize());

			this->calculatePhysicalCoverage(physicalCoverage, (*allLibraryMT)[libraryID][0].getAverageInsSize(), numThread);
			this->compensatePhysicalCoverageBasedOnGapRate(physicalCoverage, 2 * (*allLibraryMT)[libraryID][0].getAverageInsSize(), numThread);

			if (storeOnlyFlag)
				this->calculateDiffCoverage(diffCoverage, 0, numThread);
			else
				this->calculateDiffCoverage(diffCoverage, (*allLibraryMT)[libraryID][0].getAverageInsSize(), numThread);
		}

	}

	if (longLibraryFlag && longReadLibraryMT != NULL) {
		this->setTolerence((*longReadLibraryMT)[0].getAverageInsSize());

		this->calculateLongReadPhysicalCoverage(physicalCoverage, (*longReadLibraryMT)[0].getAverageInsSize(), numThread);
		this->compensatePhysicalCoverageBasedOnGapRate(physicalCoverage, 2 * (*longReadLibraryMT)[0].getAverageInsSize(), numThread);

		if (storeOnlyFlag)
			this->calculateLongReadDiffCoverage(diffCoverage, 0, numThread);
		else
			this->calculateLongReadDiffCoverage(diffCoverage, (*longReadLibraryMT)[0].getAverageInsSize(), numThread);
	}

/*
	string outFile = "out_physical.depth";
	this->dumpLongestNodeCoverage(physicalCoverage, numNode, outFile);

	outFile = "out_diff.depth";
	this->dumpLongestNodeCoverage(diffCoverage, numNode, outFile);
*/

	if (storeOnlyFlag)
		this->storeBreakpointBasedOnCoverage(physicalCoverage, diffCoverage, numThread);
	else
		this->divideNodeBasedOnCoverage(physicalCoverage, diffCoverage, bubbleFlag, numThread);

	this->setTargetLibraryIndex(currentLibraryIndex);
}


void PairedDBG::calculatePhysicalCoverage(vector<vector<unsigned> >& physicalCoverage, const long insertTolerence, const long numThread)
{
    const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
	const long readLength = (*allLibraryMT)[targetLibraryIndex][0].getAverageLength();

    omp_set_num_threads(numThread);
    omp_lock_t lock[node.size()];
    for (unsigned i = 0; i < node.size(); ++i)
        omp_init_lock(&lock[i]);

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		unsigned numPairLink;
		unsigned numReadLink;

		rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				unsigned long forwardContigIndex = id2Index(forwardResult.id);
				if (contigPositionInScaffold[forwardContigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;

				unsigned long reverseContigIndex = id2Index(reverseResult.id);
				if (contigPositionInScaffold[reverseContigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;

				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardContigIndex][0].id : -(contigPositionInScaffold[forwardContigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseContigIndex][0].id : -(contigPositionInScaffold[reverseContigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				if (forwardResult.id != -reverseResult.id)
					continue;

				unsigned long nodeIndex = id2Index(forwardResult.id);

				forwardResult.offset = contigPositionInScaffold[forwardContigIndex][0].id > 0 ? forwardResult.offset : contig[forwardContigIndex].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset += node[nodeIndex].contig[contigPositionInScaffold[forwardContigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				reverseResult.offset = contigPositionInScaffold[reverseContigIndex][0].id > 0 ? reverseResult.offset : contig[reverseContigIndex].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset += node[nodeIndex].contig[contigPositionInScaffold[reverseContigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				long insertLength = 0;
				if (forwardResult.id > 0 && forwardResult.offset < reverseResult.offset)
					insertLength = static_cast<long>(reverseResult.offset - forwardResult.offset + 1);
				else if (reverseResult.id > 0 && reverseResult.offset < forwardResult.offset)
					insertLength = static_cast<long>(forwardResult.offset - reverseResult.offset + 1);
				else
					continue;

				if (abs(insertLength - averageInsSize) > insertTolerence)
					continue;

				std::pair<int, int> range = std::minmax(
					std::min(std::max(forwardResult.offset, 0), static_cast<int>(node[nodeIndex].length - 1)),
					std::min(std::max(reverseResult.offset, 0), static_cast<int>(node[nodeIndex].length - 1))
				);
				range.first = std::min(range.first + readLength, node[nodeIndex].length - 1);
				range.second = std::max(range.second - readLength, 0L);

				omp_set_lock(&lock[nodeIndex]);
				for (long i = range.first; i <= range.second; ++i) {
					++(physicalCoverage[nodeIndex][i]);
				}
				omp_unset_lock(&lock[nodeIndex]);

				fseek((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP, 2 * (numPairLink - pairIndex - 1) * sizeof(platanus::Position), SEEK_CUR);
				break;
			}

			fseek((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP, 2 * numReadLink * sizeof(platanus::Position), SEEK_CUR);
        }


		platanus::Position leftPosition;
		long insertLength;
        rewind((*allLibraryMT)[targetLibraryIndex][threadID].mappedSameContigFP);
        while (fread(&leftPosition, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedSameContigFP)) {
			fread(&insertLength, sizeof(long), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedSameContigFP);

			if (abs(insertLength - averageInsSize) > insertTolerence)
				continue;

			unsigned long contigIndex = id2Index(leftPosition.id);
			if (contigPositionInScaffold[contigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
				continue;

			long nodeIndex = id2Index(contigPositionInScaffold[contigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
			long nodeLeftOffset = contigPositionInScaffold[contigIndex][0].id > 0 ? leftPosition.offset : contig[contigIndex].length - forwardResult.offset - insertLength; //edited by ouchi (contigPositionInScaffold)
			nodeLeftOffset += node[nodeIndex].contig[contigPositionInScaffold[contigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

			std::pair<long, long> range = std::make_pair(
				std::min(std::max(nodeLeftOffset, 0L), node[nodeIndex].length - 1),
				std::min(std::max(nodeLeftOffset + insertLength, 0L), node[nodeIndex].length - 1)
			);
			range.first = std::min(range.first + readLength, node[nodeIndex].length - 1);
			range.second = std::max(range.second - readLength, 0L);

			omp_set_lock(&lock[nodeIndex]);
			for (long i = range.first; i <= range.second; ++i) {
				++(physicalCoverage[nodeIndex][i]);
			}
			omp_unset_lock(&lock[nodeIndex]);
       }
    }
}


void PairedDBG::calculateLongReadPhysicalCoverage(vector<vector<unsigned> >& physicalCoverage, const long insertTolerence, const long numThread)
{
    omp_set_num_threads(numThread);
    omp_lock_t lock[node.size()];
    for (unsigned i = 0; i < node.size(); ++i)
        omp_init_lock(&lock[i]);

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numResult;
		int targetStart;
		int targetEnd;

        rewind((*longReadLibraryMT)[threadID].mappedReadFP);
		while (fread(&numResult, sizeof(unsigned), 1, (*longReadLibraryMT)[threadID].mappedReadFP)) {
			std::unordered_map<int, std::pair<int, int> > index2Range;
			platanus::Position position;

			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&(position), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&targetStart, sizeof(int), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&targetEnd, sizeof(int), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fseek((*longReadLibraryMT)[threadID].mappedReadFP, sizeof(long), SEEK_CUR);

				unsigned long contigIndex = id2Index(position.id);
				if (contigPositionInScaffold[contigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;

				int nodeID = position.id > 0 ? contigPositionInScaffold[contigIndex][0].id : -(contigPositionInScaffold[contigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				int nodeIndex = id2Index(nodeID);

				targetStart = contigPositionInScaffold[contigIndex][0].id > 0 ? targetStart : contig[contigIndex].length - targetStart - 1; //edited by ouchi (contigPositionInScaffold)
				targetStart += node[nodeIndex].contig[contigPositionInScaffold[contigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				targetEnd = contigPositionInScaffold[contigIndex][0].id > 0 ? targetEnd : contig[contigIndex].length - targetEnd - 1; //edited by ouchi (contigPositionInScaffold)
				targetEnd += node[nodeIndex].contig[contigPositionInScaffold[contigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				if (targetStart >targetEnd)
					std::swap(targetStart, targetEnd);

				targetStart = std::min(static_cast<long>(std::max(targetStart, 0)), node[nodeIndex].length - 1);
				targetEnd = std::min(static_cast<long>(std::max(targetEnd, 0)), node[nodeIndex].length - 1);

				auto mapItr = index2Range.find(nodeIndex);
				if (mapItr == index2Range.end()) {
					index2Range.insert(std::make_pair(nodeIndex, std::make_pair(targetStart, targetEnd)));
				}
				else {
					mapItr->second.first = std::min(mapItr->second.first, targetStart);
					mapItr->second.second = std::max(mapItr->second.second, targetEnd);
				}
			}

			for (auto mapItr = index2Range.begin(); mapItr != index2Range.end(); ++mapItr) {
				omp_set_lock(&lock[mapItr->first]);
				for (long i = mapItr->second.first; i <= mapItr->second.second; ++i) {
					++(physicalCoverage[mapItr->first][i]);
				}
				omp_unset_lock(&lock[mapItr->first]);
			}

        }
	}
}


void PairedDBG::compensatePhysicalCoverageBasedOnGapRate(vector<vector<unsigned> >& physicalCoverage, const long windowSize, const long numThread)
{
    omp_set_num_threads(numThread);

	# pragma omp parallel for schedule(dynamic)
	for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
		vector<char> gapBaseFlag;
		node2GapFlagsUnmappableContig(node[nodeIndex], gapBaseFlag);

		vector<char> gapFlag(gapBaseFlag.size() + windowSize, 0);
		std::fill(gapFlag.begin(), gapFlag.begin() + windowSize/2, 1);
		std::fill(gapFlag.begin() + gapBaseFlag.size() + windowSize/2, gapFlag.end(), 1);
		for (long baseIndex = 0; baseIndex < gapBaseFlag.size(); ++baseIndex) {
			if (gapBaseFlag[baseIndex] == 1)
				gapFlag[baseIndex + windowSize/2] = 1;
		}

		long numGapSite = std::count(gapFlag.begin(), gapFlag.begin() + windowSize, 1);
		long baseIndex = 0;
		if (numGapSite > 0)
			physicalCoverage[nodeIndex][baseIndex] = physicalCoverage[nodeIndex][baseIndex] * (windowSize / numGapSite);
		else
			physicalCoverage[nodeIndex][baseIndex] = 0;

		for (baseIndex = 1; baseIndex < gapBaseFlag.size(); ++baseIndex) {
			numGapSite -= gapFlag[baseIndex - 1];
			numGapSite += gapFlag[baseIndex + windowSize - 1];

			physicalCoverage[nodeIndex][baseIndex] = physicalCoverage[nodeIndex][baseIndex] * (windowSize + 1) / (windowSize - numGapSite + 1);
		}
	}
}


void PairedDBG::calculateDiffCoverage(vector<vector<unsigned> >& diffCoverage, const long lengthThreshold, const long numThread)
{
	const long averageInsSize = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
	const int readLength = (*allLibraryMT)[targetLibraryIndex][0].getAverageLength();

    omp_set_num_threads(numThread);
    omp_lock_t lock[node.size()];
    for (unsigned i = 0; i < node.size(); ++i)
        omp_init_lock(&lock[i]);

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		platanus::Position forwardResult;
		platanus::Position reverseResult;
		unsigned numPairLink;
		unsigned numReadLink;
		std::unordered_set<long> forwardRedundancyCheckSet;
		std::unordered_set<long> reverseRedundancyCheckSet;

		FILE *mappedFP = (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP;
        rewind(mappedFP);
		while (fread(&numPairLink, sizeof(unsigned), 1, mappedFP)) {
			fread(&numReadLink, sizeof(unsigned), 1, mappedFP);

			forwardRedundancyCheckSet.clear();
			reverseRedundancyCheckSet.clear();

			for (unsigned pairIndex = 0; pairIndex < numPairLink; ++pairIndex) {
				fread(&forwardResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);
				fread(&reverseResult, sizeof(platanus::Position), 1, (*allLibraryMT)[targetLibraryIndex][threadID].mappedFP);

				unsigned long forwardContigIndex = id2Index(forwardResult.id);
				if (contigPositionInScaffold[forwardContigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;

				unsigned long reverseContigIndex = id2Index(reverseResult.id);
				if (contigPositionInScaffold[reverseContigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;

				forwardResult.id = forwardResult.id > 0 ? contigPositionInScaffold[forwardContigIndex][0].id : -(contigPositionInScaffold[forwardContigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				reverseResult.id = reverseResult.id > 0 ? contigPositionInScaffold[reverseContigIndex][0].id : -(contigPositionInScaffold[reverseContigIndex][0].id); //edited by ouchi (contigPositionInScaffold)

				if (forwardResult.id == -reverseResult.id)
					continue;

				unsigned long forwardNodeIndex = id2Index(forwardResult.id);
				forwardResult.offset = contigPositionInScaffold[forwardContigIndex][0].id > 0 ? forwardResult.offset : contig[forwardContigIndex].length - forwardResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				forwardResult.offset += node[forwardNodeIndex].contig[contigPositionInScaffold[forwardContigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				unsigned long reverseNodeIndex = id2Index(reverseResult.id);
				reverseResult.offset = contigPositionInScaffold[reverseContigIndex][0].id > 0 ? reverseResult.offset : contig[reverseContigIndex].length - reverseResult.offset - 1; //edited by ouchi (contigPositionInScaffold)
				reverseResult.offset += node[reverseNodeIndex].contig[contigPositionInScaffold[reverseContigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)

				if (node[reverseNodeIndex].length >= lengthThreshold && forwardRedundancyCheckSet.find(forwardNodeIndex) == forwardRedundancyCheckSet.end()) {
					forwardRedundancyCheckSet.insert(forwardNodeIndex);
					std::pair<int, int> range;
					if (forwardResult.id > 0) {
						range.first = std::min(std::max(forwardResult.offset + readLength, 0), static_cast<int>(node[forwardNodeIndex].length - 1));
						range.second = std::min(std::max(forwardResult.offset + readLength + averageInsSize, 0L), node[forwardNodeIndex].length - 1);
					}
					else {
						range.first = std::min(std::max(forwardResult.offset - readLength - averageInsSize + 1, 0L), node[forwardNodeIndex].length - 1);
						range.second = std::min(std::max(forwardResult.offset - readLength, 0), static_cast<int>(node[forwardNodeIndex].length - 1));
					}

					omp_set_lock(&lock[forwardNodeIndex]);
					for (long i = range.first; i <= range.second; ++i) {
						++(diffCoverage[forwardNodeIndex][i]);
					}
					omp_unset_lock(&lock[forwardNodeIndex]);
				}

				if (node[forwardNodeIndex].length >= lengthThreshold && reverseRedundancyCheckSet.find(reverseNodeIndex) == reverseRedundancyCheckSet.end()) {
					reverseRedundancyCheckSet.insert(reverseNodeIndex);
					std::pair<int, int> range;
					if (reverseResult.id > 0) {
						range.first = std::min(std::max(reverseResult.offset + readLength, 0), static_cast<int>(node[reverseNodeIndex].length - 1));
						range.second = std::min(std::max(reverseResult.offset + readLength + averageInsSize, 0L), node[reverseNodeIndex].length - 1);
					}
					else {
						range.first = std::min(std::max(reverseResult.offset - readLength - averageInsSize + 1, 0L), node[reverseNodeIndex].length - 1);
						range.second = std::min(std::max(reverseResult.offset - readLength, 0), static_cast<int>(node[reverseNodeIndex].length - 1));
					}

					omp_set_lock(&lock[reverseNodeIndex]);
					for (long i = range.first; i <= range.second; ++i) {
						++(diffCoverage[reverseNodeIndex][i]);
					}
					omp_unset_lock(&lock[reverseNodeIndex]);
				}
			}

			fseek((*allLibraryMT)[targetLibraryIndex][threadID].mappedFP, 2 * numReadLink * sizeof(platanus::Position), SEEK_CUR);
		}
	}
}


void PairedDBG::calculateLongReadDiffCoverage(vector<vector<unsigned> >& diffCoverage, const long lengthThreshold, const long numThread)
{
	const long averageInsSize = (*longReadLibraryMT)[0].getAverageInsSize() / 2;

    omp_set_num_threads(numThread);
    omp_lock_t lock[node.size()];
    for (unsigned i = 0; i < node.size(); ++i)
        omp_init_lock(&lock[i]);

	# pragma omp parallel for schedule(static, 1)
    for (long threadID = 0; threadID < numThread; ++threadID) {
		unsigned numResult;
		int targetStart;
		long score;

        rewind((*longReadLibraryMT)[threadID].mappedReadFP);
		while (fread(&numResult, sizeof(unsigned), 1, (*longReadLibraryMT)[threadID].mappedReadFP)) {
			std::unordered_map<int, std::pair<int, long> > id2OffsetScore;
			platanus::Position position;

			for (unsigned long i = 0; i < numResult; ++i) {
				fread(&(position), sizeof(platanus::Position), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fread(&targetStart, sizeof(int), 1, (*longReadLibraryMT)[threadID].mappedReadFP);
				fseek((*longReadLibraryMT)[threadID].mappedReadFP, sizeof(int), SEEK_CUR);
				fread(&score, sizeof(long), 1, (*longReadLibraryMT)[threadID].mappedReadFP);

				unsigned long contigIndex = id2Index(position.id);
				if (contigPositionInScaffold[contigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
					continue;

				int nodeID = position.id > 0 ? contigPositionInScaffold[contigIndex][0].id : -(contigPositionInScaffold[contigIndex][0].id); //edited by ouchi (contigPositionInScaffold)
				int nodeIndex = id2Index(nodeID);

				targetStart = contigPositionInScaffold[contigIndex][0].id > 0 ? targetStart : contig[contigIndex].length - targetStart - 1; //edited by ouchi (contigPositionInScaffold)
				targetStart += node[nodeIndex].contig[contigPositionInScaffold[contigIndex][0].offset].start; //edited by ouchi (contigPositionInScaffold)
				targetStart = std::min(static_cast<long>(std::max(targetStart, 0)), node[nodeIndex].length - 1);

				auto mapItr = id2OffsetScore.find(nodeID);
				if (mapItr == id2OffsetScore.end()) {
					id2OffsetScore.insert(std::make_pair(nodeID, std::make_pair(targetStart, score)));
				}
				else if (mapItr->second.second < score) {
					mapItr->second = std::make_pair(targetStart, score);
				}
			}

			if (id2OffsetScore.size() <= 1)
				continue;

			vector<platanus::Position> positionBuffer;
			for (auto mapItr = id2OffsetScore.begin(); mapItr != id2OffsetScore.end(); ++mapItr) {
				positionBuffer.emplace_back(mapItr->second.first, mapItr->first);
			}


			for (long i = 0; i < positionBuffer.size() - 1; ++i) {
				for (long j = i + 1; j < positionBuffer.size(); ++j) {
					platanus::Position &forwardResult = positionBuffer[i];
					platanus::Position &reverseResult = positionBuffer[j];

					if (forwardResult.id == -reverseResult.id)
						continue;

					unsigned long forwardNodeIndex = id2Index(forwardResult.id);
					unsigned long reverseNodeIndex = id2Index(reverseResult.id);

					if (node[reverseNodeIndex].length >= lengthThreshold) {
						std::pair<int, int> range;
						if (forwardResult.id > 0) {
							range.first = std::min(std::max(forwardResult.offset, 0), static_cast<int>(node[forwardNodeIndex].length - 1));
							range.second = std::min(std::max(forwardResult.offset + averageInsSize, 0L), node[forwardNodeIndex].length - 1);
						}
						else {
							range.first = std::min(std::max(forwardResult.offset - averageInsSize + 1, 0L), node[forwardNodeIndex].length - 1);
							range.second = std::min(std::max(forwardResult.offset, 0), static_cast<int>(node[forwardNodeIndex].length - 1));
						}

						omp_set_lock(&lock[forwardNodeIndex]);
						for (long i = range.first; i <= range.second; ++i) {
							++(diffCoverage[forwardNodeIndex][i]);
						}
						omp_unset_lock(&lock[forwardNodeIndex]);
					}

					if (node[forwardNodeIndex].length >= lengthThreshold) {
						std::pair<int, int> range;
						if (reverseResult.id > 0) {
							range.first = std::min(std::max(reverseResult.offset, 0), static_cast<int>(node[reverseNodeIndex].length - 1));
							range.second = std::min(std::max(reverseResult.offset + averageInsSize, 0L), node[reverseNodeIndex].length - 1);
						}
						else {
							range.first = std::min(std::max(reverseResult.offset - averageInsSize + 1, 0L), node[reverseNodeIndex].length - 1);
							range.second = std::min(std::max(reverseResult.offset, 0), static_cast<int>(node[reverseNodeIndex].length - 1));
						}

						omp_set_lock(&lock[reverseNodeIndex]);
						for (long i = range.first; i <= range.second; ++i) {
							++(diffCoverage[reverseNodeIndex][i]);
						}
						omp_unset_lock(&lock[reverseNodeIndex]);
					}
				}
			}

		}
	}
}


void PairedDBG::dumpLongestNodeCoverage(vector<vector<unsigned> >& coverage, const long numOutputNode, string &outputFilename)
{
    std::ofstream out(outputFilename.c_str());

	vector<std::pair<int, int> > buffer;
    for (long i = 0 ; i < coverage.size(); ++i)
		buffer.push_back(std::make_pair(node[i].length, i));

	std::stable_sort(buffer.begin(), buffer.end(), platanus::PairFirstGreater());

	for (long i = 0; i < std::min(numOutputNode, static_cast<long>(coverage.size())); ++i) {
//		out << '>' << buffer[i].second << endl;
		out << '>' << contigName[buffer[i].second] << endl;
		for (auto itr = coverage[buffer[i].second].begin(); itr != coverage[buffer[i].second].end(); ++itr)
			out << *itr << '\n';
	}

    out.close();
}


void PairedDBG::detectBreakpointBasedOnCoverage(const vector<unsigned>& physicalCoverage, const vector<unsigned>& diffCoverage, const long edgeLength, const double minCoverageRate, const double maxDiffCoverageRate, const long minCoverage, vector<char>& breakpoint)
{
	breakpoint.assign(physicalCoverage.size(), 0);

	if (breakpoint.size() <= 2*edgeLength)
		return;

	vector<unsigned> buffer(physicalCoverage.begin() + edgeLength, physicalCoverage.end() - edgeLength);
	std::partial_sort(buffer.begin(), buffer.begin() + buffer.size() / 2 + 1, buffer.end());
	const unsigned medianPhysicalCoverage = buffer[buffer.size() / 2];

/*
	buffer.assign(diffCoverage.begin() + edgeLength, diffCoverage.end() - edgeLength);
	std::partial_sort(buffer.begin(), buffer.begin() + buffer.size() / 2, buffer.end());
	const unsigned medianDiffCoverage = buffer[buffer.size() / 2];
*/

	for (long i = edgeLength; i < physicalCoverage.size() - edgeLength; ++i) {
		if (physicalCoverage[i] < minCoverageRate*medianPhysicalCoverage && diffCoverage[i] > maxDiffCoverageRate*physicalCoverage[i]) {
			breakpoint[i] = 1;
		}
	}


}


void PairedDBG::deleteShortRunOfBreakpoints(const long minRunLength, vector<char>& breakpoint)
{
    long i  = 0;
    while (i < breakpoint.size()) {
        while (breakpoint[i] == 0 && i < breakpoint.size())
            ++i;

        long j = i;
        while (breakpoint[i] == 1 && i < breakpoint.size())
			++i;

		if (i - j < minRunLength)
			std::fill(breakpoint.begin() + j, breakpoint.begin() + i, 0);
    }
}


bool PairedDBG::detectContigBoundaryBreakpoints(const long edgeLength, const GraphNode &targetNode, const vector<char>& baseBreakpoint, vector<char>& contigBreakpoint)
{
	contigBreakpoint.assign(targetNode.numContig + 1, 0);
	contigBreakpoint.back() = 1;

	bool isBroken = false;
	for (long i = 1; i < targetNode.numContig; ++i) {
		long checkStart = std::max(targetNode.contig[i-1].end - std::min(edgeLength, (targetNode.contig[i-1].end - targetNode.contig[i-1].start) / 2), 0L);
		long checkEnd = std::min(targetNode.contig[i].start + std::min(edgeLength, (targetNode.contig[i].end - targetNode.contig[i].start) / 2), static_cast<long>(baseBreakpoint.size()));

		if (std::find(baseBreakpoint.begin() + checkStart, baseBreakpoint.begin() + checkEnd, 1) != baseBreakpoint.begin() + checkEnd) {
			contigBreakpoint[i] = 1;
			isBroken = true;
		}
	}

	return isBroken;
}


void PairedDBG::node2GapFlagsUnmappableContig(const GraphNode &node, vector<char> &ret)
{
	static const char BASE2FLAG[] = {0, 0, 0, 0, 1};

    ret.assign(node.length, 1);

	for (long i = 0; i < node.numContig; ++i) {
		long contigIndex = id2Index(node.contig[i].id);
		if (contigPositionInScaffold[contigIndex].size() != 1) //edited by ouchi (contigPositionInScaffold)
			continue;

		if (node.contig[i].id > 0) {
			for (long j = 0; j < contig[contigIndex].length; ++j)
				ret[node.contig[i].start + j] = BASE2FLAG[static_cast<unsigned>(contig[contigIndex].base[j])];
		}
		else {
			for (long j = 0; j < contig[contigIndex].length; ++j)
				ret[node.contig[i].start + j] =  BASE2FLAG[static_cast<unsigned>(contig[contigIndex].base[contig[contigIndex].length - j - 1])];
		}
	}
}


void PairedDBG::divideNodeBasedOnCoverage(const vector<vector<unsigned> >& physicalCoverage, const vector<vector<unsigned> >& diffCoverage, const bool bubbleFlag, const long numThread)
{
	const long EDGE_LENGTH = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
	const double MIN_COVERAGE_RATE = 0.1;
	const double MAX_DIFF_COVERAGE_RATE = 10.0;
	const unsigned MIN_COVERAGE = 10;
//	const unsigned MIN_BREAKPOINT_RUN_LENGTH = 100;
	const long MIN_OVERLAP_TO_JOIN = 32;

	if (bubbleFlag)
		setOppositeBubbleNodeIDAndStateForEachNode();

    long defaultMinOverlap = this->minOverlap;
    this->minOverlap = MIN_OVERLAP_TO_JOIN;

	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		this->contigPreviousParentNodeID.resize(this->contig.size(), 0);

	long numDivided = 0;
	long numNewNode = 0;
	long newContigPoolSize = 0;
	vector<FILE*> scaffoldFP(numThread);
	vector<vector<std::pair<int, int> > > unlinkListBuffer(numThread);

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(static, 1) reduction(+: numDivided, numNewNode, newContigPoolSize)
	for (unsigned threadIndex = 0; threadIndex < numThread; ++threadIndex) {
		scaffoldFP[threadIndex] = platanus::makeTemporaryFile();
		vector<char> baseBreakpoint;
		vector<char> contigBreakpoint;

		for (long i = threadIndex; i < numNode; i += numThread) {
			if (node[i].numContig == 1) {
				fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP[threadIndex]);
				fwrite(&(node[i].contig[0]), sizeof(ScaffoldPart), 1, scaffoldFP[threadIndex]);
				++numNewNode;
				++newContigPoolSize;
				continue;
			}
			else if (physicalCoverage[i].size() <= 2*EDGE_LENGTH || (bubbleFlag && (node[i].state & (DBG_PRIMARY_BUBBLE | DBG_SECONDARY_BUBBLE)))) {
				fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP[threadIndex]);
				for (auto itr = node[i].contig.begin(); itr != node[i].contig.end(); ++itr)
					fwrite(&(*itr), sizeof(ScaffoldPart), 1, scaffoldFP[threadIndex]);
				++numNewNode;
				newContigPoolSize += node[i].numContig;
				continue;
			}

			detectBreakpointBasedOnCoverage(physicalCoverage[i], diffCoverage[i], EDGE_LENGTH, MIN_COVERAGE_RATE, MAX_DIFF_COVERAGE_RATE, MIN_COVERAGE, baseBreakpoint);
//			deleteShortRunOfBreakpoints(MIN_BREAKPOINT_RUN_LENGTH, baseBreakpoint);

			bool isDivided = detectContigBoundaryBreakpoints(EDGE_LENGTH, node[i], baseBreakpoint, contigBreakpoint);

			if (isDivided) {
				++numDivided;

				long j = 0;
				while (j < node[i].numContig) {
					long start = node[i].contig[j].start;
					long k = j;
					while (contigBreakpoint[j + 1] == 0 && j < node[i].numContig - 1 && node[i].contig[j + 1].end >= start) {
						node[i].contig[j].start -= start;
						node[i].contig[j].end -= start;
						++j;
					}
					node[i].contig[j].start -= start;
					node[i].contig[j].end -= start;
					++j;

					long numContig = j - k;
					fwrite(&numContig, sizeof(long), 1, scaffoldFP[threadIndex]);
					for (long m = k; m < j; ++m) {
						fwrite(&(node[i].contig[m]), sizeof(ScaffoldPart), 1, scaffoldFP[threadIndex]);

						int contigIndex1 = 0; //edited by ouchi (contigPositionInScaffold)
						if (contigPositionInScaffold[id2Index(node[i].contig[m].id)].size() != 1) //edited by ouchi (contigPositionInScaffold)
							contigIndex1 = contigPositionInScaffold[id2Index(node[i].contig[m].id)][0].id; //edited by ouchi (contigPositionInScaffold)
						if (contigIndex1 == 0)
							continue;

						for (long n = 0; n < m; ++n) {
							int contigIndex2 = 0; //edited by ouchi (contigPositionInScaffold)
							if (contigPositionInScaffold[id2Index(node[i].contig[n].id)].size() != 1) //edited by ouchi (contigPositionInScaffold)
								contigIndex2 = id2Index(contigPositionInScaffold[id2Index(node[i].contig[n].id)][0].id); //edited by ouchi (contigPositionInScaffold)
							if (contigIndex2 == 0)
								continue;

							unlinkListBuffer[threadIndex].push_back(contigIndex1 < contigIndex2 ? std::make_pair(contigIndex1, contigIndex2) : std::make_pair(contigIndex2, contigIndex1));
						}

						for (long n = m + numContig; n < node[i].numContig; ++n) {
							int contigIndex2 = 0; //edited by ouchi (contigPositionInScaffold)
							if (contigPositionInScaffold[id2Index(node[i].contig[n].id)].size() != 1) //edited by ouchi (contigPositionInScaffold)
								contigIndex2 = id2Index(contigPositionInScaffold[id2Index(node[i].contig[n].id)][0].id); //edited by ouchi (contigPositionInScaffold)
							if (contigIndex2 == 0)
								continue;

							unlinkListBuffer[threadIndex].push_back(contigIndex1 < contigIndex2 ? std::make_pair(contigIndex1, contigIndex2) : std::make_pair(contigIndex2, contigIndex1));
						}
					}

					++numNewNode;
					newContigPoolSize += numContig;

				}
			}
			else {
				fwrite(&(node[i].numContig), sizeof(long), 1, scaffoldFP[threadIndex]);
				for (auto itr = node[i].contig.begin(); itr != node[i].contig.end(); ++itr)
					fwrite(&(*itr), sizeof(ScaffoldPart), 1, scaffoldFP[threadIndex]);
				++numNewNode;
				newContigPoolSize += node[i].numContig;
			}
		}
	}


    for (long i = 1; i < numThread; ++i) {
		char tmpChar;
        rewind(scaffoldFP[i]);
        while (fread(&tmpChar, sizeof(char), 1, scaffoldFP[i]))
            putc(tmpChar, scaffoldFP[0]);
        fclose(scaffoldFP[i]);

		for (auto itr = unlinkListBuffer[i].begin(); itr != unlinkListBuffer[i].end(); ++itr)
			this->contigUnlinkSet.insert(*itr);
		unlinkListBuffer[i].clear();
    }

    this->remake(numNewNode, newContigPoolSize, scaffoldFP[0]);
    fclose(scaffoldFP[0]);

    this->minOverlap = defaultMinOverlap;

    std::cerr << "NUM_DIVIDED_ERROR_CANDIDATES_BASE_LEVEL = " << numDivided << endl;
}


void PairedDBG::storeBreakpointBasedOnCoverage(const vector<vector<unsigned> >& physicalCoverage, const vector<vector<unsigned> >& diffCoverage, const long numThread)
{
	const long EDGE_LENGTH = (*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
	const double MIN_COVERAGE_RATE = 0.5;
	const double MAX_DIFF_COVERAGE_RATE = 1.0;
	const unsigned MIN_COVERAGE = 10;

	long numDivided = 0;
	nodeBreakpointPosition.resize(numNode);

	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		this->contigPreviousParentNodeID.resize(this->contig.size(), 0);

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(static, 1) reduction(+: numDivided)
	for (unsigned threadIndex = 0; threadIndex < numThread; ++threadIndex) {
		vector<char> baseBreakpoint;
		for (long i = threadIndex; i < numNode; i += numThread) {
			if (node[i].numContig != 1 || physicalCoverage[i].size() <= 2*EDGE_LENGTH)
				continue;

			detectBreakpointBasedOnCoverage(physicalCoverage[i], diffCoverage[i], EDGE_LENGTH, MIN_COVERAGE_RATE, MAX_DIFF_COVERAGE_RATE, MIN_COVERAGE, baseBreakpoint);

			bool divideFlag = false;
			for (long j = EDGE_LENGTH; j < baseBreakpoint.size() - EDGE_LENGTH; ++j) {
				if (baseBreakpoint[j] == 1) {
					nodeBreakpointPosition[i].push_back(j);
					divideFlag = true;
				}
			}

			if (divideFlag) {
				++numDivided;
			}
		}
	}

    std::cerr << "NUM_DIVIDED_ERROR_CANDIDATES_BASE_LEVEL = " << numDivided << endl;
}

//added by ouchi
void PairedDBG::divideNodeBasedOnHiC(long numThread)
{
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        vector<vector<int> > contactmap;
        makeContactmap(nodeID+1, contactmap);
        vector<long> locate;
        detectContactChange(contactmap,0, contactmap.size(), locate);
        if (locate.size() > 0)
            std::cerr << "divideNode: " << nodeID + 1 << std::endl;
        for (long i = 0; i < locate.size(); ++i) {
            std::cerr << locate[i] << "\t";
        }
        if (locate.size() > 0) std::cerr << std::endl;
    }
}
//

//added by ouchi
bool PairedDBG::checkHiCLinkBetweenNodePair(const long nodeID1, const char direction, const long nodeID2) {
    vector<vector<int> > contactmap;
    long boundary = makeContactmap(nodeID1, direction, nodeID2, contactmap);
    if (boundary == -1) return 1;
    std::array<double, 6> result;
    calcSeparation(contactmap, 0, boundary, contactmap.size(), result, true);
    return (result[0]/result[1] + result[3]/result[4] < HiC_CONTACTMAP_S_THRESHOLD || result[2] + result[5] < HiC_CONTACTMAP_T_THRESHOLD);
}
//

//added by ouchi
//////////////////////////////////////////////////////////////////////////////////////
// calc IdealContact by HiC links inside node
//////////////////////////////////////////////////////////////////////////////////////
void PairedDBG::calcIdealContact(const long numThread)
{
    makeHiCLink(numThread);

    idealcontact.clear();

    vector<int> contact;
    vector<int> num;

    omp_set_num_threads(numThread);
//    #pragma omp parallel for schedule(dynamic)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        vector<vector<int> > contactmap;
        this->makeContactmap(nodeID+1, contactmap);
//        #pragma omp critical (calc_contact)
        {
            if (contactmap.size() > contact.size()) {
                contact.resize(contactmap.size());
                num.resize(contactmap.size());
            }
            for (unsigned int i = 0; i < contactmap.size(); ++i) {
                for (int j = i; j < contactmap.size(); ++j) {
                    contact[i] += contactmap[j][j-i];
                    ++num[i];
                }
            }
        }
    }

    idealcontact.resize(contact.size());
    for (unsigned int i = 0; i < contact.size(); ++i) {
        idealcontact[i] = static_cast<double>(contact[i]) / num[i];
		std::cerr << idealcontact[i] << std::endl; //added by ouchi for test
    }
}
//

//added by ouchi
//////////////////////////////////////////////////////////////////////////////////////
// make HiC Contactmap of nodeID
//////////////////////////////////////////////////////////////////////////////////////
void PairedDBG::makeContactmap(const long nodeID, vector<vector<int> > &contactmap)
{
    contactmap.clear();
    const long binsize = 100000;
    GraphNode &node = this->node[id2Index(nodeID)];
    long numbin = node.length/binsize;
    if (numbin < 1)
        return;

    long lower = 0;
    long upper = node.length - node.length % binsize;

    contactmap.resize(numbin);
    for (long i = 0; i < numbin; ++i) {
        contactmap[i].resize(numbin);
    }

    for (long HiCLinkID = 0; HiCLinkID < node.HiCLinks.size(); ++HiCLinkID) {
        const HiCLink &hiclink = node.HiCLinks[HiCLinkID];
        if (abs(hiclink.end) == abs(nodeID)) {
            if (lower<=hiclink.offset1 && hiclink.offset1<upper && lower<=hiclink.offset2 && hiclink.offset2<upper)
                ++contactmap[(hiclink.offset1-lower)/binsize][(hiclink.offset2-lower)/binsize];
        }
    }

	std::cerr << "contactmap " << nodeID << std::endl;
	for (long i = 0; i < contactmap.size(); ++i) {
		for (long j = 0; j < contactmap[i].size(); ++j) {
			std::cerr << contactmap[i][j] << "\t";
		}
		std::cerr << std::endl;
	}

}
//

//added by ouchi
//////////////////////////////////////////////////////////////////////////////////////
// make HiC Contactmap between nodeID1 and nodeID2
//////////////////////////////////////////////////////////////////////////////////////
long PairedDBG::makeContactmap(const long nodeID1, const char direction, const long nodeID2, vector<vector<int> > &contactmap)
{
    const long binsize = 100000;
    const GraphNode &node1 = this->node[id2Index(nodeID1)];
    const GraphNode &node2 = this->node[id2Index(nodeID2)];

    if (node1.length/binsize < 1 || node2.length/binsize < 1)
        return -1;

    long numbin = node1.length/binsize + node2.length/binsize;
    long lower, upper, boundary;
    if (direction > 0) {
        lower = node1.length % binsize;
        upper = node1.length + node2.length - node2.length % binsize;
        boundary = node1.length/binsize;
    } else {
        lower = node2.length % binsize;
        upper = node2.length + node1.length - node1.length % binsize;
        boundary = node2.length / binsize;
    }


    contactmap.resize(numbin);
    for (long i = 0; i < numbin; ++i) {
        contactmap[i].resize(numbin);
    }

    long offset1, offset2;
    for (long HiCLinkID = 0; HiCLinkID < node1.HiCLinks.size(); ++HiCLinkID) {
        const HiCLink &hiclink = node1.HiCLinks[HiCLinkID];
        if (abs(hiclink.end) == abs(nodeID1)) {
            offset1 = hiclink.offset1;
            offset2 = hiclink.offset2;
            if (direction < 0) {
                offset1 += node2.length;
                offset2 += node2.length;
            }
        } else if (abs(hiclink.end) == abs(nodeID2)) {
            if (nodeID2 > 0)
                offset2 = hiclink.offset2;
            else
                offset2 = node2.length - hiclink.offset2 - 1;
            if (direction > 0) {
                offset1 = hiclink.offset1;
                offset2 += node1.length;
            } else {
                offset1 = hiclink.offset1 + node2.length;
            }
        } else {
            continue;
        }
        if (lower<=offset1 && offset1<upper && lower<=offset2 && offset2<upper) {
            ++contactmap[(offset1-lower)/binsize][(offset2-lower)/binsize];
        }
    }
    for (long HiCLinkID = 0; HiCLinkID < node2.HiCLinks.size(); ++HiCLinkID) {
        const HiCLink &hiclink = node2.HiCLinks[HiCLinkID];
        if (abs(hiclink.end) == abs(nodeID2)) {
            if (nodeID2 > 0) {
                offset1 = hiclink.offset1;
                offset2 = hiclink.offset2;
            } else {
                offset1 = node2.length - hiclink.offset1 - 1;
                offset2 = node2.length - hiclink.offset2 - 1;
            }
            if (direction > 0) {
                offset1 += node1.length;
                offset2 += node1.length;
            }
        } else if (abs(hiclink.end) == abs(nodeID1)) {
            if (nodeID2 > 0) {
                offset1 = hiclink.offset1;
            } else {
                offset1 = node2.length - hiclink.offset1 - 1;
            }
            if (direction > 0) {
                offset2 = hiclink.offset2;
                offset1 += node1.length;
            } else {
                offset2 = hiclink.offset2 + node2.length;
            }
        } else {
            continue;
        }
        if (lower<=offset1 && offset1<upper && lower<=offset2 && offset2<upper) {
            ++contactmap[(offset1-lower)/binsize][(offset2-lower)/binsize];
        }
    }

	/*std::cerr << "makeContactmap\t" << nodeID1 << "\t" << direction << "\t" << nodeID2 <<std::endl;
	for (long i = 0; i < contactmap.size(); ++i) {
		if (i == boundary) {
			for (long j = 0; j < contactmap[i].size(); ++j) {
				std::cerr << "_\t";
			}
			std::cerr << std::endl;
		}
		for (long j = 0; j < contactmap[i].size(); ++j) {
			if (j == boundary) std::cerr << "|\t";
			std::cerr << contactmap[i][j] << "\t";
		}
		std::cerr << std::endl;
	}*/


    return boundary;
}
//

//added by ouchi
//////////////////////////////////////////////////////////////////////////////////////
// calc separation degree (downS,downS0,downT,upS,upS0,upT)
//////////////////////////////////////////////////////////////////////////////////////
void PairedDBG::calcSeparation(const vector<vector<int> > &contactmap, const long n1, const long n2, const long n3, std::array<double, 6> &result, bool divided = false)
{
    double peakk; // = percentilvalue;
    peakk = 100000;

    result.fill(0);

    vector<vector<int> > in,out;
    long num_under_av;
    long num_out;
    double sum;
    long i,j,k;

    //calculation under side
    in.clear(); out.clear();
    if (n2-n1 == 1) {
        double av_in_k = contactmap[n1][n1];
        double av_out_k = contactmap[n1+1][n1];
        result[0] = (av_in_k-av_out_k) * (av_in_k-av_out_k) / 4;
        result[1] = (av_in_k-0) * (av_in_k-0) / 4;
        if (av_out_k < av_in_k || av_out_k == 0)
            result[2] = 1;
    } else {
        in.resize(n2-n1-1); out.resize(n2-n1-1);
        for (i = 1; i < n2-n1; ++i) {
            for (j = 0, k = i; j < n2-n1; ++j,++k) {
                if (n1+k >= n3) break;
                int t = contactmap[n1+k][n1+j];
                if (t > peakk) t = peakk;
                if (n1+k < n2)
                    in[i-1].push_back(t);
                else
                    out[i-1].push_back(t);
            }
        }

        num_under_av = 0;
        num_out = 0;
        if (divided)
            k = 1;
        else
            k = 2;
        for (; k < n2-n1; ++k) {
            sum = 0;
            for (i = 0; i < in[k-1].size(); ++i) {
                sum += in[k-1][i];
            }
            double av_in_k = sum / in[k-1].size();
            sum = 0;
            for (i = 0; i < out[k-1].size(); ++i) {
                sum += out[k-1][i];
                if (out[k-1][i] < av_in_k || out[k-1][i] == 0)
                    ++num_under_av;
                ++num_out;
            }
            double av_out_k = sum / out[k-1].size();
            if (av_in_k > av_out_k)
                result[0] += in[k-1].size() * out[k-1].size() * (av_in_k-av_out_k) * (av_in_k-av_out_k) / (in[k-1].size()+out[k-1].size()) / (in[k-1].size()+out[k-1].size());
            if (av_in_k > 0)
                result[1] += in[k-1].size() * out[k-1].size() * av_in_k * av_in_k / (in[k-1].size()+out[k-1].size()) / (in[k-1].size()+out[k-1].size());
        }
        result[2] = double(num_under_av) / num_out;
    }

    //calculation upper side
    in.clear(); out.clear();
    if (n3-n2 == 1) {
        double av_in_k = contactmap[n2][n2];
        double av_out_k = contactmap[n2][n2-1];
        result[3] = (av_in_k-av_out_k) * (av_in_k-av_out_k) / 4;
        result[4] = (av_in_k-0) * (av_in_k-0) / 4;
        if (av_out_k < av_in_k || av_out_k == 0)
            result[5] = 1;
    } else {
        in.resize(n3-n2-1); out.resize(n3-n2-1);
        for (i = 1; i < n3-n2; ++i) {
            for (j = 0, k = i; j < n3-n2; ++j,--k) {
                if (n2-k < n1) continue;
                int t = contactmap[n2+j][n2-k];
                if (t > peakk) t = peakk;
                if (k <= 0)
                    in[i-1].push_back(t);
                else
                    out[i-1].push_back(t);
            }
        }

        num_under_av = 0;
        num_out = 0;
        if (divided)
            k = 1;
        else
            k = 2;
        for (; k < n3-n2; ++k) {
            sum = 0;
            for (i = 0; i < in[k-1].size(); ++i) {
                sum += in[k-1][i];
            }
            double av_in_k = sum / in[k-1].size();
            sum = 0;
            for (i = 0; i < out[k-1].size(); ++i) {
                sum += out[k-1][i];
                if (out[k-1][i] < av_in_k || out[k-1][i] == 0)
                    ++num_under_av;
                ++num_out;
            }
            double av_out_k = sum / out[k-1].size();
            if (av_in_k > av_out_k)
                result[3] += in[k-1].size() * out[k-1].size() * (av_in_k-av_out_k) * (av_in_k-av_out_k) / (in[k-1].size()+out[k-1].size()) / (in[k-1].size()+out[k-1].size());
            if (av_in_k > 0)
                result[4] += in[k-1].size() * out[k-1].size() * av_in_k * av_in_k / (in[k-1].size()+out[k-1].size()) / (in[k-1].size()+out[k-1].size());
        }
        result[5] = double(num_under_av) / num_out;
    }

}
//

//added by ouchi
template <typename T>
void PairedDBG::detectPeak(const vector<T> &v, long delta, long bin, long maxh, vector<long> &locate)
{
    locate.clear();

    long mn = v[0]; long mx = v[0];
    long mnpos = -1; long mxpos = -1;
    bool lookformax = true;

    for (long i = 0; i < v.size(); ++i) {
        T tmp = v[i];
        if (tmp > mx) {
            mx = tmp;
            mxpos = i;
        }
        if (tmp < mn) {
            mn = tmp;
            mnpos = i;
        }

        if (lookformax) {
            if (tmp < mx -delta) {
                if (i - mxpos > bin) {
                    if (mx > maxh) {
                        locate.push_back(mxpos);
                    }
                    mn = tmp;
                    mnpos = i;
                    lookformax = false;
                }
            }
        } else {
            if (tmp > mn + delta) {
                if (i - mnpos > bin) {
                    mx = tmp;
                    mxpos = i;
                    lookformax = true;
                }
            }
        }
    }
}
//

//added by ouchi
//////////////////////////////////////////////////////////////////////////////////////
// detect Contactpattern change point
//////////////////////////////////////////////////////////////////////////////////////
void PairedDBG::detectContactChange(const vector<vector<int> > &contactmap, const long start, const long end, vector<long> &locate)
{
    long b = std::min(10, static_cast<int>(idealcontact.size() - 1)); //correction_L
    double peakk = 100000; // = percentilvalue;

    if (end-start <=1)
        return;

    long i,j,k,n1,n2,n3;

    double max_change = 0;
    vector<double> contact_change(contactmap.size()+4*b, 0);
    for (i = start-b; i < end+b+1; ++i) {
        vector<double> tmp(b);
        for (j = 0; j < b; ++j) {
            for (k = 0; k < j+1; ++k) {
                double t = 0;
                if (i+k < end && i-1-j+k >= start)
                    t = contactmap[i+k][i-1-j+k];
                if (t < peakk)
                    tmp[j] += t;
                else
                    tmp[j] += peakk;
            }
        }
        for (j = 0; j < b ; ++j) {
            if (tmp[j]/(j+1) < idealcontact[j+1])
                contact_change[i+b] += (tmp[j]/(j+1) - idealcontact[j+1]) * (tmp[j]/(j+1) - idealcontact[j+1]) / idealcontact[j+1];
        }
        if (contact_change[i+b] > max_change)
            max_change = contact_change[i+b];
    }

    vector<long> peaks;
    std::array<double, 6> result;

    double maxod = 0;
    long maxj = -1;
    unsigned int num_contactpeak = 0;
    for (j = static_cast<long>(max_change); j > 0; --j) {
        detectPeak(contact_change, 0, 0, j, peaks);
        if (num_contactpeak == peaks.size() || peaks.size() <= 2) continue;
        num_contactpeak = peaks.size();
        peaks.erase(peaks.begin());
        peaks.pop_back();
        double od = 0;
        for (i = 0; i < peaks.size(); ++i) {
            n1 = start; n2 = peaks[i]-b; n3 = end;
            if (i >= 1) n1 = peaks[i-1] - b;
            if (i+1 < peaks.size()) n3 = peaks[i+1] - b;
            calcSeparation(contactmap, n1, n2, n3, result);
            od += result[0] + result[3];
        }
        od = od / peaks.size();
        if (od > maxod) {
            maxod = od; maxj = j;
        }
    }

    if (maxj == -1) return;

    detectPeak(contact_change, 0, 0, maxj, peaks);
    peaks.erase(peaks.begin());
    peaks.pop_back();
    for (i = 0; i < peaks.size(); ++i) {
        peaks[i] -= b;
    }
    for (i = 0; i < peaks.size();) {
        n1 = start; n2 = peaks[i]; n3 = end;
        if (i >= 1) n1 = peaks[i-1];
        if (i+1 < peaks.size()) n3 = peaks[i+1];
        calcSeparation(contactmap, n1 ,n2 ,n3, result);
        if (result[0]/result[1] + result[3]/result[4] < HiC_CONTACTMAP_S_THRESHOLD || result[2] + result [5] < HiC_CONTACTMAP_T_THRESHOLD) {
            peaks.erase(peaks.begin()+i);
            continue;
        }
        detectContactChange(contactmap, n1, n2, locate);
        if (i+1 == peaks.size())
            detectContactChange(contactmap, n2, n3, locate);
        ++i;
    }

    locate.insert(locate.end(), peaks.begin(), peaks.end());
    std::sort(locate.begin(), locate.end());

    for (i = 0; i < locate.size();) {
        n1 = start; n2 = locate[i]; n3 = end;
        if (i >= 1) n1 = locate[i-1];
        if (i+1 < locate.size()) n3 = locate[i+1];
        calcSeparation(contactmap, n1, n2, n3, result);
        if (result[0]/result[1] + result[3]/result[4] < HiC_CONTACTMAP_S_THRESHOLD || result[2] + result[5] < HiC_CONTACTMAP_T_THRESHOLD) {
            locate.erase(locate.begin()+i);
            continue;
        }
        ++i;
    }

}
//

//added by ouchi
//////////////////////////////////////////////////////////////////////////////////////
// calculate PercentilValue
//////////////////////////////////////////////////////////////////////////////////////
void PairedDBG::calcPercentilValue(const long percent, const long numThread)
{
    long binsize = 100000;

    cerr << "calcParcentilValue" << endl;

    makeHiCLink(numThread);

    vector<vector<int> > contact(numThread);
    //# pragma omp parallel for schedule(static, 1)
    for (long t = 0; t < numThread; ++t) {
        for (long nodeID = t; nodeID < numNode; nodeID += numThread) {
            vector<vector<int> > contacts(numNode);
            for (long HiCLinkID = 0; HiCLinkID < node[id2Index(nodeID)].HiCLinks.size(); ++HiCLinkID) {
                const HiCLink &hiclink = node[id2Index(nodeID)].HiCLinks[HiCLinkID];
                contacts[id2Index(hiclink.end)].push_back(hiclink.offset2/binsize);
            }
            for (unsigned int i = 0; i < contacts.size(); ++i) {
                std::sort(contacts[i].begin(), contacts[i].end());
                int tmp = 0; int num = 0;
                for (unsigned int j = 0; j < contacts[i].size(); ++j) {
                    if (contacts[i][j] == tmp) {
                        ++num;
                    } else {
                        if (num > 0)
                            contact[t].push_back(num);
                        tmp = contacts[i][j];
                        num = 1;
                    }
                }
            }
        }
    }
	cerr << "finish" << endl;

    vector<int> mergedcontact;
	cerr << "merge" << endl;
    mergeAndClearMultiThreadedVector(contact, mergedcontact);
	cerr << "sort" << endl;
    std::partial_sort(mergedcontact.begin(), mergedcontact.begin() + mergedcontact.size() * percent / 100 + 1, mergedcontact.end());
    //percentilvalue = mergedcontact[mergedcontact.size() * percent / 100];
    //cerr << percent << "PERCENTILVALUE=" << percentilvalue << endl;

}
//

//added by ouchi
void PairedDBG::deleteGapRegionFromContactmap(vector<vector<int> > &contactmap, long &boundary) {
    vector<vector<int> > originContactmap = contactmap;
    contactmap.clear();
    long numDelete = 0;
    for (long i = 0; i < originContactmap.size(); ++i) {
        if (originContactmap[i][i] == 0) {
            if (i < boundary)
                ++numDelete;
            continue;
        }
        contactmap.push_back({});
        for (long j = 0; j < originContactmap[i].size(); ++j) {
            if (originContactmap[j][j] == 0)
                continue;
            contactmap.back().push_back(originContactmap[i][j]);
        }
    }
    boundary = boundary - numDelete;
}
//

//added by ouchi
double PairedDBG::getHiCLinkScoreBetweenNodePair(const long nodeID1, const char targetDirection, const long nodeID2, const long leftLimit, const long rightLimit)
{
    const GraphNode &node1 = node[id2Index(nodeID1)];
    const GraphNode &node2 = node[id2Index(nodeID2)];

    double score = 0; long offset1 = 0; long offset2 = 0;
    for (long HiCLinkID = 0; HiCLinkID < node1.HiCLinks.size(); ++HiCLinkID) {
        const HiCLink &hiclink = node1.HiCLinks[HiCLinkID];
        if (abs(hiclink.end) == abs(nodeID2)) {
            if (nodeID1 > 0)
                offset1 = node1.length - hiclink.offset1 - 1;
            else
                offset1 = hiclink.offset1;
            if (nodeID2 > 0)
                offset2 = hiclink.offset2;
            else
                offset2 = node2.length - hiclink.offset2 - 1;
            if (offset1 <= leftLimit && offset2 <= rightLimit)
                score += static_cast<double>(leftLimit - offset1) / leftLimit * (rightLimit - offset2) / rightLimit;
        }
    }
    for (long HiCLinkID = 0; HiCLinkID < node2.HiCLinks.size(); ++HiCLinkID) {
        const HiCLink &hiclink = node2.HiCLinks[HiCLinkID];
        if (abs(hiclink.end) == abs(nodeID1)) {
            if (nodeID1 > 0)
                offset1 = node1.length - hiclink.offset2 - 1;
            else
                offset1 = hiclink.offset2;
            if (nodeID2 > 0)
                offset2 = hiclink.offset1;
            else
                offset2 = node2.length - hiclink.offset1 - 1;
            if (offset1 <= leftLimit && offset2 <= rightLimit)
                score += static_cast<double>(leftLimit - offset1) / leftLimit * (rightLimit - offset2) / rightLimit;
        }
   }
   return score;
}
//

//added by ouchi
void PairedDBG::calculateHiCReadPhysicalCoverage(const long nodeID, vector<unsigned>& physicalCoverage)
{
    const GraphNode targetnode = node[id2Index(nodeID)];
    physicalCoverage.resize(targetnode.length);
    int offset1; int offset2;
    for (long HiCLinkID = 0; HiCLinkID < targetnode.HiCLinks.size(); ++HiCLinkID) {
        const HiCLink hiclink = targetnode.HiCLinks[HiCLinkID];
        if (abs(hiclink.end) == abs(nodeID)) {
            offset1 = hiclink.offset1;
            offset2 = hiclink.offset2;
            std::pair<int, int> range = std::minmax(
                std::min(std::max(offset1, 0), static_cast<int>(targetnode.length-1)),
                std::min(std::max(offset2, 0), static_cast<int>(targetnode.length-1))
            );
            for (long i = range.first; i <= range.second; ++i) {
                ++(physicalCoverage[i]);
            }
        }
    }
}
//

//added by ouchi
bool PairedDBG::searchPathBetweenNodePair(const long nodeID1, const char targetDirection, const long nodeID2, const long limitLength)
{
    const GraphNode &sourceNode = node[id2Index(nodeID1)];
    for (long edgeIndex = 0; edgeIndex < sourceNode.numEdge; ++edgeIndex) {
        const GraphEdge &edge = sourceNode.edge[edgeIndex];
        if (edge.direction == targetDirection && edge.numLink >= minLink && edge.length < limitLength) {
            if (edge.end == nodeID2) {
                return 1;
            } else {
                if (!(node[id2Index(edge.end)].state & SC_DEL)) {
                    if (edge.length + node[id2Index(edge.end)].length < limitLength) {
                        searchPathBetweenNodePair(edge.end,  targetDirection * sign(edge.end), nodeID2, limitLength - edge.length - node[id2Index(edge.end)].length);
                    }
                }
            }
        }
    }
    return 0;
}
//

//added by ouchi
void PairedDBG::getConflictRegionBetweenNodePair(const GraphEdge &edge1, const GraphEdge &edge2, const GraphNode &node1, const GraphNode &node2, vector<char> &conflictpoint1, vector<char> &conflictpoint2)
{
    conflictpoint1.resize(node1.length, 1);
    conflictpoint2.resize(node2.length, 1);
    char sign1 = edge1.direction * sign(edge1.end);
    char sign2 = edge2.direction * sign(edge2.end);

    if (edge1.length < 0) {
        for (long i = 0; i < -edge1.length && i < node1.length; ++i) {
            if (sign1 > 0)
                conflictpoint1[i] = 0;
            else
                conflictpoint1[node1.length-1-i] = 0;
        }
    }
    if (edge1.length < edge2.length) {
        for (long i = 0; i < edge2.length - edge1.length && i < node1.length; ++i) {
            if (sign1 > 0)
                conflictpoint1[i] = 0;
            else
                conflictpoint1[node1.length-1-i] = 0;
        }
    }
    if (edge1.length + node1.length > edge2.length + node2.length) {
        for (long i = edge2.length + node2.length - edge1.length; i >= 0 && i < node1.length; ++i) {
            if (sign1 > 0)
                conflictpoint1[i] = 0;
            else
                conflictpoint1[node1.length-1-i] = 0;
        }
    }
    if (edge2.length < 0) {
        for (int i = 0; i < -edge2.length && i < node2.length; ++i) {
            if (sign2 > 0)
                conflictpoint2[i] = 0;
            else
                conflictpoint2[node2.length-1-i] = 0;
        }
    }
    if (edge2.length < edge1.length) {
        for (int i = 0; i < edge1.length - edge2.length && i < node2.length; ++i) {
            if (sign2 > 0)
                conflictpoint2[i] = 0;
            else
                conflictpoint2[node2.length-1-i] = 0;
        }
    }
    if (edge2.length + node2.length > edge1.length + node1.length) {
        for (long i = edge1.length + node1.length - edge2.length; i >= 0 && i< node2.length; ++i) {
            if (sign2 > 0)
                conflictpoint2[i] = 0;
            else
                conflictpoint2[node2.length-1-i] = 0;
        }
    }

    for (long contigIndex = 0; contigIndex+1 < node1.contig.size(); ++contigIndex) {
        long start = node1.contig[contigIndex].end;
        long end = node1.contig[contigIndex+1].start;
        for (long i = start; i < end && i < node1.length; ++i) {
            if (sign1 > 0)
                conflictpoint1[i] = 0;
            else
                conflictpoint1[node1.length-1-i] = 0;
            if (i-(edge2.length-edge1.length) >= 0 && i-(edge2.length-edge1.length) < node2.length) {
                if (sign2 > 0)
                    conflictpoint2[i-(edge2.length-edge1.length)] = 0;
                else
                    conflictpoint2[node2.length-1-(i-(edge2.length-edge1.length))] = 0;
            }
        }
    }
    for (long contigIndex = 0; contigIndex+1 < node2.contig.size(); ++contigIndex) {
        long start = node2.contig[contigIndex].end;
        long end = node2.contig[contigIndex+1].start;
        for (long i = start; i < end && i < node2.length; ++i) {
            if (sign2 > 0)
                conflictpoint2[i] = 0;
            else
                conflictpoint2[node2.length-1-i] = 0;
            if (i-(edge1.length-edge2.length) >=0 && i-(edge1.length-edge2.length) < node1.length) {
                if (sign1 > 0)
                    conflictpoint1[i-(edge1.length-edge2.length)] = 0;
                else
                    conflictpoint1[node1.length-1-(i-(edge1.length-edge2.length))] = 0;
            }
        }
    }

}
//

//added by ouchi
void PairedDBG::getHeteroMPLinkNum(long centerNodeIndex, const std::array<std::array<NodeIDWithGap ,2>, 2> &externalNodeID, std::array<double, 2> &sumLinkForHaplotype, std::array<double, 4> &MPScores)
{
    GraphNode &centerNode = this->node[centerNodeIndex];
    std::array<std::array<GraphNode, 2>, 2> externalNode;
    externalNode[0][0] = this->node[id2Index(externalNodeID[0][0].id)];
    externalNode[0][1] = this->node[id2Index(externalNodeID[0][1].id)];
    externalNode[1][0] = this->node[id2Index(externalNodeID[1][0].id)];
    externalNode[1][1] = this->node[id2Index(externalNodeID[1][1].id)];

    std::array<std::array<GraphEdge, 2>, 2> externalEdge;
    externalEdge[0][0].direction = -1;
    externalEdge[0][0].end = externalNodeID[0][0].id;
    externalEdge[0][0].length = externalNodeID[0][0].gap;
    externalEdge[0][1].direction = -1;
    externalEdge[0][1].end = externalNodeID[0][1].id;
    externalEdge[0][1].length = externalNodeID[0][1].gap;
    externalEdge[1][0].direction = 1;
    externalEdge[1][0].end = externalNodeID[1][0].id;
    externalEdge[1][0].length = externalNodeID[1][0].gap;
    externalEdge[1][1].direction = 1;
    externalEdge[1][1].end = externalNodeID[1][1].id;
    externalEdge[1][1].length = externalNodeID[1][1].gap;

    std::array<std::array<vector<char>, 2>, 2> conflictpoints;
    getConflictRegionBetweenNodePair(externalEdge[0][0], externalEdge[0][1], this->node[id2Index(externalNodeID[0][0].id)], this->node[id2Index(externalNodeID[0][1].id)], conflictpoints[0][0], conflictpoints[0][1]);
    getConflictRegionBetweenNodePair(externalEdge[1][0], externalEdge[1][1], this->node[id2Index(externalNodeID[1][0].id)], this->node[id2Index(externalNodeID[1][1].id)], conflictpoints[1][0], conflictpoints[1][1]);

    long offset1 = 0; long offset2 = 0;
    for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
        for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
            for (long libraryID = 0; libraryID < (*allLibraryMT).size(); ++libraryID) {
                for (long MPLinkID = 0; MPLinkID < externalNode[0][leftIndex].MPLinks[libraryID].size(); ++MPLinkID) {
                    const MPLink &mpLink = externalNode[0][leftIndex].MPLinks[libraryID][MPLinkID];
                    if (mpLink.direction * externalNodeID[0][leftIndex].id > 0) {
                        if (mpLink.end == externalNodeID[1][rightIndex].id * sign(externalNodeID[0][leftIndex].id)) {
                            offset1 = 0; offset2 = 0;
                            offset1 = std::min(std::max(mpLink.offset1, offset1), externalNode[0][leftIndex].length-1);
                            offset2 = std::min(std::max(mpLink.offset2, offset2), externalNode[1][rightIndex].length-1);
                            if (conflictpoints[0][leftIndex][offset1] == 0 || conflictpoints[1][rightIndex][offset2] == 0)
                                continue;
                            if (externalNodeID[0][leftIndex].id > 0)
                                offset1 = externalNode[0][leftIndex].length - offset1 - 1;
                            if (externalNodeID[1][rightIndex].id < 0)
                                offset2 = externalNode[1][rightIndex].length - offset2 - 1;
                            if (centerNode.length + externalEdge[0][leftIndex].length + externalEdge[1][rightIndex].length + offset1 + offset2 > (*allLibraryMT)[libraryID][0].getAverageInsSize() + tolerence
                             || centerNode.length + externalEdge[0][leftIndex].length + externalEdge[1][rightIndex].length + offset1 + offset2 < (*allLibraryMT)[libraryID][0].getAverageInsSize() - tolerence)
                                continue;
                            offset1 += centerNode.length/2 + externalEdge[0][leftIndex].length;
                            offset2 += centerNode.length/2 + externalEdge[1][rightIndex].length;
                            sumLinkForHaplotype[leftIndex == rightIndex] += 1;
                            MPScores[2*leftIndex + rightIndex] += 1;
                        }
                    }
                }
            }
        }
    }

}
//

//added by ouchi
void PairedDBG::getHeteroHiCLinkScore(long centerNodeIndex, const std::array<std::array<NodeIDWithGap ,2>, 2> &externalNodeID, std::array<double, 2> &sumLinkForHaplotype, std::array<double, 4> &HiCScores)
{
    GraphNode &centerNode = this->node[centerNodeIndex];
    std::array<std::array<GraphNode, 2>, 2> externalNode;
    externalNode[0][0] = this->node[id2Index(externalNodeID[0][0].id)];
    externalNode[0][1] = this->node[id2Index(externalNodeID[0][1].id)];
    externalNode[1][0] = this->node[id2Index(externalNodeID[1][0].id)];
    externalNode[1][1] = this->node[id2Index(externalNodeID[1][1].id)];

    std::array<std::array<GraphEdge, 2>, 2> externalEdge;
    externalEdge[0][0].direction = -1;
    externalEdge[0][0].end = externalNodeID[0][0].id;
    externalEdge[0][0].length = externalNodeID[0][0].gap;
    externalEdge[0][1].direction = -1;
    externalEdge[0][1].end = externalNodeID[0][1].id;
    externalEdge[0][1].length = externalNodeID[0][1].gap;
    externalEdge[1][0].direction = 1;
    externalEdge[1][0].end = externalNodeID[1][0].id;
    externalEdge[1][0].length = externalNodeID[1][0].gap;
    externalEdge[1][1].direction = 1;
    externalEdge[1][1].end = externalNodeID[1][1].id;
    externalEdge[1][1].length = externalNodeID[1][1].gap;

    std::array<std::array<vector<char>, 2>, 2> conflictpoints;
    getConflictRegionBetweenNodePair(externalEdge[0][0], externalEdge[0][1], this->node[id2Index(externalNodeID[0][0].id)], this->node[id2Index(externalNodeID[0][1].id)], conflictpoints[0][0], conflictpoints[0][1]);
    getConflictRegionBetweenNodePair(externalEdge[1][0], externalEdge[1][1], this->node[id2Index(externalNodeID[1][0].id)], this->node[id2Index(externalNodeID[1][1].id)], conflictpoints[1][0], conflictpoints[1][1]);

    std::array<long, 2> maxlengths;
    maxlengths[0] = centerNode.length / 2 + std::min(externalNodeID[0][0].gap + externalNode[0][0].length, externalNodeID[0][1].gap + externalNode[0][1].length);
    maxlengths[1] = centerNode.length / 2 + std::min(externalNodeID[1][0].gap + externalNode[1][0].length, externalNodeID[1][1].gap + externalNode[1][1].length);

    long offset1 = 0; long offset2 = 0;
    for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
        for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
            for (long HiCLinkID = 0; HiCLinkID < externalNode[0][leftIndex].HiCLinks.size(); ++HiCLinkID) {
                const HiCLink &hiclink = externalNode[0][leftIndex].HiCLinks[HiCLinkID];
                if (abs(hiclink.end) == abs(externalNodeID[1][rightIndex].id)) {
                    offset1 = 0; offset2 = 0;
                    offset1 = std::min(std::max(hiclink.offset1, offset1), externalNode[0][leftIndex].length-1);
                    offset2 = std::min(std::max(hiclink.offset2, offset2), externalNode[1][rightIndex].length-1);
                    if (conflictpoints[0][leftIndex][offset1] == 0 || conflictpoints[1][rightIndex][offset2] == 0)
                        continue;
                    if (externalNodeID[0][leftIndex].id > 0)
                        offset1 = externalNode[0][leftIndex].length - offset1 - 1;
                    if (externalNodeID[1][rightIndex].id < 0)
                        offset2 = externalNode[1][rightIndex].length - offset2 - 1;
                    offset1 += centerNode.length/2 + externalEdge[0][leftIndex].length;
                    offset2 += centerNode.length/2 + externalEdge[1][rightIndex].length;
                    sumLinkForHaplotype[leftIndex == rightIndex] += 1;//static_cast<double>(maxlengths[0] - offset1) / maxlengths[0] * (maxlengths[1] - offset2) / maxlengths[1];
                    HiCScores[2*leftIndex + rightIndex] += 1;//static_cast<double>(maxlengths[0] - offset1) / maxlengths[0] * (maxlengths[1] - offset2) / maxlengths[1];
//					sumLinkForHaplotype[leftIndex == rightIndex] += 1;
                }
            }

/*            for (long HiCLinkID = 0; HiCLinkID < externalNode[1][rightIndex].HiCLinks.size(); ++HiCLinkID) {
                const HiCLink &hiclink = externalNode[1][rightIndex].HiCLinks[HiCLinkID];
                if (abs(hiclink.end) == abs(externalNodeID[0][leftIndex].id)) {
                    if (conflictpoints[0][leftIndex][hiclink.offset2] == 0 || conflictpoints[1][rightIndex][hiclink.offset1] == 0)
                        continue;
                    if (externalNodeID[0][leftIndex].id > 0)
                        offset1 = externalNode[0][leftIndex].length - hiclink.offset2 - 1;
                    else
                        offset1 = hiclink.offset2;
                    if (externalNodeID[1][rightIndex].id > 0)
                        offset2 = hiclink.offset1;
                    else
                        offset2 = externalNode[1][rightIndex].length - hiclink.offset1 - 1;
                    offset1 += centerNode.length/2 + externalEdge[0][leftIndex].length;
                    offset2 += centerNode.length/2 + externalEdge[1][rightIndex].length;
                    sumLinkForHaplotype[leftIndex == rightIndex] += 1;//static_cast<double>(maxlengths[0] - offset1) / maxlengths[0] * (maxlengths[1] - offset2) / maxlengths[1];
                    HiCScores[2*leftIndex + rightIndex] += 1;//static_cast<double>(maxlengths[0] - offset1) / maxlengths[0] * (maxlengths[1] - offset2) / maxlengths[1];
					//sumLinkForHaplotype[leftIndex == rightIndex] += 1;
                }
            }
*/
        }
    }

}
//


//added by ouchi
bool PairedDBG::checkNodeConflict(const GraphEdge &edge1, const GraphEdge &edge2, const GraphNode &node1, const GraphNode &node2)
{
    if (edge1.direction * edge2.direction < 0)
        return 0;
    if (edge1.length + node1.length <= edge2.length) {
        if (searchPathBetweenNodePair(edge1.end, edge1.direction * sign(edge1.end), edge2.end, edge2.length - edge1.length - node1.length + tolerence))
            return 0;
    } else if (edge2.length + node2.length <= edge1.length) {
        if (searchPathBetweenNodePair(edge2.end, edge2.direction * sign(edge2.end), edge1.end, edge1.length - edge2.length - node2.length + tolerence))
            return 0;
    } else {
        long edgeEnd1, edgeEnd2;
        if (edge1.isForward()) {
            edgeEnd1 = edge1.end;
            edgeEnd2 = edge2.end;
        } else {
            edgeEnd1 = edge2.end;
            edgeEnd2 = edge1.end;
        }
        if ((abs(edge1.length + node1.length - edge2.length) <= tolerence + this->getScaffoldOverlap(edgeEnd1, edgeEnd2)
        || abs(edge2.length + node2.length - edge1.length) <= tolerence + this->getScaffoldOverlap(edgeEnd2, edgeEnd1))
        && abs(edge1.end) != abs(node2.oppositeBubbleNodeID)) return 0;
    }

    return 1;
}

//added by ouchi
void PairedDBG::detectRepeatSeverely(const double averageCoverage)
{
    const double coverageThreshold = averageCoverage * 1.75;

    // # pragma omp for schedule(dynamic)
    for (long nodeID = 0; nodeID < numNode; ++nodeID) {
        if (node[nodeID].numContig == 1 && calcNodeCoverage(node[nodeID]) > coverageThreshold) {
            node[nodeID].state |= SC_REP;
            continue;
        }
        for (long edgeID1 = 0; edgeID1 < node[nodeID].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < node[nodeID].numEdge; ++edgeID2) {
                const GraphEdge &edge1 = (node[nodeID].edge[edgeID1]);
                const GraphEdge &edge2 = (node[nodeID].edge[edgeID2]);

                if (edge1.direction * edge2.direction < 0) continue;
                GraphNode &node1 = (node[abs(edge1.end) - 1]);
                GraphNode &node2 = (node[abs(edge2.end) - 1]);

                if (this->checkNodeConflict(edge1, edge2, node1, node2)) {
                    node[nodeID].state |= SC_REP;
                    edgeID1 = node[nodeID].numEdge;
                    break;
                }
            }
        }
    }
}
//

//added by ouchi
long PairedDBG::joinUnambiguousUniqueNodePairGapped(const long numThread)
{

	vector<vector<GraphPathGapped> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<NodeIDWithGap> nodeBuffer;
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & (SC_DEL)) //| SC_REP))
				continue;

			for (long i = 0; i < 2; ++i) {
				char initialDirection = 2*i - 1;
//				getLinkedNode(initialNodeIndex, initialDirection, nodeBuffer);
				getLinkedUniqueNode(initialNodeIndex, initialDirection, nodeBuffer);
				if (nodeBuffer.size() != 1)
					continue;

				long altNodeID = nodeBuffer[0].id;
				if (node[abs(altNodeID)-1].state & (SC_DEL)) //| SC_REP))
					continue;

				getLinkedUniqueNode(abs(altNodeID) - 1, -(initialDirection * sign(altNodeID)), nodeBuffer);
				if (nodeBuffer.size() != 1)
					continue;

				pathBuffer[t].emplace_back(initialNodeIndex, 2);
				pathBuffer[t].back().node[0].id = (initialNodeIndex + 1) * initialDirection;
				pathBuffer[t].back().node[1].id = altNodeID * initialDirection;
				pathBuffer[t].back().node[0].gap = 0;
				pathBuffer[t].back().node[1].gap = nodeBuffer[0].gap;

			}
		}
	}

	vector<GraphPathGapped> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	long i; long j;
	for (i = 0; i < mergedPathBuffer.size(); ++i) {
		for (j = 0; j < mergedPathBuffer.size(); ++j) {
			if (j == i) continue;
			if (mergedPathBuffer[i].node.back().id == mergedPathBuffer[j].node.front().id) {
				for (long k = 1; k < mergedPathBuffer[j].node.size(); ++k) {
					mergedPathBuffer[i].node.push_back(mergedPathBuffer[j].node[k]);
				}
				mergedPathBuffer.erase(mergedPathBuffer.begin() + j);
				--i;
				if (j < i) --i;
				break;
			}
		}
	}
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), GraphPathGappedSelfIDLess());

	long numJoin = remakeGraphAccordingToGappedPath(mergedPathBuffer);

	return numJoin;
}
//

//added by ouchi
void PairedDBG::joinUnambiguousUniqueNodePairGappedIterative(const long numThread)
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
		num = joinUnambiguousUniqueNodePairGapped(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);
	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}
//added by ouchi

//added by ouchi
void PairedDBG::joinNoBranchNodePairGappedIterative(const long numThread)
{
	long total = 0;
	long num;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "joining nobranch node pair in a scaffold graph.." << endl;
	do {
		setMinLink(1);
		makeGraph(numThread);
		setMinLink(currentMinLink);
		num = joinNoBranchNodePairGapped(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);
	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}
//

//added by ouchi
void PairedDBG::joinNoBranchNodePairGappedIterativeAllLibraries(const long numThread)
{
	long total = 0;
	long num;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "joining nobranch node pair in a scaffold graph using all libraries.." << endl;
	do {
		setMinLink(1);
		makeGraphAllLibraries(numThread);
		setMinLink(currentMinLink);
		num = joinNoBranchNodePairGapped(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);
	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}
//

//added by ouchi
long PairedDBG::joinNoBranchNodePairGapped(const long numThread)
{

	const unsigned long currentMinLink = this->minLink;
	setMinLink(1);
	vector<vector<GraphPathGapped> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<NodeIDWithGap> nodeBuffer;
		GraphEdge edge1, edge2;
		for (long initialNodeIndex = t; initialNodeIndex < numNode; initialNodeIndex += numThread) {
			if (node[initialNodeIndex].state & (SC_DEL))
				continue;

			for (long i = 0; i < 2; ++i) {
				char initialDirection = 2*i - 1;
				getLinkedNode(initialNodeIndex, initialDirection, nodeBuffer);
				if (nodeBuffer.size() < 1) continue;
				std::sort(nodeBuffer.begin(), nodeBuffer.end(), [](NodeIDWithGap const& left, NodeIDWithGap const& right) {
					return left.gap < right.gap;
				});
				if (initialDirection < 0) {
					if (getNumLinkFromIDPair(nodeBuffer[0].id, initialNodeIndex + 1) < currentMinLink)
						continue;
				} else {
					if (getNumLinkFromIDPair(initialNodeIndex + 1, nodeBuffer[0].id) < currentMinLink)
						continue;
				}
				long j = 1;
				for (; j < nodeBuffer.size(); ++j) {
					edge1.direction = initialDirection;
					edge1.end = nodeBuffer[0].id;
					edge1.length = nodeBuffer[0].gap;
					edge2.direction = initialDirection;
					edge2.end = nodeBuffer[j].id;
					edge2.length = nodeBuffer[j].gap;
					if (checkNodeConflict(edge1, edge2, node[id2Index(nodeBuffer[0].id)], node[id2Index(nodeBuffer[j].id)]))
						break;
				}
				if (j < nodeBuffer.size())
					continue;

				long altNodeID = nodeBuffer[0].id;
				if (node[abs(altNodeID)-1].state & (SC_DEL))
					continue;

				getLinkedNode(abs(altNodeID) - 1, -(initialDirection * sign(altNodeID)), nodeBuffer);
				std::sort(nodeBuffer.begin(), nodeBuffer.end(), [](NodeIDWithGap const& left, NodeIDWithGap const& right) {
					return left.gap < right.gap;
				});
				j = 1;
				for (; j < nodeBuffer.size(); ++j) {
					edge1.direction = -(initialDirection * sign(altNodeID));
					edge1.end = nodeBuffer[0].id;
					edge1.length = nodeBuffer[0].gap;
					edge2.direction = -(initialDirection * sign(altNodeID));
					edge2.end = nodeBuffer[j].id;
					edge2.length = nodeBuffer[j].gap;
					if (checkNodeConflict(edge1, edge2, node[id2Index(nodeBuffer[0].id)], node[id2Index(nodeBuffer[j].id)]))
						break;
				}
				if (j < nodeBuffer.size())
					continue;

				pathBuffer[t].emplace_back(initialNodeIndex, 2);
				pathBuffer[t].back().node[0].id = (initialNodeIndex + 1) * initialDirection;
				pathBuffer[t].back().node[1].id = altNodeID * initialDirection;
				pathBuffer[t].back().node[0].gap = 0;
				pathBuffer[t].back().node[1].gap = nodeBuffer[0].gap;

			}
		}
	}

	setMinLink(currentMinLink);

	vector<GraphPathGapped> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	long i; long j;
	for (i = 0; i < mergedPathBuffer.size(); ++i) {
		for (j = 0; j < mergedPathBuffer.size(); ++j) {
			if (j == i) continue;
			if (mergedPathBuffer[i].node.back().id == mergedPathBuffer[j].node.front().id) {
				for (long k = 1; k < mergedPathBuffer[j].node.size(); ++k) {
					mergedPathBuffer[i].node.push_back(mergedPathBuffer[j].node[k]);
				}
				mergedPathBuffer.erase(mergedPathBuffer.begin() + j);
				--i;
				if (j < i) --i;
				break;
			}
		}
	}

	long numJoin = remakeGraphAccordingToGappedPath(mergedPathBuffer);

	return numJoin;
}
//

//added by ouchi
long PairedDBG::joinUnambiguousNodePairGapped(const long numThread)
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

				//added by ouchi for test
				//if (abs(calcNodeCoverage(this->node[initialNodeIndex]) - calcNodeCoverage(this->node[id2Index(altNodeID)])) > this->heteroCoverage)
				//	continue;
				//

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
//

/*
//added by ouchi
void PairedDBG::solveRepeatStructure(const long numThread)
{
    long numNewContig = 0;
    long numNewNode = 0;
    long newContigPoolSize = 0;

    cerr << "solve RepeatStructure" << endl;

    for (long i = 0; i < numNode; ++i) {
        if (node[i].state & (SC_INC | SC_REP | SC_DEL)) continue;

    }

}
//

//added by ouchi
void PairedDBG::searchUniqueNodePath(const long nodeID, const char targetDirection, vector<GraphPathGapped>& pathBuffer)
{
    const GraphNode &sourceNode = node[id2Index(nodeID)];
    for (long edgeIndex = 0; edgeIndex < sourceNode.numEdge; ++edgeIndex) {
        const GraphEdge &edge = sourceNode.edge[edgeIndex];
        if (edge.direction == targetDirection && edge.numLink >= minLink) {
            const GraphNode &tmpNode = node[id2Index(edge.end)];
            if (tmpNode.state & SC_DEL) {
                continue;
            }
            if (tmpNode.state & SC_REP) {
                pathBuffer.back().node.emplace_back(edge.end, edge.length);
                searchUniqueNodePath(edge.end,  targetDirection * sign(edge.end), pathBuffer);
            } else {
                if (pathBuffer.back().node.size() > 0) {
                    pathBuffer.push_back(pathBuffer.back());
                    pathBuffer[pathBuffer.size()-2].node.emplace_back(edge.end, edge.length);
                } else {
                    pathBuffer.back().node.emplace_back(edge.end, edge.length);
                    pathBuffer.push_back();
                }
            }
        }
    }
    pathBuffer.pop();

    return;
}
// */

/*void PairedDBG::divideErroneousNodeBaseLevel(const long numThread)
{
	long currectLibraryIndex = this->targetLibraryIndex;

	cerr << "dividing erroneous scaffolds based on base-level coverages ..." << endl;

    vector<vector<unsigned> > physicalCoverage(numNode);
    vector<vector<unsigned> > diffcoverage(numNode);
    for (unsigned i = 0; i < numNode; ++i) {
        physicalCoverage[i].resize(node[i].length, 0);
        diffCoverage[i].resize(node[i].length, 0);
    }

    omp_set_num_threads(numThread);

    if (allLibraryMT != NULL) {
        for (libraryIndex = 0; libraryIndex < (*allLibraryMT).size(); ++libraryIndex) {
            this->setTargetLibraryIndex(libraryIndex);
            this->setTolerence(minTolerenceFactor * (*allLibraryMT)[libraryIndex][0].getSDInsize());

            this->calculatePhysicalCoverage(physicalCoverage, (*allLibraryMT)[libraryIndex][0].getAverageInsSize(), numThread);
            this->compensatePhysicalCoverageBasedOnGapRate(physicalCoverage, 2 * (*allLibraryMT)[libraryIndex][0].getAverageInsSize(), numThread);
            this->calculateDiffCoverage(diffCoverage, 0, numThread);
            this->divideNodeBasedOnCoverage(physicalCoverage, diffCoverage, numThread);
        }
    }

	this->setTargetLibraryIndex(currentLibraryIndex);
}

void PairedDBG::divdiNodeBasedOnCoverage(const vector<vector<unsigned> >& pyhsicalCoverage, const vector<vector<unsigned> >& diffCoverage, const long numThread)
{
	const long EDGE_LENGTH = (*allLibaryMT)[targetLibraryIndex][0].getAverageInsSize();
	const long Min_COVERAGE_RATE = 0.5;

	long numDivided = 0;

	if (this->contigPreviousParentNodeID.size() != this->contig.size())
		this->contigPreviousParentNodeID.resize(this->contig.size(), 0);

    omp_set_num_threads(numThread);

    #pragma omp parallel for schedule(static, 1) reduction(+: numDivided)
	for (unsigned threadIndex = 0; threadIndex < numThread; ++threadIndex) {
		vector<char> baseBreakpoint;
		vector<long> diffPeak;
		for (long i = threadIndex; i < numNode; i+= numThread) {
			if (physicalCoverage[i].size() <= 2+EDGE_LENGTH)
				continue;

			vector<unsigned> buffer(physicalCoverage.begin() + edgeLength, physicalCoverage.end() - edgeLength);
			std::partial_sort(buffer.begin(), buffer.begin() + buffer.size() / 2 + 1, buffer.end());
			const unsigned medianPhysicalCoverage = buffer[buffer.size() / 2];

			detectPeak(diffCoverage,0,0,0,diffPeak); //after check
			for (unsigned int j = 0; j < diffPeak.size(); ++j) {
				if (physicalCoverage[diffPeak[j]] < MIN_COVERAGE_RATE * medianPhysicalCoverage) {
					baseBreakpoint[diffPeak[j]] = 1;
				} else {

				}
			}
		}
	}

}*/

//added by ouchi
void PairedDBG::initConsensusNode(const long numThread)
{
    const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;

    numConsensusNode = numNode;
    consensusnode.resize(numConsensusNode);
    nodePositionInConsensus.resize(numNode);

    for (int i = 0; i < numConsensusNode; ++i) {
        consensusnode[i].length = this->node[i].length;
        consensusnode[i].numNode[0] = 1;
        consensusnode[i].node[0].resize(consensusnode[i].numNode[0]);
        consensusnode[i].node[0][0].end = this->node[i].length;
        consensusnode[i].node[0][0].id = i + 1;
        if (calcNodeCoverage(this->node[i]) <= COVERAGE_THRESHOLD) {
            consensusnode[i].numNode[1] = 0;
            consensusnode[i].node[1].resize(consensusnode[i].numNode[1]);
            nodePositionInConsensus[i].resize(1);
            nodePositionInConsensus[i][0].first = 0;
            nodePositionInConsensus[i][0].second.id = i + 1;
            nodePositionInConsensus[i][0].second.offset = 0;
        } else {
            consensusnode[i].numNode[1] = 1;
            consensusnode[i].node[1].resize(consensusnode[i].numNode[1]);
            consensusnode[i].node[1][0].end = this->node[i].length;
            consensusnode[i].node[1][0].id = i + 1;
            nodePositionInConsensus[i].resize(2);
            nodePositionInConsensus[i][0].first = 0;
            nodePositionInConsensus[i][0].second.id = i + 1;
            nodePositionInConsensus[i][0].second.offset = 0;
            nodePositionInConsensus[i][1].first = 1;
            nodePositionInConsensus[i][1].second.id = i + 1;
            nodePositionInConsensus[i][1].second.offset = 0;
        }
    }

}
//

//added by ouchi
void PairedDBG::remakeConsensus(const long numNewConsensusNode, const long newNodePoolSize, FILE *scaffoldFP)
{
    consensusnode.clear();
    numConsensusNode = numNewConsensusNode;
    consensusnode.resize(numConsensusNode);

    for (long i = 0; i < numNode; ++i) {
        nodePositionInConsensus[i].clear();
    }

    rewind(scaffoldFP);
    for (long i = 0; i < numConsensusNode; ++i) {
        fread(&(consensusnode[i].numNode[0]), sizeof(long), 1, scaffoldFP);
        consensusnode[i].node[0].resize(consensusnode[i].numNode[0]);
        for (long j = 0; j< consensusnode[i].numNode[0]; ++j)
            fread(&(consensusnode[i].node[0][j]), sizeof(ScaffoldPart), 1, scaffoldFP);
        fread(&(consensusnode[i].numNode[1]), sizeof(long), 1, scaffoldFP);
        consensusnode[i].node[1].resize(consensusnode[i].numNode[1]);
        for (long j = 0; j< consensusnode[i].numNode[1]; ++j)
            fread(&(consensusnode[i].node[1][j]), sizeof(ScaffoldPart), 1, scaffoldFP);

        long start = consensusnode[i].node[0].front().start;
        long end = consensusnode[i].node[0].back().end;
        for (long j = 0; j < consensusnode[i].numNode[0]; ++j) {
            if (start > consensusnode[i].node[0][j].start)
                start = consensusnode[i].node[0][j].start;
            if (end < consensusnode[i].node[0][j].end)
                end = consensusnode[i].node[0][j].end;
        }
        long min_start = start;
        long max_end = end;
        if (consensusnode[i].numNode[1] > 0) {
            long start2 = consensusnode[i].node[1].front().start;
            long end2 = consensusnode[i].node[1].back().end;
            for (long j = 0; j < consensusnode[i].numNode[1]; ++j) {
                if (start2 > consensusnode[i].node[1][j].start)
                    start2 = consensusnode[i].node[1][j].start;
                if (end2 < consensusnode[i].node[1][j].end)
                    end2 = consensusnode[i].node[1][j].end;
            }
            min_start = std::min(start, start2);
            max_end = std::max(end, end2);
        }
        if (min_start != 0) {
            for (long j = 0; j < consensusnode[i].numNode[0]; ++j) {
                consensusnode[i].node[0][j].start -= min_start;
                consensusnode[i].node[0][j].end -= min_start;
            }
            for (long j = 0; j < consensusnode[i].numNode[1]; ++j) {
                consensusnode[i].node[1][j].start -= min_start;
                consensusnode[i].node[1][j].end -= min_start;
            }
        }

        consensusnode[i].length = max_end - min_start;

        consensusnode[i].state = 0;

        for (long j = 0; j < consensusnode[i].numNode[0]; ++j) {
            long tmp = abs(consensusnode[i].node[0][j].id) - 1;
            if (consensusnode[i].node[0][j].id > 0)
                nodePositionInConsensus[tmp].emplace_back(0, platanus::Position(j, i + 1));
            else
                nodePositionInConsensus[tmp].emplace_back(0, platanus::Position(j, -(i + 1)));
        }
        for (long j = 0; j < consensusnode[i].numNode[1]; ++j) {
            long tmp = abs(consensusnode[i].node[1][j].id) - 1;
            if (consensusnode[i].node[1][j].id > 0)
                nodePositionInConsensus[tmp].emplace_back(1, platanus::Position(j, i + 1));
            else
                nodePositionInConsensus[tmp].emplace_back(1, platanus::Position(j, -(i + 1)));
        }

    }

}
//

//added by ouchi
void PairedDBG::writeSingletonConsensusNode(long &numNewConsensusNode, long &newNodePoolSize, FILE *storeFP)
{
    for (long i = 0; i < numConsensusNode; ++i) {
        if (!(consensusnode[i].state & SC_INC) && !(consensusnode[i].state & SC_DEL)) {
            long numNode1 = consensusnode[i].numNode[0];
            fwrite(&(consensusnode[i].numNode[0]), sizeof(long), 1, storeFP);
            for (long j = 0; j < numNode1; ++j)
                fwrite(&(consensusnode[i].node[0][j]), sizeof(ScaffoldPart), 1, storeFP);
            newNodePoolSize += numNode1;
            long numNode2 = consensusnode[i].numNode[1];
            fwrite(&(consensusnode[i].numNode[1]), sizeof(long), 1, storeFP);
            for (long j = 0; j < numNode2; ++j)
                fwrite(&(consensusnode[i].node[1][j]), sizeof(ScaffoldPart), 1, storeFP);
            newNodePoolSize += numNode2;
            ++numNewConsensusNode;
        }
    }
}
//

//added by ouchi
long PairedDBG::writeAndMarkGappedConsensusNodes(const vector<ConsensusIDWithGap> &gappedPath, FILE *storeFP)
{
	long sumNumNode = 0;
	for (char h = 0; h < 2; ++h) {
		long numNode = 0;
		for (auto it = gappedPath.begin(); it != gappedPath.end(); ++it)
			numNode += consensusnode[abs(it->id)-1].numNode[(it->h + h) % 2];
		fwrite(&numNode, sizeof(long), 1, storeFP);
		sumNumNode += numNode;

		long start = 0;
		for (auto it = gappedPath.begin(); it != gappedPath.end(); ++it) {
			start += it->gap;
			long index;
			if (it->id > 0) {
				index = it->id -1;
				consensusnode[index].state |= SC_INC;
				for (long i = 0; i < consensusnode[index].numNode[(it->h + h) % 2]; ++i) {
					ScaffoldPart scaffoldPartForWriting(consensusnode[index].node[(it->h + h) % 2][i].id, start + consensusnode[index].node[(it->h + h) % 2][i].start, start + consensusnode[index].node[(it->h + h) % 2][i].end);
					fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
				}
			} else {
				index = -(it->id) -1;
				consensusnode[index].state |= SC_INC;
				for (long i = consensusnode[index].numNode[(it->h + h) % 2] - 1; i >= 0; --i) {
					ScaffoldPart scaffoldPartForWriting(-(consensusnode[index].node[(it->h + h) % 2][i].id), start + consensusnode[index].length - consensusnode[index].node[(it->h + h) % 2][i].end, start + consensusnode[index].length - consensusnode[index].node[(it->h + h) % 2][i].start);
					fwrite(&scaffoldPartForWriting, sizeof(ScaffoldPart), 1, storeFP);
				}
			}
			start += consensusnode[index].length;
		}
	}
	return sumNumNode;
}
//

//added by ouchi
long PairedDBG::remakeConsensusGraphAccordingToGappedPath(vector<ConsensusPathGapped> &pathBuffer)
{
	long numNewConsensusNode = 0;
	long newNodePoolSize = 0;
	long numJoin = 0;

	FILE *newContigFP = platanus::makeTemporaryFile();

	for (auto pathItr = pathBuffer.begin(); pathItr != pathBuffer.end(); ++pathItr) {
		auto consensusnodeItr = pathItr->consensusnode.begin();
		for (; consensusnodeItr != pathItr->consensusnode.end(); ++consensusnodeItr) {
			if (consensusnode[abs(consensusnodeItr->id)-1].state & SC_INC)
				break;
		}
		if (consensusnodeItr != pathItr->consensusnode.end())
			continue;

		numJoin += pathItr->consensusnode.size() - 1;
		newNodePoolSize += writeAndMarkGappedConsensusNodes(pathItr->consensusnode, newContigFP);
		++numNewConsensusNode;
	}

	writeSingletonConsensusNode(numNewConsensusNode, newNodePoolSize, newContigFP);
	remakeConsensus(numNewConsensusNode, newNodePoolSize, newContigFP);
	fclose(newContigFP);

	return numJoin;
}
//

//added by ouchi
bool PairedDBG::judgeConflictConsensusNode(const long offset1, const long offset2, const long consensusID1, const long consensusID2, const char h1, const char h2, long tol)
{
    const ConsensusNode &consensusnode1 = this->consensusnode[id2Index(consensusID1)];
    const ConsensusNode &consensusnode2 = this->consensusnode[id2Index(consensusID2)];
    char d1 = sign(consensusID1);
    char d2 = sign(consensusID2);
    for (long nodeIndex1 = 0; nodeIndex1 < consensusnode1.numNode[h1]; ++nodeIndex1) {
        long start, end;
        if (d1 > 0) {
            start = consensusnode1.node[h1][nodeIndex1].start + offset1;
            end = consensusnode1.node[h1][nodeIndex1].end + offset1;
        } else {
            start = (consensusnode1.length - consensusnode1.node[h1][nodeIndex1].end) + offset1;
            end = (consensusnode1.length - consensusnode1.node[h1][nodeIndex1].start) + offset1;
        }
        for (long nodeIndex2 = 0; nodeIndex2 < consensusnode2.numNode[h2]; ++nodeIndex2) {
            long start2, end2;
            if (d2 > 0) {
                start2 = consensusnode2.node[h2][nodeIndex2].start + offset2;
                end2 = consensusnode2.node[h2][nodeIndex2].end + offset2;
            } else {
                start2 = (consensusnode2.length - consensusnode2.node[h2][nodeIndex2].end) + offset2;
                end2 = (consensusnode2.length - consensusnode2.node[h2][nodeIndex2].start) + offset2;
            }
            if (start < start2) {
                if (end < start2) {
                    continue;
                } else if (end2 < end) {
                    return true;
                } else if (start2 < end - (tol + this->getScaffoldOverlap(d1*consensusnode1.node[h1][nodeIndex1].id, d2*consensusnode2.node[h2][nodeIndex2].id))) {
                    return true;
                }
            } else {
                if (end2 < start) {
                    continue;
                } else if (end < end2) {
                    return true;
                } else if (start < end2 - (tol + this->getScaffoldOverlap(d2*consensusnode2.node[h2][nodeIndex2].id, d1*consensusnode1.node[h1][nodeIndex1].id))) {
                    return true;
                }

            }
        }
    }
    for (long nodeIndex1 = 0; nodeIndex1 < consensusnode1.numNode[1-h1]; ++nodeIndex1) {
        long start, end;
        if (d1 > 0) {
            start = consensusnode1.node[1-h1][nodeIndex1].start + offset1;
            end = consensusnode1.node[1-h1][nodeIndex1].end + offset1;
        } else {
            start = (consensusnode1.length - consensusnode1.node[1-h1][nodeIndex1].end) + offset1;
            end = (consensusnode1.length - consensusnode1.node[1-h1][nodeIndex1].start) + offset1;
        }
        for (long nodeIndex2 = 0; nodeIndex2 < consensusnode2.numNode[1-h2]; ++nodeIndex2) {
            long start2, end2;
            if (d2 > 0) {
                start2 = consensusnode2.node[1-h2][nodeIndex2].start + offset2;
                end2 = consensusnode2.node[1-h2][nodeIndex2].end + offset2;
            } else {
                start2 = (consensusnode2.length - consensusnode2.node[1-h2][nodeIndex2].end) + offset2;
                end2 = (consensusnode2.length - consensusnode2.node[1-h2][nodeIndex2].start) + offset2;
            }
            if (start < start2) {
                if (end < start2) {
                    continue;
                } else if (end2 < end) {
                    return true;
                } else if (start2 < end - (tol + this->getScaffoldOverlap(d1*consensusnode1.node[1-h1][nodeIndex1].id, d2*consensusnode2.node[1-h2][nodeIndex2].id))) {
                    return true;
                }
            } else {
                if (end2 < start) {
                    continue;
                } else if (end < end2) {
                    return true;
                } else if (start < end2 - (tol + this->getScaffoldOverlap(d2*consensusnode2.node[1-h2][nodeIndex2].id, d1*consensusnode1.node[1-h1][nodeIndex1].id))) {
                    return true;
                }

            }
        }
    }

    return false;
}
//

//added by ouchi
void PairedDBG::setOppositeBubbleConsensusNode(const long consensusnodeIndex, std::array<long, 4>  &oppositeConsensusInfo)
{
    const double COVERAGE_THRESHOLD = HETERO_COVERAGE_THRESHOLD_FACTOR * this->heteroCoverage;
    oppositeConsensusInfo = {0,0,0,0};

    const ConsensusNode &tmpConsensus = this->consensusnode[consensusnodeIndex];
    std::unordered_map<std::pair<long, char>, std::array<long, 4>, platanus::PairHash, platanus::PairEqual> countMap;
    for (char h = 0; h < 2; ++h) {
        for (long nodeIndex = 0; nodeIndex < tmpConsensus.numNode[h]; ++nodeIndex) {
            const GraphNode &tmpNode = this->node[id2Index(this->consensusnode[consensusnodeIndex].node[h][nodeIndex].id)];
            const char nodeDirection = sign(consensusnode[consensusnodeIndex].node[h][nodeIndex].id);
            if (calcNodeCoverage(tmpNode) > COVERAGE_THRESHOLD)
                continue;
            for (auto partItr = tmpNode.contig.begin(); partItr != tmpNode.contig.end(); ++partItr) {
                const char contigDirection = nodeDirection * sign(partItr->id);
                if (this->contigPositionInScaffold[id2Index(partItr->id)].size() != 1)
                    continue;
                std::array<char, 2> oppositeDirection;
                std::array<platanus::Position, 2> tmpPosition;
                for (long i = 0; i < 2; ++i) {
                    long j = partItr->id > 0 ? i : (i + 1)%2;
                    if (this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i] != 0) {
                         if (this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])].size() == 1) {
                             tmpPosition[j] = this->contigPositionInScaffold[id2Index(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i])][0];
                         }
                         oppositeDirection[j] = contigDirection * sign(this->contigBubbleInfo[id2Index(partItr->id)].oppositeContigID[i]);
                    }
                }
                if (tmpPosition[0].id != 0 && tmpPosition[0] == tmpPosition[1] && oppositeDirection[0] == oppositeDirection[1]) {
                    if (calcNodeCoverage(this->node[id2Index(tmpPosition[0].id)]) > COVERAGE_THRESHOLD)
                        continue;
                    const char oppositeNodeDirection = oppositeDirection[0] * sign(tmpPosition[0].id);
                    if (nodePositionInConsensus[id2Index(tmpPosition[0].id)].size() == 1) {
                        long oppositeConsensusID = nodePositionInConsensus[id2Index(tmpPosition[0].id)][0].second.id;
                        const char oppositeConsensusDirection = oppositeNodeDirection * sign(oppositeConsensusID);
                        char oh = nodePositionInConsensus[id2Index(tmpPosition[0].id)][0].first;

                        if (abs(oppositeConsensusID) == consensusnodeIndex + 1)
                            continue;

                        long startContig = this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start;
                        if (nodeDirection > 0)
                            startContig += partItr->start;
                        else
                            startContig += tmpNode.length - partItr->end;
                        long startOppositeContig;
                        if (oppositeConsensusDirection > 0)
                            startOppositeContig = this->consensusnode[id2Index(oppositeConsensusID)].node[oh][nodePositionInConsensus[id2Index(tmpPosition[0].id)][0].second.offset].start;
                        else
                            startOppositeContig = this->consensusnode[id2Index(oppositeConsensusID)].length - this->consensusnode[id2Index(oppositeConsensusID)].node[oh][nodePositionInConsensus[id2Index(tmpPosition[0].id)][0].second.offset].end;
                        if (oppositeNodeDirection > 0)
                            startOppositeContig += this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].start;
                        else
                            startOppositeContig += this->node[id2Index(tmpPosition[0].id)].length - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].end;

                        auto countItr = countMap.find(std::make_pair(abs(oppositeConsensusID) * oppositeConsensusDirection, oh));
                        if (countItr != countMap.end()) {
                            countItr->second[2*h] += 1;
                            countItr->second[2*h+1] += startContig - startOppositeContig;
                        } else {
                            countMap[std::make_pair(abs(oppositeConsensusID) * oppositeConsensusDirection, oh)][2*h] = 1;
                            countMap[std::make_pair(abs(oppositeConsensusID) * oppositeConsensusDirection, oh)][2*h+1] = startContig - startOppositeContig;
                        }
                    }


                    /*auto countItr = countMap.find(tmpPosition[0].id);
                    if (countItr != countMap.end()) {
                        countItr->second[2*h] += 1;
                        if (this->consensusnode[consensusnodeIndex].node[h][nodeIndex].id > 0) {
                            if (tmpPosition[0].id > 0) {
                                countItr->second[2*h+1] += this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + partItr->start - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].start; //
                            } else {
                                countItr->second[2*h+1] += this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + partItr->start - (this->node[id2Index(tmpPosition[0].id)].length - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].end); //
                            }
                        } else {
                            if (tmpPosition[0].id > 0) {
                                countItr->second[2*h+1] += this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + (tmpNode.length - partItr->end) - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].start; //
                            } else {
                                countItr->second[2*h+1] += this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + (tmpNode.length - partItr->end) - (this->node[id2Index(tmpPosition[0].id)].length - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].end); //
                            }
                        }
                    } else {
                        if (calcNodeCoverage(this->node[id2Index(tmpPosition[0].id)]) <= COVERAGE_THRESHOLD) {
							std::cerr << "h:" << h << "\tnodeID:" << tmpPosition[0].id << std::endl;
                            countMap[tmpPosition[0].id][2*h] = 1;
                            if (this->consensusnode[consensusnodeIndex].node[h][nodeIndex].id > 0) {
                                if (tmpPosition[0].id > 0) {
                                    countMap[tmpPosition[0].id][2*h+1] = this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + partItr->start - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].start; //
                                } else {
                                    countMap[tmpPosition[0].id][2*h+1] = this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + partItr->start - (this->node[id2Index(tmpPosition[0].id)].length - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].end); //
                                }
                            } else {
                                if (tmpPosition[0].id > 0) {
                                    countMap[tmpPosition[0].id][2*h+1] = this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + (tmpNode.length - partItr->end) - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].start; //
                                } else {
                                    countMap[tmpPosition[0].id][2*h+1] = this->consensusnode[consensusnodeIndex].node[h][nodeIndex].start + (tmpNode.length - partItr->end) - (this->node[id2Index(tmpPosition[0].id)].length - this->node[id2Index(tmpPosition[0].id)].contig[tmpPosition[0].offset].end); //
                                }
                            }
                        }
                    }*/
                }
            }
        }
    }

    while(!countMap.empty()) {
        long maxID = 0;
        long maxCount = 0;
        long maxOffset = 0;
        char h = -1;
        char oh = -1;
        for (auto countItr = countMap.begin(); countItr != countMap.end(); ++countItr) {
            if (maxCount < countItr->second[0]) {
                maxID = countItr->first.first;
                maxCount = countItr->second[0];
                maxOffset = countItr->second[1];
                oh = countItr->first.second;
                h = 0;
            }
            if (maxCount < countItr->second[2]) {
                maxID = countItr->first.first;
                maxCount = countItr->second[2];
                maxOffset = countItr->second[3];
                oh = countItr->first.second;
                h = 1;
            }
        }
        if (maxID == 0)
            break;
        if (judgeConflictConsensusNode(0, maxOffset / maxCount, consensusnodeIndex+1, maxID, 1-h, oh, 0)) {
            countMap[std::make_pair(maxID, oh)][2*h] = 0;
            continue;
        } else {
            oppositeConsensusInfo[0] = h;
            oppositeConsensusInfo[1] = maxID;
            oppositeConsensusInfo[2] = oh;
            oppositeConsensusInfo[3] = maxOffset / maxCount;
            break;
        }
        /*if (this->nodePositionInConsensus[id2Index(maxID)].size() != 1) {
            countMap[maxID][2*h] = 0;
            continue;
        }
        const ConsensusNode &oppositeConsensus = this->consensusnode[id2Index(this->nodePositionInConsensus[id2Index(maxID)][0].second.id)];
        long startConsensusNode = countMap[maxID][2*h+1] / maxCount;
        if (this->nodePositionInConsensus[id2Index(maxID)][0].second.id > 0)
            startConsensusNode -= oppositeConsensus.node[this->nodePositionInConsensus[id2Index(maxID)][0].first][this->nodePositionInConsensus[id2Index(maxID)][0].second.offset].start;
        else
            startConsensusNode -= oppositeConsensus.length - oppositeConsensus.node[this->nodePositionInConsensus[id2Index(maxID)][0].first][this->nodePositionInConsensus[id2Index(maxID)][0].second.offset].end;
        if (judgeConflictConsensusNode(0, startConsensusNode, consensusnodeIndex+1, this->nodePositionInConsensus[id2Index(maxID)][0].second.id, 1-h, this->nodePositionInConsensus[id2Index(maxID)][0].first)) {
            countMap[maxID][2*h] = 0;
            continue;
        } else {
            oppositeConsensusInfo[0] = h;
            oppositeConsensusInfo[1] = this->nodePositionInConsensus[id2Index(maxID)][0].second.id;
            oppositeConsensusInfo[2] = this->nodePositionInConsensus[id2Index(maxID)][0].first;
            oppositeConsensusInfo[3] = startConsensusNode;
            break;
        }*/
    }
}
//

//added by ouchi
void PairedDBG::mergeHeteroNode(const long numThread)
{
    vector<std::array<long, 4> > oppositeConsensusInfo;
    oppositeConsensusInfo.resize(numConsensusNode);

    for (long consensusnodeIndex = 0; consensusnodeIndex < static_cast<long>(this->consensusnode.size()); ++consensusnodeIndex) {
        setOppositeBubbleConsensusNode(consensusnodeIndex, oppositeConsensusInfo[consensusnodeIndex]);
    }

    for (long consensusnodeIndex = 0; consensusnodeIndex < numConsensusNode; ++consensusnodeIndex) {
        if (this->consensusnode[consensusnodeIndex].state & SC_INC)
            continue;
        long oppositeConsensusID = oppositeConsensusInfo[consensusnodeIndex][1];
        if (oppositeConsensusID == 0)
            continue;
        if (this->consensusnode[id2Index(oppositeConsensusID)].state & SC_INC) {
            setOppositeBubbleConsensusNode(consensusnodeIndex, oppositeConsensusInfo[consensusnodeIndex]);
            oppositeConsensusID = oppositeConsensusInfo[consensusnodeIndex][1];
            if (oppositeConsensusID == 0)
                continue;
        }
        if (abs(oppositeConsensusInfo[id2Index(oppositeConsensusID)][1]) != consensusnodeIndex + 1)
            continue;

        ++numConsensusNode;
        consensusnode.resize(numConsensusNode);

        char h1 = oppositeConsensusInfo[consensusnodeIndex][0];
        char h2 = oppositeConsensusInfo[consensusnodeIndex][2];
        long startConsensusNode = oppositeConsensusInfo[consensusnodeIndex][3];
        ConsensusNode& tmpConsensus1 = this->consensusnode[consensusnodeIndex];
        ConsensusNode& tmpConsensus2 = this->consensusnode[id2Index(oppositeConsensusID)];
        tmpConsensus1.state |= SC_INC;
        tmpConsensus2.state |= SC_INC;
        for (char h = 0; h < 2; ++h) {
            for (long nodeIndex = 0; nodeIndex < tmpConsensus1.numNode[h]; ++nodeIndex) {
                long tmp = id2Index(tmpConsensus1.node[h][nodeIndex].id);
                for (long i = 0; i < nodePositionInConsensus[tmp].size(); ++i) {
                    if (nodePositionInConsensus[tmp][i].first == h
                     && nodePositionInConsensus[tmp][i].second.offset == nodeIndex
                     && abs(nodePositionInConsensus[tmp][i].second.id) == consensusnodeIndex + 1)
                    nodePositionInConsensus[tmp].erase(nodePositionInConsensus[tmp].begin()+i);
                    break;
                }
            }
            for (long nodeIndex = 0; nodeIndex < tmpConsensus2.numNode[h]; ++nodeIndex) {
                long tmp = id2Index(tmpConsensus2.node[h][nodeIndex].id);
                for (long i = 0; i < nodePositionInConsensus[tmp].size(); ++i) {
                    if (nodePositionInConsensus[tmp][i].first == h
                     && nodePositionInConsensus[tmp][i].second.offset == nodeIndex
                     && abs(nodePositionInConsensus[tmp][i].second.id) == abs(oppositeConsensusID))
                    nodePositionInConsensus[tmp].erase(nodePositionInConsensus[tmp].begin()+i);
                    break;
                }
            }
        }

        ConsensusNode& newConsensus = consensusnode[numConsensusNode-1];
        newConsensus.numNode[0] = tmpConsensus1.numNode[h1] + tmpConsensus2.numNode[1-h2];
        newConsensus.node[0].resize(newConsensus.numNode[0]);
        for (long nodeIndex = 0; nodeIndex < tmpConsensus1.numNode[h1]; ++nodeIndex) {
            newConsensus.node[0][nodeIndex] = tmpConsensus1.node[h1][nodeIndex];
        }
        for (long nodeIndex = 0; nodeIndex < tmpConsensus2.numNode[1-h2]; ++nodeIndex) {
            newConsensus.node[0][tmpConsensus1.numNode[h1] + nodeIndex].id = tmpConsensus2.node[1-h2][nodeIndex].id * sign(oppositeConsensusID);
            if (oppositeConsensusID > 0) {
                newConsensus.node[0][tmpConsensus1.numNode[h1] + nodeIndex].start = tmpConsensus2.node[1-h2][nodeIndex].start + startConsensusNode;
                newConsensus.node[0][tmpConsensus1.numNode[h1] + nodeIndex].end = tmpConsensus2.node[1-h2][nodeIndex].end + startConsensusNode;
            } else {
                newConsensus.node[0][tmpConsensus1.numNode[h1] + nodeIndex].start = (tmpConsensus2.length - tmpConsensus2.node[1-h2][nodeIndex].end) + startConsensusNode;
                newConsensus.node[0][tmpConsensus1.numNode[h1] + nodeIndex].end = (tmpConsensus2.length - tmpConsensus2.node[1-h2][nodeIndex].start) + startConsensusNode;
            }
        }
        std::sort(newConsensus.node[0].begin(), newConsensus.node[0].end(), [](const ScaffoldPart& lhs, const ScaffoldPart& rhs) {
            return (lhs.start - rhs.end - (rhs.start - lhs.end)) < 0;
        });
        newConsensus.numNode[1] = tmpConsensus1.numNode[1-h1] + tmpConsensus2.numNode[h2];
        newConsensus.node[1].resize(newConsensus.numNode[1]);
        for (long nodeIndex = 0; nodeIndex < tmpConsensus1.numNode[1-h1]; ++nodeIndex) {
            newConsensus.node[1][nodeIndex] = tmpConsensus1.node[1-h1][nodeIndex];
        }
        for (long nodeIndex = 0; nodeIndex < tmpConsensus2.numNode[h2]; ++nodeIndex) {
            newConsensus.node[1][tmpConsensus1.numNode[1-h1] + nodeIndex].id = tmpConsensus2.node[h2][nodeIndex].id * sign(oppositeConsensusID);
            if (oppositeConsensusID > 0) {
                newConsensus.node[1][tmpConsensus1.numNode[1-h1] + nodeIndex].start = tmpConsensus2.node[h2][nodeIndex].start + startConsensusNode;
                newConsensus.node[1][tmpConsensus1.numNode[1-h1] + nodeIndex].end = tmpConsensus2.node[h2][nodeIndex].end + startConsensusNode;
            } else {
                newConsensus.node[1][tmpConsensus1.numNode[1-h1] + nodeIndex].start = (tmpConsensus2.length - tmpConsensus2.node[h2][nodeIndex].end) + startConsensusNode;
                newConsensus.node[1][tmpConsensus1.numNode[1-h1] + nodeIndex].end = (tmpConsensus2.length - tmpConsensus2.node[h2][nodeIndex].start) + startConsensusNode;
            }
        }
        std::sort(newConsensus.node[1].begin(), newConsensus.node[1].end(), [](const ScaffoldPart& lhs, const ScaffoldPart& rhs) {
            return (lhs.start - rhs.end - (rhs.start - lhs.end)) < 0;
        });
        long min_start, max_end;
        if (newConsensus.numNode[0] < 1) {
            if (newConsensus.numNode[1] < 1) {
                std::cerr << "error" << std::endl; //
            } else {
                min_start = newConsensus.node[1].front().start;
                max_end = newConsensus.node[1].back().end;
            }
        } else{
            if (newConsensus.numNode[1] < 1) {
                min_start = newConsensus.node[0].front().start;
                max_end = newConsensus.node[0].back().end;
            } else {
                min_start = std::min(newConsensus.node[0].front().start, newConsensus.node[1].front().start);
                max_end = std::max(newConsensus.node[0].back().end, newConsensus.node[1].back().end);
            }
        }
        newConsensus.length = max_end - min_start;
        newConsensus.state = 0;

        long i = numConsensusNode - 1;
        for (long j = 0; j < consensusnode[i].numNode[0]; ++j) {
            long tmp = abs(consensusnode[i].node[0][j].id) - 1;
            if (consensusnode[i].node[0][j].id > 0)
                nodePositionInConsensus[tmp].emplace_back(0, platanus::Position(j, i + 1));
            else
                nodePositionInConsensus[tmp].emplace_back(0, platanus::Position(j, -(i + 1)));
        }
        for (long j = 0; j < consensusnode[i].numNode[1]; ++j) {
            long tmp = abs(consensusnode[i].node[1][j].id) - 1;
            if (consensusnode[i].node[1][j].id > 0)
                nodePositionInConsensus[tmp].emplace_back(1, platanus::Position(j, i + 1));
            else
                nodePositionInConsensus[tmp].emplace_back(1, platanus::Position(j, -(i + 1)));
        }

        oppositeConsensusInfo.resize(numConsensusNode);
        setOppositeBubbleConsensusNode(numConsensusNode-1, oppositeConsensusInfo[numConsensusNode-1]);
    }

    FILE* ScaffoldFP = platanus::makeTemporaryFile();
    long numNewConsensusNode = 0;
    long newNodePoolSize = 0;
    this->writeSingletonConsensusNode(numNewConsensusNode, newNodePoolSize, ScaffoldFP);
    this->remakeConsensus(numNewConsensusNode, newNodePoolSize, ScaffoldFP);
    fclose(ScaffoldFP);

}
//

//added by ouchi
void PairedDBG::makeConsensusGraph(const long numThread, bool allLib)
{
    if (allLib) {
        makeGraphAllLibraries(numThread);
    } else {
        makeGraph(numThread);
    }

    for (int i = 0; i < numConsensusNode; ++i) {
        consensusnode[i].edge.clear();
        consensusnode[i].numEdge = 0;
    }

    for (long nodeIndex = 0; nodeIndex < numNode; ++nodeIndex) {
        for (long edgeID = 0; edgeID < node[nodeIndex].numEdge; ++edgeID) {
            long numLink = node[nodeIndex].edge[edgeID].numLink;
            long nodeID2 = node[nodeIndex].edge[edgeID].end;
            for (long i = 0; i < nodePositionInConsensus[nodeIndex].size(); ++i) {
                long consensusID1 = nodePositionInConsensus[nodeIndex][i].second.id;
                long offset1 = nodePositionInConsensus[nodeIndex][i].second.offset;
                char h1 = nodePositionInConsensus[nodeIndex][i].first;
                for (long j = 0; j < nodePositionInConsensus[id2Index(nodeID2)].size(); ++j) {
                    long consensusID2 = nodePositionInConsensus[id2Index(nodeID2)][j].second.id;
                    long offset2 = nodePositionInConsensus[id2Index(nodeID2)][j].second.offset;
                    char h2 = nodePositionInConsensus[id2Index(nodeID2)][j].first;

                    if (abs(consensusID1) == abs(consensusID2))
                        continue;

                    long gap = node[nodeIndex].edge[edgeID].length;
                    long x,y;
                    if (node[nodeIndex].edge[edgeID].direction > 0) {
                        if (sign(consensusID1) > 0) {
                            x = consensusnode[id2Index(consensusID1)].length - consensusnode[id2Index(consensusID1)].node[h1][offset1].end;
                        } else {
                            x = consensusnode[id2Index(consensusID1)].node[h1][offset1].start;
                        }
                        if (sign(consensusID2) * sign(nodeID2)> 0) {
                            y = consensusnode[id2Index(consensusID2)].node[h2][offset2].start;
                        } else {
                            y = consensusnode[id2Index(consensusID2)].length - consensusnode[id2Index(consensusID2)].node[h2][offset2].end;
                        }
                    }else {
                        if (sign(consensusID1) > 0) {
                            x = consensusnode[id2Index(consensusID1)].node[h1][offset1].start;
                        } else {
                            x = consensusnode[id2Index(consensusID1)].length - consensusnode[id2Index(consensusID1)].node[h1][offset1].end;
                        }
                        if (sign(consensusID2) * sign(nodeID2) > 0) {
                            y = consensusnode[id2Index(consensusID2)].length - consensusnode[id2Index(consensusID2)].node[h2][offset2].end;
                        } else {
                            y = consensusnode[id2Index(consensusID2)].node[h2][offset2].start;
                        }
                    }
                    gap -= x + y;
                    //if (judgeConflictConsensusNode(0, gap + consensusnode[id2Index(consensusID1)].length, consensusID1, consensusID2, h1, h2, tolerence)) //
                    //    continue;
                    if (node[nodeIndex].edge[edgeID].direction > 0) {
                        if (judgeConflictConsensusNode(0, gap + consensusnode[id2Index(consensusID1)].length, consensusID1, sign(nodeID2) * consensusID2, h1, h2, tolerence))
                            continue;
                    } else {
                        if (judgeConflictConsensusNode(0, gap + consensusnode[id2Index(consensusID2)].length, sign(nodeID2) * consensusID2, consensusID1, h2, h1, tolerence))
                            continue;
                    }

                    long index1 = id2Index(consensusID1);
                    long m = 0;
                    for (; m < consensusnode[index1].numEdge; ++m) {
                        if (consensusnode[index1].edge[m].end == sign(consensusID1) * sign(nodeID2) * consensusID2 && consensusnode[index1].edge[m].direction == sign(consensusID1) * node[nodeIndex].edge[edgeID].direction)
                            break;
                    }
                    if (m == consensusnode[index1].numEdge) {
                        consensusnode[index1].edge.resize(m + 1);
                        ++(consensusnode[index1].numEdge);
                    }

                    consensusnode[index1].edge[m].numLink[h1 == h2] += numLink;
                    consensusnode[index1].edge[m].length += gap * numLink;
                    consensusnode[index1].edge[m].direction = sign(consensusID1) * node[nodeIndex].edge[edgeID].direction;
                    consensusnode[index1].edge[m].end = sign(consensusID1) * sign(nodeID2) * consensusID2;
                }
            }
        }
    }

    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        ConsensusNode &tmp = consensusnode[consensusIndex];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);

        for (long edgeIndex = 0; edgeIndex < consensusnode[consensusIndex].numEdge; ++edgeIndex) {
            consensusnode[consensusIndex].edge[edgeIndex].length /= consensusnode[consensusIndex].edge[edgeIndex].numLink[0] + consensusnode[consensusIndex].edge[edgeIndex].numLink[1];
        }
    }

}
//

//added by ouchi
void PairedDBG::getLinkedConsensusNode(const long sourceConsensusIndex, const char targetDirection, vector<ConsensusIDWithGap> &consensusIDBuffer)
{
	consensusIDBuffer.clear();
	ConsensusNode &consensusnode = this->consensusnode[sourceConsensusIndex];
	for (long edgeIndex = 0; edgeIndex < consensusnode.numEdge; ++ edgeIndex) {
		ConsensusEdge &edge = consensusnode.edge[edgeIndex];
		if (edge.direction == targetDirection) {
			consensusIDBuffer.emplace_back(edge.end, edge.length, edge.numLink[0] > edge.numLink[1], edge.numLink[0] + edge.numLink[1]);
		}
	}
}
//

//added by ouchi
long PairedDBG::joinUnambiguousConsensusPair(const long numThread)
{
    vector<vector<ConsensusPathGapped> > pathBuffer(numThread);

    # pragma omp parallel for schedule(static, 1)
	for (long t = 0; t < numThread; ++t) {
		vector<ConsensusIDWithGap> consensusBuffer;
		for (long initialConsensusIndex = t; initialConsensusIndex < numConsensusNode; initialConsensusIndex += numThread) {
			for (long i = 0; i < 2; ++i) {
				char initialDirection = 2*i -1;
				getLinkedConsensusNode(initialConsensusIndex, initialDirection, consensusBuffer);
				if (consensusBuffer.size() != 1)
					continue;
				if (consensusBuffer[0].numLink < minLink)
					continue;
				long altConsensusID = consensusBuffer[0].id;
				getLinkedConsensusNode(abs(altConsensusID) - 1, - (initialDirection * sign(altConsensusID)), consensusBuffer);
				if (consensusBuffer.size() != 1)
					continue;
				if (consensusBuffer[0].numLink < minLink)
					continue;

				pathBuffer[t].emplace_back(initialConsensusIndex, 2);
				if (i == 0) {
					pathBuffer[t].back().consensusnode[0].id = altConsensusID;
					pathBuffer[t].back().consensusnode[1].id = initialConsensusIndex + 1;
				} else {
					pathBuffer[t].back().consensusnode[0].id = initialConsensusIndex + 1;
					pathBuffer[t].back().consensusnode[1].id = altConsensusID;
				}
				pathBuffer[t].back().consensusnode[0].gap = 0;
				pathBuffer[t].back().consensusnode[1].gap = consensusBuffer[0].gap;
				pathBuffer[t].back().consensusnode[0].h = 0;
				pathBuffer[t].back().consensusnode[1].h = consensusBuffer[0].h;

				break;
				}
		}
	}

	vector<ConsensusPathGapped> mergedPathBuffer;
	mergeAndClearMultiThreadedVector(pathBuffer, mergedPathBuffer);
	std::stable_sort(mergedPathBuffer.begin(), mergedPathBuffer.end(), [](ConsensusPathGapped const& left, ConsensusPathGapped const& right) {
		return left.selfID < right.selfID;
	});

	long numJoin = remakeConsensusGraphAccordingToGappedPath(mergedPathBuffer);

	return numJoin;
}
//

//added by ouchi
void PairedDBG::joinUnambiguousConsensusPairIterative(const long numThread)
{
	long total = 0;
	long num;
	const unsigned long currentMinLink = this->minLink;
	cerr << endl << "joining unambiguous pair of nodes in a consensus graph.." << endl;
	do {
		setMinLink(1);
		makeConsensusGraph(numThread);
		setMinLink(currentMinLink);
		num = joinUnambiguousConsensusPair(numThread);
		total += num;
		cerr << "NUM_JOINED = " << num << endl;
	} while (num > 0);
	cerr << "TOTAL_NUM_NODES_IN_JOINED_PAIRS =" << total << endl << endl;
}
//

//added by ouchi
void PairedDBG::detectConflictConsensus(void) {
    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        for (long edgeID1 = 0; edgeID1 < consensusnode[consensusIndex].numEdge - 1; ++edgeID1) {
            for (long edgeID2 = edgeID1 + 1; edgeID2 < consensusnode[consensusIndex].numEdge; ++ edgeID2) {
                const ConsensusEdge &edge1 = (consensusnode[consensusIndex].edge[edgeID1]);
                const ConsensusEdge &edge2 = (consensusnode[consensusIndex].edge[edgeID2]);
                if (edge1.direction * edge2.direction < 0) continue;
                if (edge1.direction > 0 && (consensusnode[consensusIndex].state & SC_RR))
                    continue;
                if (edge1.direction < 0 && (consensusnode[consensusIndex].state & SC_RL))
                    continue;
                if (judgeConflictConsensusNode(edge1.length, edge2.length, edge1.end, edge2.end, edge1.numLink[0] > edge1.numLink[1], edge2.numLink[0] > edge2.numLink[1], tolerence)) {
                    if (edge1.direction > 0)
                        consensusnode[consensusIndex].state |= SC_RR;
                    else
                        consensusnode[consensusIndex].state |= SC_RL;
                }
            }
        }
    }
}
//

//added by ouchi
void PairedDBG::consensusScaffolding(void)
{
    vector<GraphLayout> include, candidate;
    long numNewConsensusNode = 0;
    long newNodePoolSize = 0;
    cerr << "consensus scaffolding..." << endl;

    FILE *scaffoldFP = platanus::makeTemporaryFile();

    for (long i = 0; i < numConsensusNode; ++i) {
//		cerr << "consensusIndex: " << i << endl;
        if (consensusnode[i].state & SC_INC) continue;
        include.clear();
        candidate.clear();
        include.push_back(GraphLayout());
        include[0].start = include[0].distance = 0;
        include[0].end = consensusnode[i].length;
        include[0].id = i + 1;
        include[0].h = 0;
        consensusnode[i].state |= SC_INC;
        for (long j = 0; j < consensusnode[i].numEdge; ++j) {
            if (consensusnode[i].edge[j].direction > 0 && (consensusnode[i].state & SC_RR))
                continue;
            if (consensusnode[i].edge[j].direction < 0 && (consensusnode[i].state & SC_RL))
                continue;
            long tmpConsensusIndex = abs(consensusnode[i].edge[j].end) - 1;
            if (consensusnode[tmpConsensusIndex].state & SC_INC)
                continue;
            if (-1 * consensusnode[i].edge[j].direction * consensusnode[i].edge[j].end > 0 && (consensusnode[tmpConsensusIndex].state & SC_RR))
                continue;
            if (-1 * consensusnode[i].edge[j].direction * consensusnode[i].edge[j].end < 0 && (consensusnode[tmpConsensusIndex].state & SC_RL))
                continue;
            candidate.push_back(GraphLayout());
            unsigned long long candidateEnd = candidate.size() - 1;
            if (consensusnode[i].edge[j].direction > 0) {
                candidate[candidateEnd].start = include[0].end + consensusnode[i].edge[j].length;
                candidate[candidateEnd].end = candidate[candidateEnd].start + consensusnode[tmpConsensusIndex].length;
            } else {
                candidate[candidateEnd].end = -(consensusnode[i].edge[j].length);
                candidate[candidateEnd].start = candidate[candidateEnd].end - consensusnode[tmpConsensusIndex].length;
            }
            candidate[candidateEnd].id = consensusnode[i].edge[j].end;
            candidate[candidateEnd].distance = 1;
            candidate[candidateEnd].h = consensusnode[i].edge[j].numLink[0] > consensusnode[i].edge[j].numLink[1];
        }

        while (candidate.size() > 0) {
            long minDistance = candidate[0].distance;
            long closest = candidate[0].start;
            long minCandidateID = 0;
            for (unsigned j = 1; j < candidate.size(); ++j) {
                if (candidate[j].distance < minDistance || (candidate[j].distance == minDistance && abs(candidate[j].start) < closest)) {
                    minDistance = candidate[j].distance;
                    closest = abs(candidate[j].start);
                    minCandidateID = j;
                }
            }

//		cerr << "minCandidateID: " << candidate[minCandidateID].id << endl;

            long tmpConsensusIndex = abs(candidate[minCandidateID].id) - 1;

            if (consensusnode[tmpConsensusIndex].state & SC_INC) {
                auto candidateIterator = candidate.begin();
                candidateIterator = candidate.erase(candidateIterator + minCandidateID);
                continue;
            }
            unsigned j = 0;
            for (; j < include.size(); ++j) {
                long overlapTolerence = std::min(tolerence, std::min(candidate[minCandidateID].end - candidate[minCandidateID].start, include[j].end - include[j].start) / 2);
                if (judgeConflictConsensusNode(include[j].start, candidate[minCandidateID].start, include[j].id, candidate[minCandidateID].id, include[j].h, candidate[minCandidateID].h, overlapTolerence))
                    break;
            }
            if (j == include.size()) {
                include.push_back(candidate[minCandidateID]);
                ConsensusNode &newConsensus = consensusnode[abs(include[include.size() - 1].id) - 1];
                for (long k = 0; k < newConsensus.numEdge; ++k) {
                    if (newConsensus.edge[k].direction > 0 && (newConsensus.state & SC_RR))
                        continue;
                    if (newConsensus.edge[k].direction < 0 && (newConsensus.state & SC_RL))
                        continue;
                    long tmpConsensusIndex = abs(newConsensus.edge[k].end) - 1;
                    if (consensusnode[tmpConsensusIndex].state & SC_INC)
                        continue;
                    if (-1 * newConsensus.edge[k].direction * newConsensus.edge[k].end > 0 && (consensusnode[tmpConsensusIndex].state & SC_RR))
                        continue;
                    if (-1 * newConsensus.edge[k].direction * newConsensus.edge[k].end < 0 && (consensusnode[tmpConsensusIndex].state & SC_RL))
                        continue;
                    GraphLayout tmpLayout;
                    if (include[include.size() - 1].id * newConsensus.edge[k].direction > 0) {
                        tmpLayout.start = include[include.size() - 1].end + newConsensus.edge[k].length;
                        tmpLayout.end = tmpLayout.start + consensusnode[tmpConsensusIndex].length;
                    } else {
                        tmpLayout.end = include[include.size() -1].start - newConsensus.edge[k].length;
                        tmpLayout.start = tmpLayout.end - consensusnode[tmpConsensusIndex].length;
                    }
                    tmpLayout.id = include[include.size() - 1].id > 0 ? newConsensus.edge[k].end : - (newConsensus.edge[k].end);
                    tmpLayout.distance = include[include.size() - 1].distance + 1;
                    tmpLayout.h = (include[include.size() - 1].h + (newConsensus.edge[k].numLink[0] > newConsensus.edge[k].numLink[1])) % 2;
                    candidate.push_back(tmpLayout);
                }
                newConsensus.state |= SC_INC;
            }
            auto candidateIterator = candidate.begin();
            candidateIterator = candidate.erase(candidateIterator + minCandidateID);
        }

//		cerr << "sort include" << endl;
        sort(include.begin(), include.end());
        long includeSize = include.size();
        long minStart = include[0].start;
        for (char h = 0; h < 2; ++h) {
            long numNewNode = 0;
            for (long j = 0; j < includeSize; ++j) {
                numNewNode += consensusnode[abs(include[j].id) - 1].numNode[(include[j].h + h) % 2];
            }
            fwrite(&numNewNode, sizeof(long), 1, scaffoldFP);
            newNodePoolSize += numNewNode;

            for (long j = 0; j < includeSize; ++j) {
                if (include[j].id > 0) {
                    for (long k = 0; k < consensusnode[abs(include[j].id) - 1].numNode[(include[j].h + h) % 2]; ++k) {
                        ScaffoldPart tmpScaffoldPart(consensusnode[abs(include[j].id ) - 1].node[(include[j].h + h) % 2][k].id, include[j].start - minStart + consensusnode[abs(include[j].id) - 1].node[(include[j].h + h) % 2][k].start, include[j].start - minStart + consensusnode[abs(include[j].id) - 1].node[(include[j].h + h) % 2][k].end);
                        fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                    }
                } else {
                    for (long k = consensusnode[abs(include[j].id) - 1].numNode[(include[j].h + h) % 2] - 1; k >= 0; --k) {
                        ScaffoldPart tmpScaffoldPart(-(consensusnode[abs(include[j].id ) - 1].node[(include[j].h + h) % 2][k].id), include[j].start - minStart + consensusnode[abs(include[j].id) - 1].length - consensusnode[abs(include[j].id) - 1].node[(include[j].h + h) % 2][k].end, include[j].start - minStart + consensusnode[abs(include[j].id) - 1].length - consensusnode[abs(include[j].id) - 1].node[(include[j].h + h) % 2][k].start);
                        fwrite(&tmpScaffoldPart, sizeof(ScaffoldPart), 1, scaffoldFP);
                    }
                }
            }
        }
        ++numNewConsensusNode;
    }
    include.clear();
    candidate.clear();

		cerr << "numNewConsensusNode: "<< numNewConsensusNode << endl;

    writeSingletonConsensusNode(numNewConsensusNode, newNodePoolSize, scaffoldFP);
    remakeConsensus(numNewConsensusNode, newNodePoolSize, scaffoldFP);
    fclose(scaffoldFP);
}
//

//added by ouchi
void PairedDBG::deleteConsensusEdges(vector<long> &ids)
{
    long id1, id2;
    long idsSize = ids.size();
    ConsensusNode *consensusPointer;
    for (long i = 0; i < idsSize; i+= 2) {
        id1 = ids[i];
        id2 = ids[i+1];
        if (id2 == 0)
            continue;

        consensusPointer = &(consensusnode[id1 - 1]);
        for (long j = 0; j < consensusPointer->numEdge; ++j) {
            if (consensusPointer->edge[j].end == id2) {
                consensusPointer->edge[j] = consensusPointer->edge[consensusPointer->numEdge - 1];
                --(consensusPointer->numEdge);
                break;
            }
        }

        if (id2 > 0) {
            consensusPointer = &(consensusnode[id2 - 1]);
            for (long j = 0; j < consensusPointer->numEdge; ++j) {
                if (consensusPointer->edge[j].end == id1) {
                    consensusPointer->edge[j] = consensusPointer->edge[consensusPointer->numEdge - 1];
                    --(consensusPointer->numEdge);
                    break;
                }
            }
        } else {
            consensusPointer = &(consensusnode[-id2 - 1]);
            for (long j = 0; j < consensusPointer->numEdge; ++j) {
                if (consensusPointer->edge[j].end == -id1) {
                    consensusPointer->edge[j] = consensusPointer->edge[consensusPointer->numEdge - 1];
                    --(consensusPointer->numEdge);
                    break;
                }
            }
        }
    }
    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        ConsensusNode &tmp = consensusnode[consensusIndex];
        std::sort(tmp.edge.begin(), tmp.edge.begin() + tmp.numEdge);
    }
}
//

//added by ouchi
void PairedDBG::deleteErroneousConsensusEdgebyHiC(const long numThread)
{
    vector<long> ids;
    long numDelete = 0;

    cerr << "removing erroneous edges using HiC reads ..." << endl << endl;

    makeHiCLink(numThread, 1);

    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        for (long edgeID = 0; edgeID < consensusnode[consensusIndex].numEdge; ++edgeID) {
            if (consensusnode[consensusIndex].length < 200000 || consensusnode[id2Index(consensusnode[consensusIndex].edge[edgeID].end)].length < 200000) {
                if (consensusnode[consensusIndex].edge[edgeID].numLink[0] + consensusnode[consensusIndex].edge[edgeID].numLink[1] >= minLink)
                    continue;
                {
                    ids.push_back(consensusIndex + 1);
                    ids.push_back(consensusnode[consensusIndex].edge[edgeID].end);
                }
                ++numDelete;
                continue;
            }
            if (checkHiCLinkBetweenConsensusPair(consensusIndex+1, consensusnode[consensusIndex].edge[edgeID].direction, consensusnode[consensusIndex].edge[edgeID].end))
                continue;
            {
                ids.push_back(consensusIndex + 1);
                ids.push_back(consensusnode[consensusIndex].edge[edgeID].end);
            }
            ++numDelete;
        }
    }
    deleteConsensusEdges(ids);

    cerr << "TOTAL_NUM_REMOVED_EDGES_BY_HIC =" << numDelete << endl << endl;
}
//

//added by ouchi
long PairedDBG::makeConsensusContactmap(const long consensusID1, const char direction, const long consensusID2, vector<vector<int> > &contactmap)
{
    const long binsize = 100000;
    const ConsensusNode &consensus1 = consensusnode[id2Index(consensusID1)];
    const ConsensusNode &consensus2 = consensusnode[id2Index(consensusID2)];

    if (consensus1.length / binsize < 1 || consensus2.length/binsize < 1)
        return -1;

    char cd2 = sign(consensusID2);

    long numbin = consensus1.length / binsize + consensus2.length / binsize;
    long lower, upper, boundary;
    if (direction > 0) {
        lower = consensus1.length % binsize;
        upper = consensus1.length + consensus2.length - consensus2.length % binsize;
        boundary = consensus1.length / binsize;
    } else {
        lower = consensus2.length % binsize;
        upper = consensus2.length + consensus1.length - consensus1.length % binsize;
        boundary = consensus2.length / binsize;
    }

    contactmap.resize(numbin);
    for (long i = 0; i < numbin; ++i) {
        contactmap[i].resize(numbin);
    }

    long offset1, offset2;
    for (char h = 0; h < 2; ++h) {
        for (long nodeID = 0; nodeID < consensus1.numNode[h]; ++nodeID) {
            const GraphNode &tmp = node[id2Index(consensus1.node[h][nodeID].id)];
            for (long HiCLinkID = 0; HiCLinkID < tmp.HiCLinks.size(); ++HiCLinkID) {
                const HiCLink &hiclink = tmp.HiCLinks[HiCLinkID];
                offset1 = consensus1.node[h][nodeID].start;
                if (sign(consensus1.node[h][nodeID].id) > 0) {
                    offset1 += hiclink.offset1;
                } else {
                    offset1 += tmp.length - hiclink.offset1;
                }
                if (direction < 0)
                    offset1 += consensus2.length;

                if (nodePositionInConsensus[id2Index(hiclink.end)].size() > 0) {
                    long tmpConsensusID = nodePositionInConsensus[id2Index(hiclink.end)][0].second.id;
                    long tmpOffset = nodePositionInConsensus[id2Index(hiclink.end)][0].second.offset;
                    long tmph = nodePositionInConsensus[id2Index(hiclink.end)][0].first;

                    if (nodePositionInConsensus[id2Index(consensus1.node[h][nodeID].id)].size() == 2 && nodePositionInConsensus[id2Index(hiclink.end)].size() == 1) {
                        if (h != tmph)
                            continue;
                    }

                    if (abs(tmpConsensusID) == abs(consensusID1)) {
                        offset2 = consensus1.node[tmph][tmpOffset].start;
                        if (sign(tmpConsensusID) > 0) {
                            offset2 += hiclink.offset2;
                        } else {
                            offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        }
                        if (direction < 0)
                            offset2 += consensus2.length;
                    } else if (abs(tmpConsensusID) == abs(consensusID2)) {
                        if (cd2 > 0) {
                            offset2 = consensus2.node[tmph][tmpOffset].start;
                        } else {
                            offset2 = consensus2.length - consensus2.node[tmph][tmpOffset].end;
                        }
                        if (cd2 * sign(tmpConsensusID) > 0) {
                            offset2 += hiclink.offset2;
                        } else {
                            offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        }
                        if (direction > 0)
                            offset2 += consensus1.length;
                    } else {
                        continue;
                    }
                    if (lower <= offset1 && offset1 < upper && lower <= offset2 && offset2 < upper) {
                        ++contactmap[(offset1-lower)/binsize][(offset2-lower)/binsize];
                    }
                }
            }
        }
        for (long nodeID = 0; nodeID < consensus2.numNode[h]; ++nodeID) {
            const GraphNode &tmp = node[id2Index(consensus2.node[h][nodeID].id)];
            for (long HiCLinkID = 0; HiCLinkID < tmp.HiCLinks.size(); ++HiCLinkID) {
                const HiCLink &hiclink = tmp.HiCLinks[HiCLinkID];
                if (cd2 > 0) {
                    offset1 = consensus2.node[h][nodeID].start;
                } else {
                    offset1 = consensus2.length - consensus2.node[h][nodeID].end;
                }
                if (cd2 * sign(consensus2.node[h][nodeID].id) > 0) {
                    offset1 += hiclink.offset1;
                } else {
                    offset1 += tmp.length - hiclink.offset1;
                }
                if (direction > 0)
                    offset1 += consensus1.length;

                if (nodePositionInConsensus[id2Index(hiclink.end)].size() > 0) {
                    long tmpConsensusID = nodePositionInConsensus[id2Index(hiclink.end)][0].second.id;
                    long tmpOffset = nodePositionInConsensus[id2Index(hiclink.end)][0].second.offset;
                    long tmph = nodePositionInConsensus[id2Index(hiclink.end)][0].first;

                    if (nodePositionInConsensus[id2Index(consensus2.node[h][nodeID].id)].size() == 2 && nodePositionInConsensus[id2Index(hiclink.end)].size() == 1) {
                        if (h != tmph)
                            continue;
                    }

                    if (abs(tmpConsensusID) == abs(consensusID1)) {
                        offset2 = consensus1.node[tmph][tmpOffset].start;
                        if (sign(tmpConsensusID) > 0) {
                            offset2 += hiclink.offset2;
                        } else {
                            offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        }
                        if (direction < 0)
                            offset2 += consensus2.length;
                    } else if (abs(tmpConsensusID) == abs(consensusID2)) {
                        if (cd2 > 0) {
                            offset2 = consensus2.node[tmph][tmpOffset].start;
                        } else {
                            offset2 = consensus2.length - consensus2.node[tmph][tmpOffset].end;
                        }
                        if (cd2 * sign(tmpConsensusID) > 0) {
                            offset2 += hiclink.offset2;
                        } else {
                            offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        }
                        if (direction > 0)
                            offset2 += consensus1.length;
                    } else {
                        continue;
                    }
                    if (lower <= offset1 && offset1 < upper && lower <= offset2 && offset2 < upper) {
                        ++contactmap[(offset1-lower)/binsize][(offset2-lower)/binsize];
                    }
                }
            }
        }

    }

	std::cerr << "contactmap" << std::endl;
	std::cerr << "consensusID1:" << consensusID1 << "\tdirection:" << direction << "\tconsensusID2:" << consensusID2 << std::endl;
	for (long i = 0; i < contactmap.size(); ++i) {
		if (i == boundary) {
			for (long j = 0; j < contactmap[i].size(); ++j) {
				std::cerr << "_\t";
			}
			std::cerr << std::endl;
		}
		for (long j = 0; j < contactmap[i].size(); ++j) {
			if (j == boundary) std::cerr << "|\t";
			std::cerr << contactmap[i][j] << "\t";
		}
		std::cerr << std::endl;
	}

    return boundary;
}
//

//added by ouchi
bool PairedDBG::checkHiCLinkBetweenConsensusPair(const long consensusID1, const char direction, const long consensusID2)
{
    vector<vector<int> > contactmap;
    long boundary = makeConsensusContactmap(consensusID1, direction, consensusID2, contactmap);
    if (boundary == -1)
        return 1;
    deleteGapRegionFromContactmap(contactmap, boundary);
	std::cerr << contactmap.size() << "\t" << boundary << std::endl;
    std::array<double, 6> result;
    calcSeparation(contactmap, 0, boundary, contactmap.size(), result, true);
	std::cerr << "result: " << result[0]/result[1] << "\t" << result[3]/result[4] << "\t" << result[2] << "\t" << result[5] << std::endl;
	std::cerr << (result[0]/result[1] + result[3]/result[4] < HiC_CONTACTMAP_S_THRESHOLD || result[2] + result[5] < HiC_CONTACTMAP_T_THRESHOLD) << std::endl;
    return (result[0]/result[1] + result[3]/result[4] < HiC_CONTACTMAP_S_THRESHOLD || result[2] + result[5] < HiC_CONTACTMAP_T_THRESHOLD);
}
//

//added by ouchi
double PairedDBG::calcConsensusCoverage(const ConsensusNode& consensusnode)
{
    std::array<unsigned long long, 2> sum;
    std::array<unsigned long long, 2> num;
    std::array<double, 2> coverage;
    for (char h = 0; h < 2; ++h) {
        sum[h] = 0; num[h] = 0;
        for (long i = 0; i < consensusnode.numNode[h]; ++i) {
            long nodeIndex = id2Index(consensusnode.node[h][i].id);
            double tmpCoverage = calcNodeCoverage(node[nodeIndex]) / nodePositionInConsensus[nodeIndex].size();
            num[h] += node[nodeIndex].length;
            sum[h] += tmpCoverage * node[nodeIndex].length;
        }
        if (num[h] == 0)
            coverage[h] = 0.0;
        else
            coverage[h] = static_cast<double>(sum[h] / num[h]);
    }
    return coverage[0] + coverage[1];
}
//

//added by ouchi
long PairedDBG::checkConsensusPath(long start, long end, long maxDistance)
{
	std::cerr << "checkConsensusPath " << start << " to " << end << std::endl;
    vector<long> candiate;
    ConsensusNode &tmpConsensus = consensusnode[id2Index(start)];
    for (long edgeIndex = 0; edgeIndex < tmpConsensus.numEdge; ++edgeIndex) {
        ConsensusEdge &tmpEdge = tmpConsensus.edge[edgeIndex];
        if (tmpEdge.direction == sign(start)) {
            if (tmpEdge.direction * tmpEdge.end == end) {
                candiate.push_back(std::max(tmpEdge.numLink[0], tmpEdge.numLink[1]));
				std::cerr << "candiate.push_back " << std::max(tmpEdge.numLink[0], tmpEdge.numLink[1]) << std::endl;
            } else {
                if (maxDistance > 1) {
                    long tmpNumLink = checkConsensusPath(tmpEdge.direction * tmpEdge.end, end, maxDistance-1);
					std::cerr << "tmpNumLink: " << tmpNumLink << std::endl;
                    if (tmpNumLink > 0) {
                        candiate.push_back(std::min(std::max(tmpEdge.numLink[0], tmpEdge.numLink[1]), tmpNumLink));
						std::cerr << "candiate.push_back2 " << std::min(std::max(tmpEdge.numLink[0], tmpEdge.numLink[1]), tmpNumLink) << std::endl;
                    }
                }
            }
        }
    }
    if (candiate.empty()) {
        return 0;
    } else {
        long maxLink = 0;
        for (long i = 0; i < candiate.size(); ++i) {
            if (candiate[i] > maxLink)
                maxLink = candiate[i];
        }
        return maxLink;
    }
}
//

/*
//added by ouchi
void PairedDBG::makeConsensusPathEdge(const long numThread)
{
    for (int i = 0; i < numConsensusNode; ++i) {
        consensusnode[i].pathEdge.clear();
        consensusnode[i].numPathEdge = 0;
    }

    vector<vector<long> > tmpPath;
    vector<vector<char> > tmph;
    vector<long> tmpLength;
    vector<std::array<long, 2> > tmpNumLink;
    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        searchConsensusPath(consensusIndex + 1, tmpPath, tmph, tmpLength, tmpNumLink);
        searchConsensusPath(-consensusIndex - 1, tmpPath, tmph, tmpLength, tmpNumLink);
    }

   for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        std::cerr << "consensusIndex: " << consensusIndex << std::endl;
        for (long pathEdgeIndex = 0; pathEdgeIndex < consensusnode[consensusIndex].numPathEdge; ++pathEdgeIndex) {
            std::cerr << consensusnode[consensusIndex].pathEdge[pathEdgeIndex].end << "\t" << consensusnode[consensusIndex].pathEdge[pathEdgeIndex].length << std::endl;
        }
   }

}
//

//added by ouchi
void PairedDBG::searchConsensusPath(long startID, vector<vector<long> > &path, vector<vector<char> > &h, vector<long> &lengths, vector<std::array<long, 2> > &numLinks)
{
	std::cerr << "searchConsensusPath " << startID << std::endl;
    ConsensusNode& tmpConsensus = this->consensusnode[id2Index(startID)];
    if ((tmpConsensus.searchState & 0x1) && startID > 0) {
        for (long pathEdgeIndex = 0; pathEdgeIndex < tmpConsensus.numPathEdge; ++pathEdgeIndex) {
            if (tmpConsensus.pathEdge[pathEdgeIndex].direction == sign(startID)) {
                path.push_back(tmpConsensus.pathEdge[pathEdgeIndex].path[0]);
                h.push_back(tmpConsensus.pathEdge[pathEdgeIndex].h[0]);
                lengths.push_back(tmpConsensus.pathEdge[pathEdgeIndex].length);
                numLinks.push_back(tmpConsensus.pathEdge[pathEdgeIndex].numLink);
            }
        }
        return;
    } else if ((tmpConsensus.searchState & 0x2) && startID < 0) {
        for (long pathEdgeIndex = 0; pathEdgeIndex < tmpConsensus.numPathEdge; ++pathEdgeIndex) {
            if (tmpConsensus.pathEdge[pathEdgeIndex].direction == sign(startID)) {
                path.push_back(tmpConsensus.pathEdge[pathEdgeIndex].path[0]);
                h.push_back(tmpConsensus.pathEdge[pathEdgeIndex].h[0]);
                lengths.push_back(tmpConsensus.pathEdge[pathEdgeIndex].length);
                numLinks.push_back(tmpConsensus.pathEdge[pathEdgeIndex].numLink);
            }
        }
        return;
    } else {
        if (startID > 0) {
            if (tmpConsensus.searchState & 0x4)
                return;
            else
                tmpConsensus.searchState |= 0x4;
        } else {
            if (tmpConsensus.searchState & 0x8)
                return;
            else
                tmpConsensus.searchState |= 0x8;
        }
    }


    for (long edgeIndex = 0; edgeIndex < tmpConsensus.numEdge; ++edgeIndex) {
		std::cerr << "edgeIndex: " << edgeIndex<<std::endl;
        if (sign(startID) == tmpConsensus.edge[edgeIndex].direction) {
            long endID = sign(startID) * tmpConsensus.edge[edgeIndex].end;
            long index = tmpConsensus.numPathEdge;
            tmpConsensus.pathEdge.resize(index + 1);
            ++tmpConsensus.numPathEdge;
            tmpConsensus.pathEdge[index].direction = sign(startID);
            tmpConsensus.pathEdge[index].end = endID;
            tmpConsensus.pathEdge[index].length = tmpConsensus.edge[edgeIndex].length;
            tmpConsensus.pathEdge[index].numLink = tmpConsensus.edge[edgeIndex].numLink;
            tmpConsensus.pathEdge[index].path.resize(1);
            tmpConsensus.pathEdge[index].path[0] = {startID, endID};
            tmpConsensus.pathEdge[index].h.resize(1);
            tmpConsensus.pathEdge[index].h[0] = {0, (tmpConsensus.edge[edgeIndex].numLink[0] > tmpConsensus.edge[edgeIndex].numLink[1])};
            path.push_back(tmpConsensus.pathEdge[index].path[0]);
            h.push_back(tmpConsensus.pathEdge[index].h[0]);
            lengths.push_back(tmpConsensus.pathEdge[index].length);
            numLinks.push_back(tmpConsensus.pathEdge[index].numLink);
			std::cerr << "\t" << endID << std::endl;
            vector<vector<long> > tmpPath;
            vector<vector<char> > tmph;
            vector<long> tmpLength;
            vector<std::array<long, 2> > tmpNumLink;
            searchConsensusPath(endID, tmpPath, tmph, tmpLength, tmpNumLink);
			std::cerr << "search fin\t" << endID << std::endl;
            for (long i = 0; i < tmpPath.size(); ++i) {
                if (tmpPath[i].size() == 0)
                    continue;
				std::cerr << "pathEdge.end: " << tmpPath[i].back();
				std::cerr << "\t" << "length: " << tmpLength[i];
				std::cerr << "\tnumLink:" << tmpNumLink[i][0];
				std::cerr << "\t" << tmpNumLink[i][1] << std::endl;
                index = tmpConsensus.numPathEdge;
                tmpConsensus.pathEdge.resize(index + 1);
                ++tmpConsensus.numPathEdge;
                tmpConsensus.pathEdge[index].direction = sign(startID);
                tmpConsensus.pathEdge[index].end = tmpPath[i].back();
                tmpConsensus.pathEdge[index].length = tmpLength[i] + consensusnode[id2Index(endID)].length + tmpConsensus.edge[edgeIndex].length;
                tmpConsensus.pathEdge[index].numLink[0] = std::min(tmpConsensus.edge[edgeIndex].numLink[0], tmpNumLink[i][(tmpConsensus.edge[edgeIndex].numLink[0] > tmpConsensus.edge[edgeIndex].numLink[1])]);
                tmpConsensus.pathEdge[index].numLink[1] = std::min(tmpConsensus.edge[edgeIndex].numLink[1], tmpNumLink[i][(1 + (tmpConsensus.edge[edgeIndex].numLink[0] > tmpConsensus.edge[edgeIndex].numLink[1])) % 2]);
                tmpConsensus.pathEdge[index].path.resize(1);
                tmpConsensus.pathEdge[index].path[0] = tmpPath[i];
                tmpConsensus.pathEdge[index].path[0].insert(tmpConsensus.pathEdge[index].path[0].begin(), startID);
                tmpConsensus.pathEdge[index].h.resize(1);
                tmpConsensus.pathEdge[index].h[0].push_back(0);
                for (long j = 0; j < tmph[i].size(); ++j) {
                    tmpConsensus.pathEdge[index].h[0].push_back((tmpConsensus.edge[edgeIndex].numLink[0] > tmpConsensus.edge[edgeIndex].numLink[1] + tmph[i][j]) % 2);
                }
                path.push_back(tmpConsensus.pathEdge[index].path[0]);
                h.push_back(tmpConsensus.pathEdge[index].h[0]);
                lengths.push_back(tmpConsensus.pathEdge[index].length);
                numLinks.push_back(tmpConsensus.pathEdge[index].numLink);
            }
        }
    }

    if (startID > 0) {
        tmpConsensus.searchState &= ~0x4;
        tmpConsensus.searchState |= 0x1;
    } else {
        tmpConsensus.searchState &= ~0x8;
        tmpConsensus.searchState |= 0x2;
    }

    return;
}
// */


//added by ouchi
void PairedDBG::HiC_Scaffolding(const long numThread)
{
    const unsigned long currentMinLink = this->minLink;
    const long currentCutoffLength = this->cutoffLength;
    setMinLink(1);
    setCutoffLength(0);

    makeConsensusGraph(numThread);
	std::cerr << "makeHiCLink" << std::endl;
    makeHiCLink(numThread, 2);
	std::cerr<< "finish makeHiCLink" << std::endl;

    setMinLink(currentMinLink);
    setCutoffLength(currentCutoffLength);

    //HiCNodeGraph hg;
    //hg.initGraph(numConsensusNode);

    long numHiCNode;
    std::vector<HiCNode> hicnode;
    std::vector<platanus::Position> consensusPositionInHiCNode;

    numHiCNode = numConsensusNode;
    hicnode.resize(numHiCNode);
    consensusPositionInHiCNode.resize(numConsensusNode);
	std::cerr << "numHiCNode:" << numHiCNode << std::endl;
    for (long i = 0; i < numHiCNode; ++i) {
        hicnode[i].numScaffold = 1;
        hicnode[i].scaffold.resize(1);
        hicnode[i].scaffold[0].id = i + 1;
        hicnode[i].scaffold[0].end = consensusnode[i].length;
        hicnode[i].length = consensusnode[i].length;
        consensusPositionInHiCNode[i].id = i + 1;
        consensusPositionInHiCNode[i].offset = 0;
    }

    vector<long> Llist = {10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000};
    for (long i = 0; i < Llist.size(); ++i) {
        long L = Llist[i];
        vector<HiCConsensusEdge> es;
		std::cerr << "L: " << L << std::endl;

		long numChecked = 0;
		long numContactmapChecked = 0;
		long numContactmapTrue = 0;
		long numContactmapFalse = 0;
		long numMatePairTrue = 0;
		long numMatePairFalse = 0;
		long numContactmapAndMateTrue = 0;
		long numContactmapAndMateFalse = 0;

        calcHiCLinkScore(L, es, numHiCNode, hicnode, consensusPositionInHiCNode);
        for (long ei = 0; ei < es.size(); ++ei) {
            HiCConsensusEdge &e = es[ei];

            long hicnodeID1 = consensusPositionInHiCNode[id2Index(e.id1)].id;
            long hicnodeID2 = consensusPositionInHiCNode[id2Index(e.id2)].id;
            //judgeSameScaffold
            if (abs(hicnodeID1) == abs(hicnodeID2))
                continue;
            vector<ConsensusPart> scaffolda, scaffoldb, scaffolda2, scaffoldb2;
            getConsensusPart(e.id1, 0, 1, scaffolda, numHiCNode, hicnode, consensusPositionInHiCNode);
            getConsensusPart(e.id2, e.h, 0, scaffoldb, numHiCNode, hicnode, consensusPositionInHiCNode);
            getConsensusPart(e.id1, 0, 2, scaffolda2, numHiCNode, hicnode, consensusPositionInHiCNode);
            getConsensusPart(e.id2, e.h, 3, scaffoldb2, numHiCNode, hicnode, consensusPositionInHiCNode);

            //checkEndOfScaffold
            if (scaffolda2.size() != 0 && scaffoldb2.size() != 0) continue;

            if (L < 100000) {
                if (e.numLink < 50) break;
            } else {
                if (e.numLink < 100) break;
            }

            std::cerr << "Score:" << e.numLink << "\n";
            for (long i  = 0; i < scaffolda.size(); ++i) {
                std::cerr << abs(scaffolda[i].id);
                if (scaffolda[i].id > 0)
                    std::cerr << "+";
                else
                    std::cerr << "-";
                std::cerr << abs(scaffolda[i].id);
                if (scaffolda[i].id > 0)
                    std::cerr << "-";
                else
                    std::cerr << "+";
            }
            std::cerr << "\t";
            for (long i  = 0; i < scaffoldb.size(); ++i) {
                std::cerr << abs(scaffoldb[i].id);
                if (scaffoldb[i].id > 0)
                    std::cerr << "+";
                else
                    std::cerr << "-";
                std::cerr << abs(scaffoldb[i].id);
                if (scaffoldb[i].id > 0)
                    std::cerr << "-";
                else
                    std::cerr << "+";
            }
            std::cerr << std::endl;

			//
			++numChecked;
			long length1 = 0; long length2 = 0;
			for (long i = 0; i < scaffolda.size(); ++i) {
				length1 += consensusnode[id2Index(scaffolda[i].id)].length;
			}
			for (long i = 0; i < scaffoldb.size(); ++i) {
				length2 += consensusnode[id2Index(scaffoldb[i].id)].length;
			}
			if (!(length1 / 100000 < 1 || length2 / 100000 < 1))
				++numContactmapChecked;
			//

            //checkLink id1 id2
            long numLink = checkConsensusPath(-e.id1, e.id2, 2);
            std::cerr << "MatePairLink: " << numLink << std::endl;

			//
			if (numLink == 0)
				++numMatePairFalse;
			else
				++numMatePairTrue;
			bool result = checkHiCLinkBetweenHiCNodePair(scaffolda, scaffoldb);
			if (result) {
				if (!(length1 / 100000 < 1 || length2 / 100000 < 1))
					++numContactmapTrue;
			} else {
				++numContactmapFalse;
			}
			if (numLink > 0 && result)
				if (!(length1 / 100000 < 1 || length2 / 100000 < 1))
					++numContactmapAndMateTrue;
			if (numLink == 0 && !result)
				++numContactmapAndMateFalse;
			if (!(result)) continue;
			//
//            if (!(checkHiCLinkBetweenHiCNodePair(scaffolda, scaffoldb))) continue;
            if (numLink == 0 && (length1/100000<1 || length2/100000 <1)) continue;

            //compare
            if (scaffolda2.size() != 0) {
                std::cerr << "compare: ";
                for (long i  = 0; i < scaffolda.size(); ++i) {
                    std::cerr << abs(scaffolda[i].id);
                    if (scaffolda[i].id > 0)
                        std::cerr << "+";
                    else
                        std::cerr << "-";
                    std::cerr << abs(scaffolda[i].id);
                    if (scaffolda[i].id > 0)
                        std::cerr << "-";
                    else
                        std::cerr << "+";
                }
                std::cerr << "\t";
                for (long i  = 0; i < scaffolda2.size(); ++i) {
                    std::cerr << abs(scaffolda2[i].id);
                    if (scaffolda2[i].id > 0)
                        std::cerr << "+";
                    else
                        std::cerr << "-";
                    std::cerr << abs(scaffolda2[i].id);
                    if (scaffolda2[i].id > 0)
                        std::cerr << "-";
                    else
                        std::cerr << "+";
                }
                std::cerr << "\t";
                for (long i  = 0; i < scaffoldb.size(); ++i) {
                    std::cerr << abs(scaffoldb[i].id);
                    if (scaffoldb[i].id > 0)
                        std::cerr << "+";
                    else
                        std::cerr << "-";
                    std::cerr << abs(scaffoldb[i].id);
                    if (scaffoldb[i].id > 0)
                        std::cerr << "-";
                    else
                        std::cerr << "+";
                }
                std::cerr << std::endl;
				continue;
            }
            if (scaffoldb2.size() != 0) {
                std::cerr << "sompare: ";
                for (long i  = 0; i < scaffoldb.size(); ++i) {
                    std::cerr << abs(scaffoldb[i].id);
                    if (scaffoldb[i].id > 0)
                        std::cerr << "+";
                    else
                        std::cerr << "-";
                    std::cerr << abs(scaffoldb[i].id);
                    if (scaffoldb[i].id > 0)
                        std::cerr << "-";
                    else
                        std::cerr << "+";
                }
                std::cerr << "\t";
                for (long i  = 0; i < scaffoldb2.size(); ++i) {
                    std::cerr << abs(scaffoldb2[i].id);
                    if (scaffoldb2[i].id > 0)
                        std::cerr << "+";
                    else
                        std::cerr << "-";
                    std::cerr << abs(scaffoldb2[i].id);
                    if (scaffoldb2[i].id > 0)
                        std::cerr << "-";
                    else
                        std::cerr << "+";
                }
                std::cerr << "\t";
                for (long i  = 0; i < scaffolda.size(); ++i) {
                    std::cerr << abs(scaffolda[i].id);
                    if (scaffolda[i].id > 0)
                        std::cerr << "+";
                    else
                        std::cerr << "-";
                    std::cerr << abs(scaffolda[i].id);
                    if (scaffolda[i].id > 0)
                        std::cerr << "-";
                    else
                        std::cerr << "+";
                }
                std::cerr << std::endl;
				continue;
            }

            //unite
            uniteHiCNode(e.id1, e.id2, e.h, numHiCNode, hicnode, consensusPositionInHiCNode);
        }

		//
		std::cerr << "Checked: " << numChecked << std::endl;
		std::cerr << "ContactmapChecked: " << numContactmapChecked << "\tTrue: " << numContactmapTrue << "\tFalse: " << numContactmapFalse << std::endl;
		std::cerr << "MatePair True: " << numMatePairTrue << "\tFalse: " << numMatePairFalse << std::endl;
		std::cerr << "Contactmap&MatePair True: " << numContactmapAndMateTrue << "\tFalse: " << numContactmapAndMateFalse << std::endl;
		//

        std::cerr << "Output HiCNode" << std::endl;
        std::string outputName = "HiCNode_result_" + std::to_string(L) + ".fa";
        outputHiCNode(outputName, numHiCNode, hicnode);

    }
}
//

//added by ouchi
void PairedDBG::calcHiCLinkScore(const long L, vector<HiCConsensusEdge> &es, long &numHiCNode, vector<HiCNode> &hicnode, vector<platanus::Position> &consensusPositionInHiCNode) {
    vector<vector<vector<long> > > info(numConsensusNode);

    for (long i = 0; i < numHiCNode; ++i) {
        if (hicnode[i].state & SC_DEL) continue;
        long l, id;
        long num = hicnode[i].numScaffold;
        if (num == 0) continue;
        for (long j = 0; j < num; ++j) {
            id = id2Index(hicnode[i].scaffold[j].id);
            if (consensusnode[id].length >= L) continue;
            l = L - consensusnode[id].length;
            for (long k = j-1; k >= 0; --k) {
                vector<long> tmp(4);
                tmp[0] = -1 * hicnode[i].scaffold[j].id;
                tmp[1] = l;
                tmp[2] = -1 * sign(hicnode[i].scaffold[k].id);
                tmp[3] = (hicnode[i].scaffold[j].h == hicnode[i].scaffold[k].h);
                info[id2Index(hicnode[i].scaffold[k].id)].push_back(tmp);
                l = l - consensusnode[id2Index(hicnode[i].scaffold[k].id)].length;
                if (l <= 0) break;
            }
            l = L - consensusnode[id].length;
            for (long k = j+1; k < num; ++k) {
                vector<long> tmp(4);
                tmp[0] = hicnode[i].scaffold[j].id;
                tmp[1] = l;
                tmp[2] = sign(hicnode[i].scaffold[k].id);
                tmp[3] = (hicnode[i].scaffold[j].h == hicnode[i].scaffold[k].h);
                info[id2Index(hicnode[i].scaffold[k].id)].push_back(tmp);
                l = l - consensusnode[id2Index(hicnode[i].scaffold[k].id)].length;
                if (l <= 0) break;
            }
        }
    }

    vector<ConsensusLink> graphLinkPool;

	# pragma omp parallel for schedule(dynamic)
    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
        const ConsensusNode& tmpConsensus = consensusnode[consensusIndex];
        if (tmpConsensus.length < cutoffLength)
            continue;
        for (char h = 0; h < 2; ++h) {
            for (long nodeID = 0; nodeID < tmpConsensus.numNode[h]; ++nodeID) {
                const GraphNode& tmp = node[id2Index(tmpConsensus.node[h][nodeID].id)];
                for (long HiCLinkID = 0; HiCLinkID < tmp.HiCLinks.size(); ++HiCLinkID) {
                    const HiCLink& hiclink = tmp.HiCLinks[HiCLinkID];
                    long offset1, offset2;
                    vector<long> startIDs, endIDs;
                    vector<char> starths, endhs;
                    offset1 = tmpConsensus.node[h][nodeID].start;
                    if (sign(tmpConsensus.node[h][nodeID].id) > 0) {
                        offset1 += hiclink.offset1;
                    } else {
                        offset1 += tmp.length - hiclink.offset1;
                    }

                    if (offset1 <= L) {
                        startIDs.push_back(consensusIndex + 1);
                        starths.push_back(h);
                    } else if (offset1 >= tmpConsensus.length - L) {
                        startIDs.push_back(- (consensusIndex + 1));
                        starths.push_back(h);
                    }
                    for (long j = 0; j < info[consensusIndex].size(); ++j) {
                        if (info[consensusIndex][j][2] > 0) {
                            if (offset1 <= info[consensusIndex][j][1]) {
                                startIDs.push_back(info[consensusIndex][j][0]);
                                starths.push_back((h + info[consensusIndex][j][3]) % 2);
                            }
                        } else {
                            if (offset1 >= tmpConsensus.length - info[consensusIndex][j][1]) {
                                startIDs.push_back(info[consensusIndex][j][0]);
                                starths.push_back((h + info[consensusIndex][j][3]) % 2);
                            }
                        }
                    }

                    if (startIDs.empty())
                        continue;

                    for (long i = 0; i < nodePositionInConsensus[id2Index(hiclink.end)].size(); ++i) {
                        long tmpConsensusID = nodePositionInConsensus[id2Index(hiclink.end)][i].second.id;
                        long tmpOffset = nodePositionInConsensus[id2Index(hiclink.end)][i].second.offset;
                        char tmph = nodePositionInConsensus[id2Index(hiclink.end)][i].first;

                        if (abs(tmpConsensusID) == consensusIndex + 1)
                            continue;
                        const ConsensusNode& tmpConsensus2 = consensusnode[id2Index(tmpConsensusID)];
                        if (tmpConsensus2.length < cutoffLength)
                            continue;
                        offset2 = tmpConsensus2.node[tmph][tmpOffset].start;
                        if (sign(tmpConsensus2.node[tmph][tmpOffset].id) > 0) {
                            offset2 += hiclink.offset2;
                        } else {
                            offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        }

                        if (offset2 <= L) {
                            endIDs.push_back(abs(tmpConsensusID));
                            endhs.push_back(tmph);
                        } else if (offset2 >= tmpConsensus2.length - L) {
                            endIDs.push_back(- abs(tmpConsensusID));
                            endhs.push_back(tmph);
                        }
                        for (long j = 0; j < info[id2Index(tmpConsensusID)].size(); ++j) {
                            if (info[id2Index(tmpConsensusID)][j][2] > 0) {
                                if (offset2 <= info[id2Index(tmpConsensusID)][j][1]) {
                                    endIDs.push_back(info[id2Index(tmpConsensusID)][j][0]);
                                    endhs.push_back((tmph + info[id2Index(tmpConsensusID)][j][3]) % 2);
                                }
                            } else {
                                if (offset2 >= tmpConsensus2.length - info[id2Index(tmpConsensusID)][j][1]) {
                                    endIDs.push_back(info[id2Index(tmpConsensusID)][j][0]);
                                    endhs.push_back((tmph + info[id2Index(tmpConsensusID)][j][3]) % 2);
                                }
                            }
                        }
                    }

					# pragma omp critical (push)
					{
                        for (long j = 0; j < startIDs.size(); ++j) {
                            for (long k = 0; k < endIDs.size(); ++k) {
                                ConsensusLink graphLink;
                                graphLink.id1 = startIDs[j];
                                graphLink.h1 = starths[j];
                                graphLink.id2 = endIDs[k];
                                graphLink.h2 = endhs[k];
                                if (abs(graphLink.id1) < abs(graphLink.id2)) {
                                    graphLinkPool.push_back(graphLink);
                                }
                            }
                        }
					}

                }
            }
        }
    }

    std::stable_sort(graphLinkPool.begin(), graphLinkPool.end());

    vector<ConsensusLinkPoolIndex> indexes(1, 0);
    ++indexes.back().numLinks[graphLinkPool[0].h1 == graphLinkPool[0].h2];
    for (long idx = 1; idx < graphLinkPool.size(); ++idx) {
        if (graphLinkPool[idx-1].id1 != graphLinkPool[idx].id1 || graphLinkPool[idx-1].id2 != graphLinkPool[idx].id2) {
            indexes.back().numLink = std::max(indexes.back().numLinks[0], indexes.back().numLinks[1]);
            indexes.emplace_back(idx);
        }
        ++indexes.back().numLinks[graphLinkPool[idx].h1 == graphLinkPool[idx].h2];
    }
    indexes.back().numLink = std::max(indexes.back().numLinks[0], indexes.back().numLinks[1]);
    sort(indexes.begin(), indexes.end(), ConsensusLinkPoolIndexGreater());

    es.clear();
    for (unsigned long idxIndex = 0; idxIndex < indexes.size(); ++idxIndex) {
        HiCConsensusEdge HiCEdge;
        HiCEdge.numLink = indexes[idxIndex].numLink;
        HiCEdge.id1 = graphLinkPool[indexes[idxIndex].index].id1;
        HiCEdge.id2 = graphLinkPool[indexes[idxIndex].index].id2;
        HiCEdge.h = indexes[idxIndex].numLinks[0] > indexes[idxIndex].numLinks[1];
		//std::cerr << "HiCEdge: " << HiCEdge.id1 << " " << HiCEdge.id2 << " " << HiCEdge.numLink << " " << HiCEdge.h << std::endl;
//        if (!es.empty() && es.back().id1 == HiCEdge.id2 && es.back().id2 == HiCEdge.id1)
//            continue;
        es.push_back(HiCEdge);
    }

}
//

//added by ouchi
void PairedDBG::getConsensusPart(long id, char h, char mode, vector<ConsensusPart>& scaffold, long &numHiCNode, vector<HiCNode> &hicnode, vector<platanus::Position> &consensusPositionInHiCNode)
{
    long hicnodeID = consensusPositionInHiCNode[id2Index(id)].id;
    long hicOffset = consensusPositionInHiCNode[id2Index(id)].offset;
    const HiCNode& tmpHiCNode = hicnode[id2Index(hicnodeID)];
    char tmph = (h == tmpHiCNode.scaffold[hicOffset].h);

    scaffold.clear();

    char direction = sign(id) * sign(hicnodeID);
    if (mode == 0) {
        if (direction > 0) {
            long start = tmpHiCNode.scaffold[hicOffset].start;
            for (long i = hicOffset; i < tmpHiCNode.numScaffold; ++i) {
                ConsensusPart tmpConsensusPart(tmpHiCNode.scaffold[i].id, tmpHiCNode.scaffold[i].start - start, tmpHiCNode.scaffold[i].end - start, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                scaffold.push_back(tmpConsensusPart);
            }
        } else {
            long end = tmpHiCNode.scaffold[hicOffset].end;
            for (long i = hicOffset; i >= 0; --i) {
                ConsensusPart tmpConsensusPart(-tmpHiCNode.scaffold[i].id, end - tmpHiCNode.scaffold[i].end, end - tmpHiCNode.scaffold[i].start, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                scaffold.push_back(tmpConsensusPart);
            }
        }
    } else if (mode == 1) {
        if (direction > 0) {
            long end = tmpHiCNode.length;
            for (long i = tmpHiCNode.numScaffold - 1; i >= hicOffset; --i) {
                ConsensusPart tmpConsensusPart(-tmpHiCNode.scaffold[i].id, end - tmpHiCNode.scaffold[i].end, end - tmpHiCNode.scaffold[i].start, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                scaffold.push_back(tmpConsensusPart);
            }
        } else {
            for (long i = 0; i <= hicOffset; ++i) {
                ConsensusPart tmpConsensusPart(tmpHiCNode.scaffold[i].id, tmpHiCNode.scaffold[i].start, tmpHiCNode.scaffold[i].end, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                scaffold.push_back(tmpConsensusPart);
            }
        }
    } else if (mode == 2) {
        if (direction > 0) {
            if (hicOffset - 1 >= 0) {
                long end = tmpHiCNode.scaffold[hicOffset-1].end;
                for (long i = hicOffset - 1; i >= 0; --i) {
                    ConsensusPart tmpConsensusPart(-tmpHiCNode.scaffold[i].id, end - tmpHiCNode.scaffold[i].end, end - tmpHiCNode.scaffold[i].start, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                    scaffold.push_back(tmpConsensusPart);
                }
            }
        } else {
            if (hicOffset + 1 < tmpHiCNode.numScaffold) {
                long start = tmpHiCNode.scaffold[hicOffset+1].start;
                for (long i = hicOffset+1; i < tmpHiCNode.numScaffold; ++i) {
                    ConsensusPart tmpConsensusPart(tmpHiCNode.scaffold[i].id, tmpHiCNode.scaffold[i].start - start, tmpHiCNode.scaffold[i].end - start, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                    scaffold.push_back(tmpConsensusPart);
                }
            }
        }
    } else if (mode == 3) {
        if (direction > 0) {
            for (long i = 0; i < hicOffset; ++i) {
                ConsensusPart tmpConsensusPart(tmpHiCNode.scaffold[i].id, tmpHiCNode.scaffold[i].start, tmpHiCNode.scaffold[i].end, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                scaffold.push_back(tmpConsensusPart);
            }
        } else {
            long end = tmpHiCNode.length;
            for (long i = tmpHiCNode.numScaffold - 1; i > hicOffset; --i) {
                ConsensusPart tmpConsensusPart(-tmpHiCNode.scaffold[i].id, end - tmpHiCNode.scaffold[i].end, end - tmpHiCNode.scaffold[i].start, (tmpHiCNode.scaffold[i].h + tmph) % 2);
                scaffold.push_back(tmpConsensusPart);
            }
        }
    }

}
//

//added by ouchi
void PairedDBG::uniteHiCNode(long id1, long id2, char h, long &numHiCNode, vector<HiCNode> &hicnode, vector<platanus::Position> &consensusPositionInHiCNode)
{
    vector<ConsensusPart> scaffolda, scaffoldb;
    getConsensusPart(id1, 0, 1, scaffolda, numHiCNode, hicnode, consensusPositionInHiCNode);
    getConsensusPart(id2, h, 0, scaffoldb, numHiCNode, hicnode, consensusPositionInHiCNode);

	std::cerr << scaffolda[0].id << std::endl;
	std::cerr << scaffoldb[0].id << std::endl;

    long gap = 500;
    HiCNode newHiCNode;
    newHiCNode.numScaffold = scaffolda.size() + scaffoldb.size();

	std::cerr << "numScaffold:" << newHiCNode.numScaffold << std::endl;

    newHiCNode.scaffold = scaffolda;
    for (long i = 0; i < scaffoldb.size(); ++i) {
        ConsensusPart tmpConsensusPart(scaffoldb[i].id, scaffolda.back().end + gap + scaffoldb[i].start, scaffolda.back().end + gap + scaffoldb[i].end, scaffoldb[i].h);
        newHiCNode.scaffold.push_back(tmpConsensusPart);
    }
    newHiCNode.length = newHiCNode.scaffold.back().end;

	std::cerr << newHiCNode.length << std::endl;

    hicnode[id2Index(consensusPositionInHiCNode[id2Index(id1)].id)].state |= SC_DEL;
    hicnode[id2Index(consensusPositionInHiCNode[id2Index(id2)].id)].state |= SC_DEL;
    ++numHiCNode;
    for (long i = 0; i < newHiCNode.numScaffold; ++i) {
        long id = newHiCNode.scaffold[i].id;
        consensusPositionInHiCNode[id2Index(id)].id = sign(id) * numHiCNode;
        consensusPositionInHiCNode[id2Index(id)].offset = i;
    }
    hicnode.push_back(newHiCNode);

	std::cerr << "finish unite" << std::endl;

}
//

//added by ouchi
void PairedDBG::phasing(const long numThread) {
    std::cerr << "phasing" << std::endl;

    const double linkRateThreshold = 0.25;

    makeMPLink(numThread, 0, (*allLibraryMT).size());
    makeHiCLink(numThread);

	long totalNumHetero = 0;
	long totalNumMerge = 0;

    for (long consensusIndex = 0; consensusIndex < numConsensusNode; ++consensusIndex) {
		std::cerr << "consensusNode:" << consensusIndex + 1 << std::endl;
        const ConsensusNode& tmpConsensus = consensusnode[consensusIndex];

		std::cerr << "befor phasing" << std::endl;
		for (long i = 0; i < tmpConsensus.numNode[0]; ++i) {
			std::cerr << "NODE:" << tmpConsensus.node[0][i].id << ", s:" << tmpConsensus.node[0][i].start << ", e:" << tmpConsensus.node[0][i].end << ", ";
		}
		std::cerr << std::endl;
		for (long i = 0; i < tmpConsensus.numNode[1]; ++i) {
			std::cerr << "NODE:" << tmpConsensus.node[1][i].id << ", s:" << tmpConsensus.node[1][i].start << ", e:" << tmpConsensus.node[1][i].end << ", ";
		}
		std::cerr << std::endl;

        vector<std::array<vector<ScaffoldPart>, 2> > blocks;
        blocks.resize(tmpConsensus.numNode[0]);
        for (long i = 0; i < tmpConsensus.numNode[0]; ++i) {
            blocks[i][0].push_back(tmpConsensus.node[0][i]);
        }
        long index = 0;
        for (long i = 0; i < tmpConsensus.numNode[1]; ++i) {
            long maxOverlap = 0;
            long j;
            for (j = index; j < blocks.size(); ++j) {
                if (blocks[j][0].empty()) continue;
                if (tmpConsensus.node[1][i].id == blocks[j][0][0].id) {
                    index = j;
                    maxOverlap = blocks[j][0][0].end - blocks[j][0][0].start;
                    break;
                }
                long overlap = std::min(blocks[j][0][0].end, tmpConsensus.node[1][i].end) - std::max(blocks[j][0][0].start, tmpConsensus.node[1][i].start);
                if (maxOverlap < overlap) {
                    maxOverlap = overlap;
                    index = j;
                }
                if (tmpConsensus.node[1][i].end <= blocks[j][0][0].start) {
                    if (maxOverlap == 0)
                        index = j;
                    break;
                }
            }
            if (maxOverlap == 0) {
                if (j == blocks.size()) index = j;
                std::array<vector<ScaffoldPart>, 2> tmp {{ {}, {} }};
                blocks.insert(blocks.begin()+index, tmp);
                blocks[index][1].push_back(tmpConsensus.node[1][i]);
            } else {
                blocks[index][1].push_back(tmpConsensus.node[1][i]);
            }
        }

        vector<vector<long> > heteroBlocks;
        vector<char> shortFlag;
        vector<long> blockIndex(blocks.size());

        long numHeteroBlock = 0;
        for (long i = 0; i < blocks.size(); ++i) {
            if (blocks[i][0] == blocks[i][1])
                continue;
            ++numHeteroBlock;
            heteroBlocks.resize(numHeteroBlock);
            heteroBlocks[numHeteroBlock-1].push_back(i);
            shortFlag.push_back(0);
            blockIndex[i] = numHeteroBlock-1;
        }
		std::cerr << "numHeteroBlock: " << numHeteroBlock << std::endl;
		long numMerge = 0;

/*
		//distance method
        while(1) {
            vector<std::array<long, 3> > distances;
            long l = 0;
            for (long i = 0; i+1 < heteroBlocks.size(); ++i) {
                if (heteroBlocks[i].empty() || shortFlag[i] >= 2)
                    continue;
                long j = i + 1;
                for (; j < heteroBlocks.size(); ++j) {
                    if (!heteroBlocks[j].empty() && shortFlag[j] < 2)
                        break;
                }
                if (j == heteroBlocks.size()) break;
                long start, end;
                if (blocks[heteroBlocks[j].front()][0].empty()) {
                    if (blocks[heteroBlocks[j].front()][1].empty()) {
                        continue;
                    } else {
                        start = blocks[heteroBlocks[j].front()][1].front().start;
                    }
                } else {
                    if (blocks[heteroBlocks[j].front()][1].empty()) {
                        start = blocks[heteroBlocks[j].front()][0].front().start;
                    } else {
                        start = std::min(blocks[heteroBlocks[j].front()][0].front().start, blocks[heteroBlocks[j].front()][1].front().start);
                    }
                }
                if (blocks[heteroBlocks[i].back()][0].empty()) {
                    if (blocks[heteroBlocks[i].back()][1].empty()) {
                        continue;
                    } else {
                        end = blocks[heteroBlocks[i].back()][1].back().end;
                    }
                } else {
                    if (blocks[heteroBlocks[i].back()][1].empty()) {
                        end = blocks[heteroBlocks[i].back()][0].back().end;
                    } else {
                        end = std::max(blocks[heteroBlocks[i].back()][0].back().end, blocks[heteroBlocks[i].back()][1].back().end);
                    }
                }
                distances.resize(l+1);
                distances[l][0] = heteroBlocks[i][0];
                distances[l][1] = heteroBlocks[j][0];
                distances[l][2] = start - end;
                ++l;
            }
            if (distances.empty()) break;
            std::sort(distances.begin(), distances.end(), [](std::array<long, 3> const& left, std::array<long, 3> const& right){return left[2] < right[2];});
            for (long i = 0; i < distances.size(); ++i) {
                long index1 = blockIndex[distances[i][0]];
                long index2 = blockIndex[distances[i][1]];
				std::cerr << "heteloBlock: " << index1 << "\t" << index2 << "\t" << distances[i][2] << std::endl;
                std::array<std::array<vector<ScaffoldPart>, 2>, 2> scaffolds;
                for (long j = 0; j < heteroBlocks[index1].size(); ++j) {
					std::cerr << "scaffolds[0][0]\t";
                    for (long k = 0; k < blocks[heteroBlocks[index1][j]][0].size(); ++k) {
						std::cerr << "nodeID:" << blocks[heteroBlocks[index1][j]][0][k].id << ", s:" << blocks[heteroBlocks[index1][j]][0][k].start << ", e:" << blocks[heteroBlocks[index1][j]][0][k].end << ", ";
                        scaffolds[0][0].push_back(blocks[heteroBlocks[index1][j]][0][k]);
                    }
					std::cerr << std::endl;
					std::cerr << "scaffolds[0][1]\t";
                    for (long k = 0; k < blocks[heteroBlocks[index1][j]][1].size(); ++k) {
						std::cerr << "nodeID:" << blocks[heteroBlocks[index1][j]][1][k].id << ", s:" << blocks[heteroBlocks[index1][j]][1][k].start << ", e:" << blocks[heteroBlocks[index1][j]][1][k].end << ", ";
                        scaffolds[0][1].push_back(blocks[heteroBlocks[index1][j]][1][k]);
                    }
					std::cerr << std::endl;
                }
                for (long j = 0; j < heteroBlocks[index2].size(); ++j) {
					std::cerr << "scaffolds[1][0]\t";
                    for (long k = 0; k < blocks[heteroBlocks[index2][j]][0].size(); ++k) {
						std::cerr << "nodeID:" << blocks[heteroBlocks[index2][j]][0][k].id << ", s:" << blocks[heteroBlocks[index2][j]][0][k].start << ", e:" << blocks[heteroBlocks[index2][j]][0][k].end << ", ";
                        scaffolds[1][0].push_back(blocks[heteroBlocks[index2][j]][0][k]);
                    }
					std::cerr << std::endl;
					std::cerr << "scaffolds[1][1]\t";
                    for (long k = 0; k < blocks[heteroBlocks[index2][j]][1].size(); ++k) {
						std::cerr << "nodeID:" << blocks[heteroBlocks[index2][j]][1][k].id << ", s:" << blocks[heteroBlocks[index2][j]][1][k].start << ", e:" << blocks[heteroBlocks[index2][j]][1][k].end << ", ";
                        scaffolds[1][1].push_back(blocks[heteroBlocks[index2][j]][1][k]);
                    }
					std::cerr << std::endl;
                }
                //check MatePairLink
                char MPcrossFlag = -1;
                std::array<long, 2> sumLinkForHaplotypeforMP;
                sumLinkForHaplotypeforMP.fill(0);
                getHeteroMPLinkScoreBetweenScaffold(scaffolds, sumLinkForHaplotypeforMP);
				std::cerr << "MPLink: " << sumLinkForHaplotypeforMP[0] << "\t" << sumLinkForHaplotypeforMP[1] << std::endl;
                if (linkRateThreshold * sumLinkForHaplotypeforMP[0] >= sumLinkForHaplotypeforMP[1])
                    MPcrossFlag = 1;
                else if (linkRateThreshold * sumLinkForHaplotypeforMP[1] >= sumLinkForHaplotypeforMP[0])
                    MPcrossFlag = 0;
                if (std::max(sumLinkForHaplotypeforMP[0], sumLinkForHaplotypeforMP[1]) < 3)
                    MPcrossFlag = -1;
                //check HiCLink
                char HiCcrossFlag = -1;
                std::array<long, 2> sumLinkForHaplotypeforHiC;
                sumLinkForHaplotypeforHiC.fill(0);
                getHeteroHiCLinkScoreBetweenScaffold(scaffolds, sumLinkForHaplotypeforHiC);
				std::cerr << "HiCLink: " << sumLinkForHaplotypeforHiC[0] << "\t" << sumLinkForHaplotypeforHiC[1] << std::endl;
                if (linkRateThreshold * sumLinkForHaplotypeforHiC[0] >= sumLinkForHaplotypeforHiC[1])
                    HiCcrossFlag = 1;
                else if (linkRateThreshold * sumLinkForHaplotypeforHiC[1] >= sumLinkForHaplotypeforHiC[0])
                    HiCcrossFlag = 0;
               if (std::max(sumLinkForHaplotypeforHiC[0], sumLinkForHaplotypeforHiC[1]) < 3)
                    HiCcrossFlag = -1;
                char crossFlag;
                if (MPcrossFlag == -1) {
                    if (HiCcrossFlag == -1) {
                        shortFlag[index1] += 1;
                        shortFlag[index2] += 1;
						std::cerr << "result: false" << std::endl;
                        continue;
                    } else {
                        crossFlag = HiCcrossFlag;
                    }
                } else {
                    if (HiCcrossFlag == -1) {
                        crossFlag = MPcrossFlag;
                    } else {
                        if (std::max(sumLinkForHaplotypeforMP[0], sumLinkForHaplotypeforMP[1]) > std::max(sumLinkForHaplotypeforHiC[0], sumLinkForHaplotypeforHiC[1]))
                            crossFlag = MPcrossFlag;
                        else
                            crossFlag = HiCcrossFlag;
                    }
                }
                if (crossFlag == 1) {
                    for (long j = 0; j < heteroBlocks[index2].size(); ++j) {
                        std::swap(blocks[heteroBlocks[index2][j]][0], blocks[heteroBlocks[index2][j]][1]);
                    }
					std::cerr << "result: cross" << std::endl;
                } else {
					std::cerr << "result: parallel" << std::endl;
                }
                //merge heteroBlock
				std::cerr << "merge heteroBlock" << std::endl;
				++numMerge;
                --numHeteroBlock;
                heteroBlocks[index1].insert(heteroBlocks[index1].end(), heteroBlocks[index2].begin(), heteroBlocks[index2].end());
                for (long j = 0; j < heteroBlocks[index2].size(); ++j) {
                    blockIndex[heteroBlocks[index2][j]] = index1;
                }
                heteroBlocks[index2].clear();
                shortFlag[index1] = 0;
                shortFlag[index2] = 0;
            }
        }
*/

		//Link difference method
        vector<char> result(blocks.size(), 0);
        vector<vector<std::array<long, 2> > > HiCLinkScore;
        getHeteroHiCLinkScoreBetweenBlock(blocks, heteroBlocks, HiCLinkScore);
		std::cerr << "HiCLinkScore:" << std::endl;
		for (long i = 0; i < HiCLinkScore.size(); ++i) {
			for (long j = 0; j < HiCLinkScore[i].size(); ++j) {
				std::cerr << HiCLinkScore[i][j][0] << "," << HiCLinkScore[i][j][1] << "\t";
				}
				std::cerr << std::endl;
		}
        while(1) {
            long maxDifference = 0;
            std::array<long, 4> maxLinkPair;
            for (long Index1 = 0; Index1 < numHeteroBlock; ++Index1) {
                for (long Index2 = Index1 + 1; Index2 < numHeteroBlock; ++Index2) {
                    long tmp1 = 0; long tmp2 = 0;
                    for (long i = 0; i < heteroBlocks[Index1].size(); ++i) {
                        long index1 = blockIndex[heteroBlocks[Index1][i]];
                        char h1 = result[heteroBlocks[Index1][i]];
                        for (long j = 0; j < heteroBlocks[Index2].size(); ++j) {
                            long index2 = blockIndex[heteroBlocks[Index2][j]];
                            char h2 = result[heteroBlocks[Index2][j]];
                            tmp1 += HiCLinkScore[index1][index2][h1 != h2];
                            tmp2 += HiCLinkScore[index1][index2][h1 == h2];
                        }
                    }
                    if (maxDifference < abs(tmp1 - tmp2)) {
                        maxDifference = abs(tmp1 - tmp2);
                        maxLinkPair[0] = Index1;
                        maxLinkPair[1] = Index2;
                        maxLinkPair[2] = tmp1;
                        maxLinkPair[3] = tmp2;
                    }
                }
            }
            if (maxDifference == 0) break;
            if (maxLinkPair[2] > maxLinkPair[3]) {
				std::cerr << "result: cross" << std::endl;
                for (long i = 0; i < heteroBlocks[maxLinkPair[1]].size(); ++i) {
                    result[heteroBlocks[maxLinkPair[1]][i]] += 1;
                    result[heteroBlocks[maxLinkPair[1]][i]] = result[heteroBlocks[maxLinkPair[1]][i]] % 2;
                }
            } else if (maxLinkPair[2] < maxLinkPair[3]) {
				std::cerr << "result: parallel" << std::endl;
            } else {
                break;
            }
            ++numMerge;
            --numHeteroBlock;
            heteroBlocks[maxLinkPair[0]].insert(heteroBlocks[maxLinkPair[0]].end(), heteroBlocks[maxLinkPair[1]].begin(), heteroBlocks[maxLinkPair[1]].end());
            heteroBlocks[maxLinkPair[1]].clear();
        }
        for (long i = 0; i < result.size(); ++i) {
            if (result[i] == 1) {
                std::swap(blocks[i][0], blocks[i][1]);
            }
        }


		std::cerr << "numMerge: " << numMerge << std::endl;
		totalNumHetero += numHeteroBlock;
		totalNumMerge += numMerge;

		std::cerr << "after phasing" << std::endl;
		for (long i = 0; i < blocks.size(); ++i) {
			for (long j = 0; j < blocks[i][0].size(); ++j) {
				std::cerr << "NODE:" << blocks[i][0][j].id << ", s:" <<blocks[i][0][j].start << ", e:" << blocks[i][0][j].end << ", ";
			}
			std::cerr << "|";
		}
		std::cerr << std::endl;
		for (long i = 0; i < blocks.size(); ++i) {
			for (long j = 0; j < blocks[i][1].size(); ++j) {
				std::cerr << "NODE:" << blocks[i][1][j].id << ", s:" <<blocks[i][1][j].start << ", e:" << blocks[i][1][j].end << ", ";
			}
			std::cerr << "|";
		}
		std::cerr << std::endl;
    }

	std::cerr << "totalNumHetero: " << totalNumHetero << "\ttotalNumMerge" << totalNumMerge << std::endl;

}
//

//added by ouchi
void PairedDBG::getHeteroMPLinkScoreBetweenScaffold(const std::array<std::array<vector<ScaffoldPart>, 2>, 2> &scaffolds, std::array<long, 2> &sumLinkForHaplotype)
{
    std::array<std::array<GraphNode, 2>, 2> tmpNodes;
    std::array<std::array<GraphEdge, 2>, 2> tmpEdges;
    for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
        for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
            if (scaffolds[leftIndex][rightIndex].empty())
                continue;
            tmpNodes[leftIndex][rightIndex].numContig = 0;
            for (long i = 0; i < scaffolds[leftIndex][rightIndex].size(); ++i) {
                tmpNodes[leftIndex][rightIndex].numContig += node[id2Index(scaffolds[leftIndex][rightIndex][i].id)].numContig;
                long tmpNodeID = id2Index(scaffolds[leftIndex][rightIndex][i].id);
                if (scaffolds[leftIndex][rightIndex][i].id > 0) {
                    for (long j = 0; j < node[tmpNodeID].numContig; ++j) {
                        ScaffoldPart tmpScaffoldPart(node[tmpNodeID].contig[j].id, scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].contig[j].start, scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].contig[j].end);
                        tmpNodes[leftIndex][rightIndex].contig.push_back(tmpScaffoldPart);
                    }
                } else {
                    for (long j = node[tmpNodeID].numContig - 1; j >= 0; --j) {
                        ScaffoldPart tmpScaffoldPart(-(node[tmpNodeID].contig[j].id), scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].length - node[tmpNodeID].contig[j].end, scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].length - node[tmpNodeID].contig[j].start);
                        tmpNodes[leftIndex][rightIndex].contig.push_back(tmpScaffoldPart);
                    }
                }
            }
            tmpEdges[leftIndex][rightIndex].direction = 1;
            tmpEdges[leftIndex][rightIndex].end = 1; //
            tmpEdges[leftIndex][rightIndex].length = tmpNodes[leftIndex][rightIndex].contig.front().start;

            long start = tmpNodes[leftIndex][rightIndex].contig.front().start;
            long end = tmpNodes[leftIndex][rightIndex].contig.back().end;
            for (long i = 0; i < tmpNodes[leftIndex][rightIndex].numContig; ++i) {
                if (start > tmpNodes[leftIndex][rightIndex].contig[i].start)
                    start = tmpNodes[leftIndex][rightIndex].contig[i].start;
                if (end < tmpNodes[leftIndex][rightIndex].contig[i].end)
                    end = tmpNodes[leftIndex][rightIndex].contig[i].end;
            }
            if (start != 0) {
                for (long i = 0; i < tmpNodes[leftIndex][rightIndex].numContig; ++i) {
                    tmpNodes[leftIndex][rightIndex].contig[i].start -= start;
                    tmpNodes[leftIndex][rightIndex].contig[i].end -= start;
                }
            }
            tmpNodes[leftIndex][rightIndex].length = end - start;
        }
    }

    std::array<std::array<vector<char>, 2>, 2> conflictpoints;
    if (tmpNodes[0][0].length != 0 && tmpNodes[0][1].length != 0) {
        getConflictRegionBetweenNodePair(tmpEdges[0][0], tmpEdges[0][1], tmpNodes[0][0], tmpNodes[0][1], conflictpoints[0][0], conflictpoints[0][1]);
    } else {
        conflictpoints[0][0].resize(tmpNodes[0][0].length, 1);
        conflictpoints[0][1].resize(tmpNodes[0][1].length, 1);
    }
    if (tmpNodes[1][0].length != 0 && tmpNodes[1][1].length != 0) {
        getConflictRegionBetweenNodePair(tmpEdges[1][0], tmpEdges[1][1], tmpNodes[1][0], tmpNodes[1][1], conflictpoints[1][0], conflictpoints[1][1]);
    } else {
        conflictpoints[1][0].resize(tmpNodes[1][0].length, 1);
        conflictpoints[1][1].resize(tmpNodes[1][1].length, 1);
    }

	std::cerr << "distance:\t"
			  << std::accumulate(conflictpoints[0][0].begin(), conflictpoints[0][0].end(), 0.0) << "\t"
			  << std::accumulate(conflictpoints[0][1].begin(), conflictpoints[0][1].end(), 0.0) << "\t"
			  << std::accumulate(conflictpoints[1][0].begin(), conflictpoints[1][0].end(), 0.0) << "\t"
			  << std::accumulate(conflictpoints[1][1].begin(), conflictpoints[1][1].end(), 0.0) << std::endl;

    if ((tmpNodes[0][0].length == 0 || tmpNodes[0][1].length == 0) && (tmpNodes[1][0].length == 0 || tmpNodes[1][1].length == 0))
        return;

    long offset1, offset2;
    for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
        for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
            for (long i = 0; i < scaffolds[0][leftIndex].size(); ++i) {
                const GraphNode &tmpNode = node[id2Index(scaffolds[0][leftIndex][i].id)];
                for (long libraryID = 0; libraryID < tmpNode.MPLinks.size(); ++libraryID) {
                    for (long MPLinkID = 0; MPLinkID < tmpNode.MPLinks[libraryID].size(); ++MPLinkID) {
                        const MPLink &mplink = tmpNode.MPLinks[libraryID][MPLinkID];
                        if (mplink.direction * sign(scaffolds[0][leftIndex][i].id) > 0) {
                            long j = 0;
                            for (; j < scaffolds[1][rightIndex].size(); ++j) {
                                if (sign(scaffolds[0][leftIndex][i].id) * mplink.end == scaffolds[1][rightIndex][j].id)
                                    break;
                            }
                            if (j == scaffolds[1][rightIndex].size()) continue;
                            const GraphNode &tmpNode2 = node[id2Index(scaffolds[1][rightIndex][j].id)];
                            offset1 = mplink.offset1;
                            offset2 = mplink.offset2;
                            if (scaffolds[0][leftIndex][i].id < 0)
                                offset1 = tmpNode.length - offset1 -1;
                            if (scaffolds[1][rightIndex][j].id < 0)
                                offset2 = tmpNode2.length - offset2 -1;
                            offset1 += scaffolds[0][leftIndex][i].start - scaffolds[0][leftIndex][0].start;
                            offset2 += scaffolds[1][rightIndex][j].start - scaffolds[1][rightIndex][0].start;
                            if (offset1 < 0 || offset1 > conflictpoints[0][leftIndex].size()-1)
                                continue;
                            if (offset2 < 0 || offset2 > conflictpoints[1][rightIndex].size()-1)
                                continue;
                            if (conflictpoints[0][leftIndex][offset1] == 0 || conflictpoints[1][rightIndex][offset2] == 0)
                                continue;
                            if (scaffolds[1][rightIndex][0].start + offset2 - scaffolds[0][leftIndex][0].start - offset1 <= (*allLibraryMT)[libraryID][0].getAverageInsSize() + tolerence)
                            // || scaffolds[1][rightIndex][0].start + offset2 - scaffolds[0][leftIndex][0].start - offset1 >= (*allLibraryMT)[libraryID][0].getAverageInsSize() - tolerence)
                                sumLinkForHaplotype[leftIndex == rightIndex] += 1;
                        }
                    }
                }
            }
        }
    }
}
//

//added by ouchi
void PairedDBG::getHeteroHiCLinkScoreBetweenScaffold(const std::array<std::array<vector<ScaffoldPart>, 2>, 2> &scaffolds, std::array<long, 2> &sumLinkForHaplotype)
{
    std::array<std::array<GraphNode, 2>, 2> tmpNodes;
    std::array<std::array<GraphEdge, 2>, 2> tmpEdges;
    for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
        for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
            if (scaffolds[leftIndex][rightIndex].empty())
                continue;
            tmpNodes[leftIndex][rightIndex].numContig = 0;
            for (long i = 0; i < scaffolds[leftIndex][rightIndex].size(); ++i) {
                tmpNodes[leftIndex][rightIndex].numContig += node[id2Index(scaffolds[leftIndex][rightIndex][i].id)].numContig;
                long tmpNodeID = id2Index(scaffolds[leftIndex][rightIndex][i].id);
                if (scaffolds[leftIndex][rightIndex][i].id > 0) {
                    for (long j = 0; j < node[tmpNodeID].numContig; ++j) {
                        ScaffoldPart tmpScaffoldPart(node[tmpNodeID].contig[j].id, scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].contig[j].start, scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].contig[j].end);
                        tmpNodes[leftIndex][rightIndex].contig.push_back(tmpScaffoldPart);
                    }
                } else {
                    for (long j = node[tmpNodeID].numContig - 1; j >= 0; --j) {
                        ScaffoldPart tmpScaffoldPart(-(node[tmpNodeID].contig[j].id), scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].length - node[tmpNodeID].contig[j].end, scaffolds[leftIndex][rightIndex][i].start + node[tmpNodeID].length - node[tmpNodeID].contig[j].start);
                        tmpNodes[leftIndex][rightIndex].contig.push_back(tmpScaffoldPart);
                    }
                }
            }
            tmpEdges[leftIndex][rightIndex].direction = 1;
            tmpEdges[leftIndex][rightIndex].end = 1; //
            tmpEdges[leftIndex][rightIndex].length = tmpNodes[leftIndex][rightIndex].contig.front().start;

            long start = tmpNodes[leftIndex][rightIndex].contig.front().start;
            long end = tmpNodes[leftIndex][rightIndex].contig.back().end;
            for (long i = 0; i < tmpNodes[leftIndex][rightIndex].numContig; ++i) {
                if (start > tmpNodes[leftIndex][rightIndex].contig[i].start)
                    start = tmpNodes[leftIndex][rightIndex].contig[i].start;
                if (end < tmpNodes[leftIndex][rightIndex].contig[i].end)
                    end = tmpNodes[leftIndex][rightIndex].contig[i].end;
            }
            if (start != 0) {
                for (long i = 0; i < tmpNodes[leftIndex][rightIndex].numContig; ++i) {
                    tmpNodes[leftIndex][rightIndex].contig[i].start -= start;
                    tmpNodes[leftIndex][rightIndex].contig[i].end -= start;
                }
            }
            tmpNodes[leftIndex][rightIndex].length = end - start;
        }
    }

    std::array<std::array<vector<char>, 2>, 2> conflictpoints;
    if (tmpNodes[0][0].length != 0 && tmpNodes[0][1].length != 0) {
        getConflictRegionBetweenNodePair(tmpEdges[0][0], tmpEdges[0][1], tmpNodes[0][0], tmpNodes[0][1], conflictpoints[0][0], conflictpoints[0][1]);
    } else {
        conflictpoints[0][0].resize(tmpNodes[0][0].length, 1);
        conflictpoints[0][1].resize(tmpNodes[0][1].length, 1);
    }
    if (tmpNodes[1][0].length != 0 && tmpNodes[1][1].length != 0) {
        getConflictRegionBetweenNodePair(tmpEdges[1][0], tmpEdges[1][1], tmpNodes[1][0], tmpNodes[1][1], conflictpoints[1][0], conflictpoints[1][1]);
    } else {
        conflictpoints[1][0].resize(tmpNodes[1][0].length, 1);
        conflictpoints[1][1].resize(tmpNodes[1][1].length, 1);
    }
    if ((tmpNodes[0][0].length == 0 || tmpNodes[0][1].length == 0) && (tmpNodes[1][0].length == 0 || tmpNodes[1][1].length == 0))
        return;

    long offset1, offset2;
    for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
        for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
            for (long i = 0; i < scaffolds[0][leftIndex].size(); ++i) {
                const GraphNode &tmpNode = node[id2Index(scaffolds[0][leftIndex][i].id)];
                for (long HiCLinkID = 0; HiCLinkID < tmpNode.HiCLinks.size(); ++HiCLinkID) {
                    const HiCLink &hiclink = tmpNode.HiCLinks[HiCLinkID];
                    long j = 0;
                    for (; j < scaffolds[1][rightIndex].size(); ++j) {
                        if (abs(hiclink.end) == abs(scaffolds[1][rightIndex][j].id))
                            break;
                    }
                    if (j == scaffolds[1][rightIndex].size()) continue;
                    const GraphNode &tmpNode2 = node[id2Index(scaffolds[1][rightIndex][j].id)];
                    offset1 = hiclink.offset1;
                    offset2 = hiclink.offset2;
                    if (scaffolds[0][leftIndex][i].id < 0)
                        offset1 = tmpNode.length - offset1 -1;
                    if (scaffolds[1][rightIndex][j].id < 0)
                        offset2 = tmpNode2.length - offset2 -1;
                    offset1 += scaffolds[0][leftIndex][i].start - scaffolds[0][leftIndex][0].start;
                    offset2 += scaffolds[1][rightIndex][j].start - scaffolds[1][rightIndex][0].start;
                    if (offset1 < 0 || offset1 > conflictpoints[0][leftIndex].size()-1)
                        continue;
                    if (offset2 < 0 || offset2 > conflictpoints[1][rightIndex].size()-1)
                        continue;
                    if (conflictpoints[0][leftIndex][offset1] == 0 || conflictpoints[1][rightIndex][offset2] == 0)
                        continue;
                    sumLinkForHaplotype[leftIndex == rightIndex] += 1;
                }
            }
        }
    }

}
//

//added by ouchi
void PairedDBG::getHeteroHiCLinkScoreBetweenBlock(const vector<std::array<vector<ScaffoldPart>, 2> > &blocks, const vector<vector<long> > heteroBlocks, vector<vector<std::array<long, 2> > > &HiCLinkScore)
{
    long numHetero = heteroBlocks.size();
    HiCLinkScore.clear();
    HiCLinkScore.resize(numHetero);
    for (long i = 0; i < numHetero; ++i) {
        HiCLinkScore[i].resize(numHetero);
    }
    vector<std::array<vector<char>, 2> > conflictpoints;
    conflictpoints.resize(numHetero);
    for (long Index = 0; Index < numHetero; ++Index) {
        long blockIndex = heteroBlocks[Index][0];
        std::array<GraphNode, 2> tmpNodes;
        std::array<GraphEdge, 2> tmpEdges;
        for (char h = 0; h < 2; ++h) {
            if (blocks[blockIndex][h].empty())
                continue;
            tmpNodes[h].numContig = 0;
            for (long i = 0; i < blocks[blockIndex][h].size(); ++i) {
                long tmpNodeID = id2Index(blocks[blockIndex][h][i].id);
                tmpNodes[h].numContig += node[tmpNodeID].numContig;
                if (blocks[blockIndex][h][i].id > 0) {
                    for (long j = 0; j < node[tmpNodeID].numContig; ++j) {
                        ScaffoldPart tmpScaffoldPart(node[tmpNodeID].contig[j].id, blocks[blockIndex][h][i].start + node[tmpNodeID].contig[j].start, blocks[blockIndex][h][i].start + node[tmpNodeID].contig[j].end);
                        tmpNodes[h].contig.push_back(tmpScaffoldPart);
                    }
                } else {
                    for (long j = node[tmpNodeID].numContig - 1; j >= 0; --j) {
                        ScaffoldPart tmpScaffoldPart(-(node[tmpNodeID].contig[j].id), blocks[blockIndex][h][i].start + node[tmpNodeID].length - node[tmpNodeID].contig[j].end, blocks[blockIndex][h][i].start + node[tmpNodeID].length - node[tmpNodeID].contig[j].start);
                        tmpNodes[h].contig.push_back(tmpScaffoldPart);
                    }
                }
            }
            tmpEdges[h].direction = 1;
            tmpEdges[h].end = 1; //
            tmpEdges[h].length = tmpNodes[h].contig.front().start;

            long start = tmpNodes[h].contig.front().start;
            long end = tmpNodes[h].contig.back().end;
            for (long i = 0; i < tmpNodes[h].numContig; ++i) {
                if (start > tmpNodes[h].contig[i].start)
                    start = tmpNodes[h].contig[i].start;
                if (end < tmpNodes[h].contig[i].end)
                    end = tmpNodes[h].contig[i].end;
            }
            if (start != 0) {
                for (long i = 0; i < tmpNodes[h].numContig; ++i) {
                    tmpNodes[h].contig[i].start -= start;
                    tmpNodes[h].contig[i].end -= start;
                }
            }
            tmpNodes[h].length = end - start;
        }

        if (tmpNodes[0].length != 0 && tmpNodes[1].length != 0) {
            getConflictRegionBetweenNodePair(tmpEdges[0], tmpEdges[1], tmpNodes[0], tmpNodes[1], conflictpoints[Index][0], conflictpoints[Index][1]);
        } else {
            conflictpoints[Index][0].resize(tmpNodes[0].length, 1);
            conflictpoints[Index][1].resize(tmpNodes[1].length, 1);
        }
    }

    for (long Index = 0; Index < numHetero; ++Index) {
        long blockIndex = heteroBlocks[Index][0];
        for (long leftIndex = 0; leftIndex < 2; ++leftIndex) {
            for (long i = 0; i < blocks[blockIndex][leftIndex].size(); ++i) {
                const GraphNode &tmpNode = node[id2Index(blocks[blockIndex][leftIndex][i].id)];
                for (long HiCLinkID = 0; HiCLinkID < tmpNode.HiCLinks.size(); ++HiCLinkID) {
                    const HiCLink &hiclink = tmpNode.HiCLinks[HiCLinkID];
                    for (long rightIndex = 0; rightIndex < 2; ++rightIndex) {
                        for (long Index2 = 0; Index2 < numHetero; ++Index2) {
                            long blockIndex2 = heteroBlocks[Index2][0];
                            long j = 0;
                            for (; j < blocks[blockIndex2][rightIndex].size(); ++j) {
                                if (abs(hiclink.end) == abs(blocks[blockIndex2][rightIndex][j].id))
                                    break;
                            }
                            if (j == blocks[blockIndex2][rightIndex].size()) continue;
                            const GraphNode &tmpNode2 = node[id2Index(blocks[blockIndex2][rightIndex][j].id)];
                            long offset1 = hiclink.offset1;
                            long offset2 = hiclink.offset2;
                            if (blocks[blockIndex][leftIndex][i].id < 0)
                                offset1 = tmpNode.length - offset1 - 1;
                            if (blocks[blockIndex2][rightIndex][j].id < 0)
                                offset2 = tmpNode2.length - offset2 - 1;
                            offset1 += blocks[blockIndex][leftIndex][i].start - blocks[blockIndex][leftIndex][0].start;
                            offset2 += blocks[blockIndex2][rightIndex][i].start - blocks[blockIndex2][rightIndex][0].start;
                            if (offset1 < 0 || offset1 > conflictpoints[Index][leftIndex].size()-1)
                                continue;
                            if (offset2 < 0 || offset2 > conflictpoints[Index2][rightIndex].size()-1)
                                continue;
                            if (conflictpoints[Index][leftIndex][offset1] == 0 || conflictpoints[Index2][rightIndex][offset2] == 0)
                                continue;
                            HiCLinkScore[Index][Index2][leftIndex == rightIndex] += 1;
                        }
                    }
                }
            }
        }
    }


}
//


//added by ocuhi
void PairedDBG::consensus2seq(const ConsensusNode &consensusnode, vector<char> &ret)
{
    std::array<vector<char>, 2> tmpSeq;
    vector<char> tmp;
    long k;
    ret.clear();

    for (char h = 0; h < 2; ++h) {
        for (long j = 0; j < consensusnode.numNode[h]; ++j) {
            long nodeID = consensusnode.node[h][j].id;
            node2seq(node[id2Index(nodeID)], tmp);
            if (j == 0) {
                k = 0 - consensusnode.node[h][j].start;
            } else {
                k = consensusnode.node[h][j-1].end - consensusnode.node[h][j].start;
            }
            for (; k < 0; ++k)
                tmpSeq[h].push_back(4);

            if (nodeID > 0) {
                if (k < tmp.size()) {
                    tmpSeq[h].insert(tmpSeq[h].end(), tmp.begin() + k, tmp.end());
                }
            } else {
                for (; k < tmp.size(); ++k) {
                    if (tmp[tmp.size() - k - 1 ] >= 4)
                        tmpSeq[h].push_back(4);
                    else
                        tmpSeq[h].push_back(0x3 ^ tmp[tmp.size() - k - 1]);
                }
            }
        }
    }

    for (long j = 0; j < std::max(tmpSeq[0].size(), tmpSeq[1].size()); ++j) {
        if (j < tmpSeq[0].size()) {
            if (j < tmpSeq[1].size()) {
                if (tmpSeq[0][j] != 4)
                   ret.push_back(tmpSeq[0][j]);
                else if (tmpSeq[1][j] != 4)
                   ret.push_back(tmpSeq[1][j]);
                else
                   ret.push_back(4);
            } else {
                ret.push_back(tmpSeq[0][j]);
            }
        } else {
            if (j < tmpSeq[1].size()) {
                ret.push_back(tmpSeq[1][j]);
            } else {
                ret.push_back(4);
            }
        }
    }
}
//


//added by ouchi
bool PairedDBG::checkHiCLinkBetweenHiCNodePair(vector<ConsensusPart> scaffold1, vector<ConsensusPart> scaffold2)
{
    vector<vector<int> > contactmap;
    long boundary = makeContactmapOfScaffold(scaffold1, scaffold2, contactmap);
    if (boundary == -1)
        return 1;
    deleteGapRegionFromContactmap(contactmap, boundary);
    std::array<double, 6> result;
    calcSeparation(contactmap, 0, boundary, contactmap.size(), result, true);
	std::cerr << "result: " << result[0]/result[1] << "\t" << result[3]/result[4] << "\t" << result[2] << "\t" << result[5] << std::endl;
    return (result[0]/result[1] + result[3]/result[4] < HiC_CONTACTMAP_S_THRESHOLD || result[2] + result[5] < HiC_CONTACTMAP_T_THRESHOLD);
}
//

//added by ouchi
long PairedDBG::makeContactmapOfScaffold(vector<ConsensusPart> scaffold1, vector<ConsensusPart> scaffold2, vector<vector<int> >& contactmap)
{
	std::cerr << "makeContactmapOfScaffold" << std::endl;
    const long binsize = 100000;
    long length1 = 0;
    long length2 = 0;
    for (long i = 0; i < scaffold1.size(); ++i) {
        length1 += consensusnode[id2Index(scaffold1[i].id)].length;
    }
    for (long i = 0; i < scaffold2.size(); ++i) {
        length2 += consensusnode[id2Index(scaffold2[i].id)].length;
    }

    if (length1 / binsize < 1 || length2 / binsize < 1)
        return -1;

    long numbin = length1 / binsize + length2 / binsize;
    long lower, upper, boundary;
    lower = length1 % binsize;
    upper = length1 + length2 - length2 % binsize;
    boundary = length1 / binsize;

    contactmap.resize(numbin);
    for (long i = 0; i < numbin; ++i) {
        contactmap[i].resize(numbin);
    }

    for (char h = 0; h < 2; ++h) {
        long tmpLength = 0;
        for (long consensusIndex = 0; consensusIndex < scaffold1.size(); ++consensusIndex) {
            const ConsensusNode &tmpConsensus = consensusnode[id2Index(scaffold1[consensusIndex].id)];
			# pragma omp parallel for schedule(dynamic)
            for (long nodeID = 0; nodeID < tmpConsensus.numNode[(h+scaffold1[consensusIndex].h)%2]; ++nodeID) {
                const GraphNode &tmpNode = node[id2Index(tmpConsensus.node[(h+scaffold1[consensusIndex].h)%2][nodeID].id)];
                for (long HiCLinkID = 0; HiCLinkID < tmpNode.HiCLinks.size(); ++HiCLinkID) {
                    const HiCLink &hiclink = tmpNode.HiCLinks[HiCLinkID];
                    long offset1, offset2;
                    offset1 = tmpLength;
                    if (scaffold1[consensusIndex].id > 0)
                        offset1 += tmpConsensus.node[(h+scaffold1[consensusIndex].h)%2][nodeID].start;
                    else
                        offset1 += (tmpConsensus.length - tmpConsensus.node[(h+scaffold1[consensusIndex].h)%2][nodeID].end);
                    if (scaffold1[consensusIndex].id * tmpConsensus.node[(h+scaffold1[consensusIndex].h)%2][nodeID].id > 0)
                        offset1 += hiclink.offset1;
                    else
                        offset1 += tmpNode.length - hiclink.offset1;

                    offset2 = -1;
                    if (nodePositionInConsensus[id2Index(hiclink.end)].size() > 0) {
                        long tmpConsensusID = nodePositionInConsensus[id2Index(hiclink.end)][0].second.id;
                        long tmpOffset = nodePositionInConsensus[id2Index(hiclink.end)][0].second.offset;
                        long tmph = nodePositionInConsensus[id2Index(hiclink.end)][0].first;

                        const ConsensusNode& tmpConsensus2 = consensusnode[id2Index(tmpConsensusID)];
                        long tmpLength2 = 0;
                        long i = 0;
                        for (; i < scaffold1.size(); ++i) {
                            if (abs(scaffold1[i].id) == abs(tmpConsensusID))
                                break;
                            tmpLength2 += consensusnode[id2Index(scaffold1[i].id)].length;
                        }
                        if (i < scaffold1.size()) {
                            if (nodePositionInConsensus[id2Index(tmpConsensus.node[(h+scaffold1[consensusIndex].h)%2][nodeID].id)].size() == 2 && nodePositionInConsensus[id2Index(hiclink.end)].size() == 1) {
                                if (h != (tmph+scaffold1[i].h)%2)
                                    continue;
                            }
                            offset2 = tmpLength2;
                            if (scaffold1[i].id > 0)
                                offset2 += tmpConsensus2.node[tmph][tmpOffset].start;
                            else
                                offset2 += tmpConsensus2.length - tmpConsensus2.node[tmph][tmpOffset].end;
                            if (scaffold1[i].id * tmpConsensus2.node[tmph][tmpOffset].id > 0)
                                offset2 += hiclink.offset2;
                            else
                                offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        } else {
                            i = 0;
                            for (; i < scaffold2.size(); ++i) {
                                if (abs(scaffold2[i].id) == abs(tmpConsensusID))
                                    break;
                                tmpLength2 += consensusnode[id2Index(scaffold2[i].id)].length;
                            }
                            if (i < scaffold2.size()) {
                                if (nodePositionInConsensus[id2Index(tmpConsensus.node[(h+scaffold1[consensusIndex].h)%2][nodeID].id)].size() == 2 && nodePositionInConsensus[id2Index(hiclink.end)].size() == 1) {
                                    if (h != (tmph+scaffold2[i].h)%2)
                                        continue;
                                }
                                offset2 = tmpLength2;
                                if (scaffold2[i].id > 0)
                                    offset2 += tmpConsensus2.node[tmph][tmpOffset].start;
                                else
                                    offset2 += tmpConsensus2.length - tmpConsensus2.node[tmph][tmpOffset].end;
                                if (scaffold2[i].id * tmpConsensus2.node[tmph][tmpOffset].id > 0)
                                    offset2 += hiclink.offset2;
                                else
                                    offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                            }
                        }
                    }

					#pragma omp critical (contactmap)
					{
                        if (lower <= offset1 && offset1 < upper && lower <= offset2 && offset2 < upper) {
                            ++contactmap[(offset1-lower)/binsize][(offset2-lower)/binsize];
                        }
					}

                }

            }
            tmpLength += tmpConsensus.length;
        }

        for (long consensusIndex = 0; consensusIndex < scaffold2.size(); ++consensusIndex) {
            const ConsensusNode &tmpConsensus = consensusnode[id2Index(scaffold2[consensusIndex].id)];
			# pragma omp parallel for schedule(dynamic)
            for (long nodeID = 0; nodeID < tmpConsensus.numNode[(h+scaffold2[consensusIndex].h)%2]; ++nodeID) {
                const GraphNode &tmpNode = node[id2Index(tmpConsensus.node[(h+scaffold2[consensusIndex].h)%2][nodeID].id)];
                for (long HiCLinkID = 0; HiCLinkID < tmpNode.HiCLinks.size(); ++HiCLinkID) {
                    const HiCLink &hiclink = tmpNode.HiCLinks[HiCLinkID];
                    long offset1, offset2;
                    offset1 = tmpLength;
                    if (scaffold2[consensusIndex].id > 0)
                        offset1 += tmpConsensus.node[(h+scaffold2[consensusIndex].h)%2][nodeID].start;
                    else
                        offset1 += (tmpConsensus.length - tmpConsensus.node[(h+scaffold2[consensusIndex].h)%2][nodeID].end);
                    if (scaffold2[consensusIndex].id * tmpConsensus.node[(h+scaffold2[consensusIndex].h)%2][nodeID].id > 0)
                        offset1 += hiclink.offset1;
                    else
                        offset1 += tmpNode.length - hiclink.offset1;

                    offset2 = -1;
                    if (nodePositionInConsensus[id2Index(hiclink.end)].size() > 0) {
                        long tmpConsensusID = nodePositionInConsensus[id2Index(hiclink.end)][0].second.id;
                        long tmpOffset = nodePositionInConsensus[id2Index(hiclink.end)][0].second.offset;
                        long tmph = nodePositionInConsensus[id2Index(hiclink.end)][0].first;

                        const ConsensusNode& tmpConsensus2 = consensusnode[id2Index(tmpConsensusID)];
                        long tmpLength2 = 0;
                        long i = 0;
                        for (; i < scaffold1.size(); ++i) {
                            if (abs(scaffold1[i].id) == abs(tmpConsensusID))
                                break;
                            tmpLength2 += consensusnode[id2Index(scaffold1[i].id)].length;
                        }
                        if (i < scaffold1.size()) {
                            if (nodePositionInConsensus[id2Index(tmpConsensus.node[(h+scaffold2[consensusIndex].h)%2][nodeID].id)].size() == 2 && nodePositionInConsensus[id2Index(hiclink.end)].size() == 1) {
                                if (h != (tmph+scaffold1[i].h)%2)
                                    continue;
                            }
                            offset2 = tmpLength2;
                            if (scaffold1[i].id > 0)
                                offset2 += tmpConsensus2.node[tmph][tmpOffset].start;
                            else
                                offset2 += tmpConsensus2.length - tmpConsensus2.node[tmph][tmpOffset].end;
                            if (scaffold1[i].id * tmpConsensus2.node[tmph][tmpOffset].id > 0)
                                offset2 += hiclink.offset2;
                            else
                                offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                        } else {
                            i = 0;
                            for (; i < scaffold2.size(); ++i) {
                                if (abs(scaffold2[i].id) == abs(tmpConsensusID))
                                    break;
                                tmpLength2 += consensusnode[id2Index(scaffold2[i].id)].length;
                            }
                            if (i < scaffold2.size()) {
                                if (nodePositionInConsensus[id2Index(tmpConsensus.node[(h+scaffold2[consensusIndex].h)%2][nodeID].id)].size() == 2 && nodePositionInConsensus[id2Index(hiclink.end)].size() == 1) {
                                    if (h != (tmph+scaffold2[i].h)%2)
                                        continue;
                                }
                                offset2 = tmpLength2;
                                if (scaffold2[i].id > 0)
                                    offset2 += tmpConsensus2.node[tmph][tmpOffset].start;
                                else
                                    offset2 += tmpConsensus2.length - tmpConsensus2.node[tmph][tmpOffset].end;
                                if (scaffold2[i].id * tmpConsensus2.node[tmph][tmpOffset].id > 0)
                                    offset2 += hiclink.offset2;
                                else
                                    offset2 += node[id2Index(hiclink.end)].length - hiclink.offset2;
                            }
                        }
                    }

					#pragma omp critical (contactmap)
					{
                        if (lower <= offset1 && offset1 < upper && lower <= offset2 && offset2 < upper) {
                            ++contactmap[(offset1-lower)/binsize][(offset2-lower)/binsize];
                        }
					}

                }

            }
            tmpLength += tmpConsensus.length;
        }
    }

	std::cerr << "contactmap" << std::endl;
	for (long i = 0; i < contactmap.size(); ++i) {
		if (i == boundary) {
			for (long j = 0; j < contactmap[i].size(); ++j) {
				std::cerr << "_\t";
			}
			std::cerr << std::endl;
		}
		for (long j = 0; j < contactmap[i].size(); ++j) {
			if (j == boundary) std::cerr << "|\t";
			std::cerr << contactmap[i][j] << "\t";
		}
		std::cerr << std::endl;
	}

    return boundary;

}
//

//added by ouchi
void PairedDBG::outputConsensusFastg(const string &outputFilename)
{
    for (long consensusID = 0; consensusID < numConsensusNode; ++consensusID) {
        cerr << consensusID << endl;
        for (char h = 0; h < 2; ++h) {
            for (long i = 0; i < consensusnode[consensusID].numNode[h]; ++i) {
                cerr << "NODE:" << consensusnode[consensusID].node[h][i].id << ", s:" << consensusnode[consensusID].node[h][i].start << ", e:" << consensusnode[consensusID].node[h][i].end << ", ";
            }
            cerr << endl;
        }
    }

    std::ofstream out(outputFilename.c_str());

	vector<string> nodeName;
	for (long consensusID = 0; consensusID < numConsensusNode; ++consensusID) {
        ConsensusNode &consensusnode = this->consensusnode[consensusID];
		std::ostringstream oss;
		oss << "NODE_" << consensusID + 1 << "_"
			<< "length_" << consensusnode.length << "_"
			<< "cov_" << calcConsensusCoverage(consensusnode);

		nodeName.push_back(oss.str());
    }

    for (long consensusID = 0; consensusID < numConsensusNode; ++consensusID) {
        ConsensusNode &consensusnode1 = this->consensusnode[consensusID];
		vector<char> nodeSeq;
		unsigned long i;
		unsigned currentNumEdge;

		unsigned numLeftEdge = 0;
		unsigned numRightEdge = 0;
        for (long edgeID = 0; edgeID < consensusnode1.numEdge; ++edgeID) {
			ConsensusEdge &edge = (consensusnode1.edge[edgeID]);
			if (edge.direction > 0)
				++numRightEdge;
			else
				++numLeftEdge;
		}

		currentNumEdge = 0;
		out << '>' << nodeName[consensusID];
		if (numRightEdge > 0)
			out << ':';
		for (long edgeID = 0; edgeID < consensusnode1.numEdge; ++edgeID) {
			ConsensusEdge &edge = (consensusnode1.edge[edgeID]);
			if (edge.direction < 0)
				continue;

			if (currentNumEdge > 0)
				out << ',';
			out << nodeName[abs(consensusnode1.edge[edgeID].end) - 1];
			if (consensusnode1.edge[edgeID].end < 0)
				out << '\'';
			++currentNumEdge;
		}
		out << ";\n";

		consensus2seq(consensusnode1, nodeSeq);
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
		out << '>' << nodeName[consensusID];
		out << '\'' << ':';
		for (long edgeID = 0; edgeID < consensusnode1.numEdge; ++edgeID) {
			ConsensusEdge &edge = (consensusnode1.edge[edgeID]);
			if (edge.direction > 0)
				continue;

			if (currentNumEdge > 0)
				out << ',';
			out << nodeName[abs(consensusnode1.edge[edgeID].end) - 1];
			if (consensusnode1.edge[edgeID].end > 0)
				out << '\'';
			++currentNumEdge;
		}
		out << ";\n";

		consensus2seq(consensusnode1, nodeSeq);
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

//added by ouchi
void PairedDBG::outputHiCNode(const string &outFilename, const long numHiCNode, const vector<HiCNode> &hicnode)
{
    std::ofstream out(outFilename);
    long numOutputNode = 0;
    for (long hicIndex = 0; hicIndex < numHiCNode; ++hicIndex) {
        const HiCNode &tmpHiCNode = hicnode[hicIndex];
        if (tmpHiCNode.state & SC_DEL) continue;
        ConsensusNode newConsensus;
        for (char h = 0; h < 2; ++h) {
            long numNewNode = 0;
            for (long j = 0; j < tmpHiCNode.numScaffold; ++j) {
                numNewNode += consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].numNode[(tmpHiCNode.scaffold[j].h + h) % 2];
            }
            newConsensus.numNode[h] = numNewNode;
            newConsensus.node[h].resize(newConsensus.numNode[h]);
            long num = 0;
            for (long j = 0; j < tmpHiCNode.numScaffold; ++j) {
                if (tmpHiCNode.scaffold[j].id > 0) {
                    for (long k = 0; k < consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].numNode[(tmpHiCNode.scaffold[j].h + h) % 2]; ++k) {
                        ScaffoldPart tmpScaffoldPart(consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].node[(tmpHiCNode.scaffold[j].h + h) % 2][k].id, tmpHiCNode.scaffold[j].start + consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].node[(tmpHiCNode.scaffold[j].h + h) % 2][k].start, tmpHiCNode.scaffold[j].start + consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].node[(tmpHiCNode.scaffold[j].h + h) % 2][k].end);
                        newConsensus.node[h][num] = tmpScaffoldPart;
                        ++num;
                    }
                } else {
                    for (long k = consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].numNode[(tmpHiCNode.scaffold[j].h + h) % 2] - 1; k >= 0; --k) {
                        ScaffoldPart tmpScaffoldPart(-(consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].node[(tmpHiCNode.scaffold[j].h + h) % 2][k].id), tmpHiCNode.scaffold[j].start + consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].length - consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].node[(tmpHiCNode.scaffold[j].h + h) % 2][k].end, tmpHiCNode.scaffold[j].start + consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].length - consensusnode[id2Index(tmpHiCNode.scaffold[j].id)].node[(tmpHiCNode.scaffold[j].h + h) % 2][k].start);
                        newConsensus.node[h][num] = tmpScaffoldPart;
                        ++num;
                    }
                }
            }
        }

        long min_start, max_end;
        if (newConsensus.numNode[0] > 0) {
            long start = newConsensus.node[0].front().start;
            long end = newConsensus.node[0].back().end;
            for (long j = 0; j < newConsensus.numNode[0]; ++j) {
                if (start > newConsensus.node[0][j].start)
                    start = newConsensus.node[0][j].start;
                if (end < newConsensus.node[0][j].end)
                    end = newConsensus.node[0][j].end;
            }
            min_start = start;
            max_end = end;
        }
        if (newConsensus.numNode[1] > 0) {
            long start2 = newConsensus.node[1].front().start;
            long end2 = newConsensus.node[1].back().end;
            for (long j = 0; j < newConsensus.numNode[1]; ++j) {
                if (start2 > newConsensus.node[1][j].start)
                    start2 = newConsensus.node[1][j].start;
                if (end2 < newConsensus.node[1][j].end)
                    end2 = newConsensus.node[1][j].end;
            }
            if (newConsensus.numNode[0] > 0) {
                min_start = std::min(min_start, start2);
                max_end = std::max(max_end, end2);
            } else {
                min_start = start2;
                max_end = end2;
            }
        }
        if (min_start != 0) {
            for (long j = 0; j < newConsensus.numNode[0]; ++j) {
                newConsensus.node[0][j].start -= min_start;
                newConsensus.node[0][j].end -= min_start;
            }
            for (long j = 0; j < newConsensus.numNode[1]; ++j) {
                newConsensus.node[1][j].start -= min_start;
                newConsensus.node[1][j].end -= min_start;
            }
        }
        newConsensus.length = max_end - min_start;

        vector<char> nodeSeq;
        consensus2seq(newConsensus, nodeSeq);
        ++numOutputNode;


        std::ostringstream oss;
        oss << ">seq" << numOutputNode << "_"
            << "len" << nodeSeq.size() << "_"
            << "cov" << calcConsensusCoverage(newConsensus) << endl;

        out << oss.str();

        unsigned long i;
        for (i = 0; i < nodeSeq.size(); ++i) {
            out << platanus::Bin2Char(nodeSeq[i]);
            if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
        }
        if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }

    out.close();

}
//

//
