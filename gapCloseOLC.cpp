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

#include "gapCloseOLC.h"
#include "mapper.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <list>


using std::string;
using std::vector;


//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
// not change value below parameters!!
const unsigned OverlapLayoutConsensus::HEADID = 0;
const unsigned OverlapLayoutConsensus::TAILID = 1;



//////////////////////////////////////////////////////////////////////////////////////
// make kmer table to approxiate overlap reads
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::makeKmerTable(void)
{
    for (unsigned long long i = 0; i < numSeq; ++i) {
        string kmer;
        if (gapInfo->seq[i].size() < kmerLength) continue;
        for (unsigned j = 0; j < kmerLength; ++j) {
            char tmp = gapInfo->seq[i][j] + static_cast<char>(1);
            kmer += tmp;
        }
        kmerList[kmer].emplace_back(i);
        for (unsigned j = kmerLength, n = gapInfo->seq[i].size(); j < n; ++j) {
            kmer.erase(0, 1);
            kmer += (gapInfo->seq[i][j] + static_cast<char>(1));
            kmerList[kmer].emplace_back(i);
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// Check table is written whether tow reads overlap or not
//////////////////////////////////////////////////////////////////////////////////////
inline void OverlapLayoutConsensus::updateCheckTable(vector<vector<ALIGN_STATUS> > &checkTable, const unsigned long long readID1, const unsigned long long readID2, ALIGN_STATUS align)
{
    checkTable[readID1][readID2 - readID1 - 1] = align;
}




//////////////////////////////////////////////////////////////////////////////////////
// check overlap between reads using pairwise alignment
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::makeOverlapGraph(const unsigned minOverlap, const double maxMissRate)
{
    this->initSeqLink();

    vector<vector<ALIGN_STATUS> > checkTable(numSeq);
    for (unsigned long long i = 0; i < numSeq - 1; ++i)
        checkTable[i].resize(numSeq - 1 - i, YET);

    // check seqs overlap which have same kmer
    for (auto kmerListIterator = kmerList.begin(), kmerListEnd = kmerList.end(); kmerListIterator != kmerListEnd; ++kmerListIterator) {
        for (auto pair1 = kmerListIterator->second.begin(), listEnd = kmerListIterator->second.end(); pair1 != listEnd; ++pair1) {
            auto pair2 = pair1;
            for (++pair2; pair2 != listEnd; ++pair2) {
                if (*pair1 != *pair2 && checkAlreadyAlign(checkTable, *pair1, *pair2)) {
                    const ALIGN_STATUS aligned = pairwiseAlignment(*pair1, *pair2, minOverlap, maxMissRate);
                    updateCheckTable(checkTable, *pair1, *pair2, aligned);
                }
            }
        }
    }

    // check seq overkap ownselef
    for (unsigned long long i = 0; i < numSeq; ++i) {
        pairwiseAlignmentOwnself(i, minOverlap, maxMissRate);
        seqLink[i].tailLink.sort();
    }
    deleteKmerList();
}



//////////////////////////////////////////////////////////////////////////////////////
// set parent seq information in seqLink
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::initSeqLink(void)
{
    unsigned id = 0;
    auto linkIt = seqLink.begin();
    auto linkEnd = seqLink.end();
    for(; linkIt != linkEnd; ++linkIt) {
        linkIt->seqID = id;
        linkIt->seqLength = gapInfo->seq[id].size();
        ++id;
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// exec pairwise alignment (wrapper function)
//////////////////////////////////////////////////////////////////////////////////////
inline OverlapLayoutConsensus::ALIGN_STATUS OverlapLayoutConsensus::pairwiseAlignment(const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const double maxMissRate)
{
    return pairwiseAlignmentOLC(readID1, readID2, minOverlap, maxMissRate, false);
}

inline void OverlapLayoutConsensus::pairwiseAlignmentOwnself(const unsigned long long readID, const unsigned minOverlap, const double maxMissRate)
{
    if (seqLink[readID].seqLength >= minOverlap)
        pairwiseAlignmentOLC(readID, readID, minOverlap, maxMissRate, true);
}



//////////////////////////////////////////////////////////////////////////////////////
// exec pairwise alignment
//////////////////////////////////////////////////////////////////////////////////////
OverlapLayoutConsensus::ALIGN_STATUS OverlapLayoutConsensus::pairwiseAlignmentOLC(const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const double maxMissRate, const bool ownself)
{
/*
std::cout << ">seq1" << std::endl;
for (auto it = gapInfo->seq[readID1].begin(); it != gapInfo->seq[readID1].end(); ++it)
	std::cout << platanus::Bin2Char(*it);
std::cout << std::endl;
std::cout << ">seq2" << std::endl;
for (auto it = gapInfo->seq[readID2].begin(); it != gapInfo->seq[readID2].end(); ++it)
	std::cout << platanus::Bin2Char(*it);
std::cout << std::endl;
*/
    vector<vector<Score> > alignmentTable(seqLink[readID1].seqLength + 1, vector<Score>(seqLink[readID2].seqLength + 1));
    unsigned start = 1;
    unsigned end = seqLink[readID1].seqLength - minOverlap + 1;
    for (unsigned j = 1; j < seqLink[readID2].seqLength + 1; ++j) {
        if (seqLink[readID2].seqLength - j + 1 < minOverlap) ++start;
        if (end != seqLink[readID1].seqLength + 1) ++end;
        for (unsigned i = start; i < end; ++i) {

            // fill alignment table
            // calculate each score
            if (gapInfo->seq[readID1][i-1] == gapInfo->seq[readID2][j-1] && gapInfo->seq[readID1][i-1] != 4) {
                alignmentTable[i][j].score = alignmentTable[i-1][j-1].score + 1;
                alignmentTable[i][j].miss = alignmentTable[i-1][j-1].miss;
            } else {
                alignmentTable[i][j].score = alignmentTable[i-1][j-1].score;
                alignmentTable[i][j].miss = alignmentTable[i-1][j-1].miss + 1;
            }

        }
    }
    // chekc align proper or not
    return checkAlignmentMissRate(alignmentTable, readID1, readID2, minOverlap, maxMissRate, ownself);
//    return checkAlignment(alignmentTable, readID1, readID2, minOverlap, maxEditDistance, ownself);
}



//////////////////////////////////////////////////////////////////////////////////////
// check alignment result is satisfied enougph aligned
// enough aligned means much overlapped and few missmatch
//////////////////////////////////////////////////////////////////////////////////////
OverlapLayoutConsensus::ALIGN_STATUS OverlapLayoutConsensus::checkAlignment(vector<vector<Score> > &table, const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const unsigned maxEditDistance, const bool ownself)
{
    Score maxScoreI, maxScoreJ;
    ALIGN_STATUS align = NOT_ALIGNED;
    std::pair<unsigned, unsigned> positionI, positionJ;
    unsigned lengthI = seqLink[readID1].seqLength;
    unsigned lengthJ = seqLink[readID2].seqLength;

    unsigned i, j;
    unsigned addTerminal = ownself ? 0 : 1;


    // check overlap from j-i
    // seq i                   -------------------
    // seq j -----------------------
    j = lengthJ;
    for (i = minOverlap; i < lengthI + addTerminal; ++i) {
        if (table[i][j].score >= maxScoreI.score && table[i][j].miss <= maxEditDistance) {
            maxScoreI.score = table[i][j].score;
            maxScoreI.miss  = table[i][j].miss;
            positionI = std::pair<unsigned, unsigned>(i, j);
        }
    }

    // check overlap from i-j
    // seq i ----------------------
    // seq j                 -----------------------
    i = lengthI;
    for (j = minOverlap; j < lengthJ + addTerminal; ++j) {
        if (table[i][j].score >= maxScoreJ.score && table[i][j].miss <= maxEditDistance) {
            maxScoreJ.score = table[i][j].score;
            maxScoreJ.miss  = table[i][j].miss;
            positionJ = std::pair<unsigned, unsigned>(i, j);
        }
    }
/*
std::cout << "ownself = " << (int)ownself << std::endl;
std::cout << "scoreI = " << (int)maxScoreI.score << std::endl;
std::cout << "positionI = " << (int)positionI.first << ", " << (int)positionI.second << std::endl;
std::cout << "scoreJ = " << (int)maxScoreJ.score << std::endl;
std::cout << "positionJ = " << (int)positionJ.first << ", " << (int)positionJ.second << std::endl;
*/
    short scoreSub = maxScoreI.score - maxScoreJ.score;
    if (std::max(maxScoreI.score, maxScoreJ.score) >= minOverlap) {
        align = ALIGNED;
        if ((scoreSub > 0 && positionI.first >= lengthJ)
            || (scoreSub == 0 && positionI.first == lengthI)) { }
        else if (scoreSub < 0 && positionJ.second >= lengthI) { }
        else {
            if (ownself)
                seqLink[readID1].ownOverlap = true;
            else {
                if (scoreSub >= 0) updateSeqLink(readID1, readID2, positionI, maxScoreI.miss, TO1);
                if (scoreSub <= 0) updateSeqLink(readID1, readID2, positionJ, maxScoreJ.miss, TO2);
            }
        }
    }
    return align;
}


OverlapLayoutConsensus::ALIGN_STATUS OverlapLayoutConsensus::checkAlignmentMissRate(vector<vector<Score> > &table, const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const double maxMissRate, const bool ownself)
{
    Score maxScoreI, maxScoreJ;
    ALIGN_STATUS align = NOT_ALIGNED;
    std::pair<unsigned, unsigned> positionI, positionJ;
    unsigned lengthI = seqLink[readID1].seqLength;
    unsigned lengthJ = seqLink[readID2].seqLength;

    unsigned i, j;
    unsigned addTerminal = ownself ? 0 : 1;


    // check overlap from j-i
    // seq i                   -------------------
    // seq j -----------------------
    j = lengthJ;
    for (i = minOverlap; i < lengthI + addTerminal; ++i) {
        if (table[i][j].score >= maxScoreI.score && static_cast<double>(table[i][j].miss) / i <= maxMissRate) {
            maxScoreI.score = table[i][j].score;
            maxScoreI.miss  = table[i][j].miss;
            positionI = std::pair<unsigned, unsigned>(i, j);
        }
    }

    // check overlap from i-j
    // seq i ----------------------
    // seq j                 -----------------------
    i = lengthI;
    for (j = minOverlap; j < lengthJ + addTerminal; ++j) {
        if (table[i][j].score >= maxScoreJ.score && static_cast<double>(table[i][j].miss) / j <= maxMissRate) {
            maxScoreJ.score = table[i][j].score;
            maxScoreJ.miss  = table[i][j].miss;
            positionJ = std::pair<unsigned, unsigned>(i, j);
        }
    }

/*
std::cout << "ownself = " << (int)ownself << std::endl;
std::cout << "scoreI = " << (int)maxScoreI.score << std::endl;
std::cout << "positionI = " << (int)positionI.first << ", " << (int)positionI.second << std::endl;
std::cout << "scoreJ = " << (int)maxScoreJ.score << std::endl;
std::cout << "positionJ = " << (int)positionJ.first << ", " << (int)positionJ.second << std::endl;
*/

    short scoreSub = maxScoreI.score - maxScoreJ.score;
    if (std::max(maxScoreI.score, maxScoreJ.score) >= minOverlap) {
        align = ALIGNED;
        if ((scoreSub > 0 && positionI.first >= lengthJ)
            || (scoreSub == 0 && positionI.first == lengthI)) { }
        else if (scoreSub < 0 && positionJ.second >= lengthI) { }
        else {
            if (ownself)
                seqLink[readID1].ownOverlap = true;
            else {
                if (scoreSub >= 0) updateSeqLink(readID1, readID2, positionI, maxScoreI.miss, TO1);
                if (scoreSub <= 0) updateSeqLink(readID1, readID2, positionJ, maxScoreJ.miss, TO2);
            }
        }
    }
    return align;
}



//////////////////////////////////////////////////////////////////////////////////////
// update link information
// seq information has omly tail(right) side
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::updateSeqLink(const unsigned long long readID1, const unsigned long long readID2, const std::pair<unsigned, unsigned> &position, const unsigned miss, const OverlapLayoutConsensus::OVERLAP pattern)
{
    unsigned difference = std::min(position.first, position.second);

    LinkInformation link(difference, miss);

    switch (pattern) {
        case TO1:
        {
            link.seqID = readID1;
            link.seqLength = seqLink[readID1].seqLength;
            seqLink[readID2].tailLink.emplace_back(link);
            break;
        }
        case TO2:
        {
            link.seqID = readID2;
            link.seqLength = seqLink[readID2].seqLength;
            seqLink[readID1].tailLink.emplace_back(link);
            break;
        }
        case SAME:
        {
            if (position.first > position.second) {
                link.seqID = readID2;
                link.seqLength = seqLink[readID2].seqLength;
                link.overlapLength = seqLink[readID1].seqLength - position.first - difference;
                seqLink[readID1].tailLink.emplace_back(link);
            } else {
                link.seqID = readID1;
                link.seqLength = seqLink[readID1].seqLength;
                link.overlapLength = seqLink[readID2].seqLength - position.second - difference;
                seqLink[readID2].tailLink.emplace_back(link);
            }
            break;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// layout and consensus phase




//////////////////////////////////////////////////////////////////////////////////////
// do greedy extension
//////////////////////////////////////////////////////////////////////////////////////
bool OverlapLayoutConsensus::greedyExtension(const double threshold)
{
    this->clearObject();
    const unsigned topID  = getMostOverlapHeadSeqID();
    bool closed = false;
    seqLayoutInfo[HEADID].startPosition = -1 * GapClose::HEAD_TAIL_SEQ_LEN;
    seqLayoutInfo[topID].startPosition  = -1 * seqLink[HEADID].tailLink.front().overlapLength;
    unsigned numSeqUsed = 0;
//    REASON reason;


    // if topID equals tail seq, closed gap
    if (isTailSeq(topID)) {
//        reason = CLOSED;
        seqLength = seqLayoutInfo[topID].startPosition - (seqLayoutInfo[topID].startPosition > 0 ? 0 : 1);
        seq.resize(std::max(seqLength, 0));
        closed = true;
        return closed;
    }

    // exec consensus first
    greedyFirst(topID, numSeqUsed);

    ObservedSeq observed(topID, gapInfo->seq[topID], 0);
    observedSeqInfo.emplace_back(std::move(observed));

    long endPoint = seqLength + seqLayoutInfo[topID].startPosition;
    seq.resize(endPoint);


    for (long i = seqLayoutInfo[topID].startPosition; i < endPoint; ++i) {
        vector<unsigned> newSeq = popCandidate(i);
        for (long j = 0; j < static_cast<long>(newSeq.size()); ++j) {
            // update and check layout confliction
            ACTION action = updateSeqLayoutInfo(newSeq[j], numSeqUsed, i);
            moveFromCandidateToObserved(newSeq[j]);
            switch (action) {
                // stop means conflict (duplication) is found in layout
                case STOP:
                {
//                    reason = DUPLICATION;
                    closed = false;
                    return closed;
                    break;
                }
                // redo means redo layout
                case REDO:
                {
                    if (this->hasHeadSeqLink()) {
                        bool comp = this->greedyExtension(threshold);
                        return comp;
                    } else {
                        closed = false;
                        return closed;
                    }
                    break;
                }
                case CONTINUE:
                {

                    if (isTailSeq(newSeq[j])) {
//                        reason = CLOSED;
                        seqLength = i > 0 ? i : i - 1;
                        seq.resize(std::max(seqLength, 0));
                        closed = true;
                        return closed;
                    } else if (isOwnOverlap(newSeq[j])) {
//                        reason = OWN_OVERLAP;
                        closed = false;
                        return closed;
                    }
                    break;
                }
            }
        }
        // consensus
        long nowSeqPosition = i > 0 ? i : 0;
        seq[nowSeqPosition] = consensusBase(threshold);

        if (notConsensus(seq[nowSeqPosition])) {
            if (i > 0) {
                seqLength = i;
//                reason = UNKNOWN_BASE;
                closed = false;
                return closed;
            }
        }

        updateObservedSeqInfo();

        if (notExtension()) {
            seqLength = i > 0 ? i : 0;
//            reason = NO_SEQ;
            closed = false;
            return closed;
        }
    }
    return closed;
}


//////////////////////////////////////////////////////////////////////////////////////
// consensus
//////////////////////////////////////////////////////////////////////////////////////
char OverlapLayoutConsensus::consensusBase(const double threshold) const
{
    vector<unsigned> base(5, 0);
    unsigned total = 0;
    unsigned max = 0;
    char maxBase = 0;
    auto listIt = observedSeqInfo.begin();
    auto listEnd = observedSeqInfo.end();

    for (; listIt != listEnd; ++listIt) {
        if (listIt->seq[listIt->nowPosition] > 4 || listIt->seq[listIt->nowPosition] < 0) {
            throw platanus::BaseError();
        }
        ++base[listIt->seq[listIt->nowPosition]];
        ++total;
    }
    for (char i = 0; i < 4; ++i) {
        if (base[i] > max) {
            max = base[i];
            maxBase = i;
        }
    }
    return (static_cast<double>(max) / total) >= threshold ? maxBase : 4;

}



//////////////////////////////////////////////////////////////////////////////////////
// update observedSeqInfo
// plus 1 seq position
// if seq position reaches seq end, erase this seq
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::updateObservedSeqInfo(void)
{
    auto listIt = observedSeqInfo.begin();
    auto listEnd = observedSeqInfo.end();

    for (; listIt != listEnd; ) {
        ++listIt->nowPosition;
        if (listIt->nowPosition == listIt->seq.size()) {
            listIt = observedSeqInfo.erase(listIt);
            continue;
        }
        ++listIt;
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// set first candidate seqs used greedy extension
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::greedyFirst(const unsigned topID, unsigned &numSeqUsed)
{
    auto listIt = seqLink[topID].tailLink.begin();
    auto listEnd = seqLink[topID].tailLink.end();

    seqLayoutInfo[topID].used = true;
    for (; listIt != listEnd; ++listIt) {
        addCandidateSeq(*listIt, seqLink[topID].seqLength, seqLayoutInfo[topID].startPosition);
        ++numSeqUsed;
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// set candidate seq actually
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::addCandidateSeq(const LinkInformation &link, const unsigned seqLength, const unsigned startPosition)
{
    seqLayoutInfo[link.seqID].startPosition = startPosition + seqLength - link.overlapLength;
    seqLayoutInfo[link.seqID].used = true;
    candidateSeqList.emplace_back(link.seqID);
}




//////////////////////////////////////////////////////////////////////////////////////
// update seqLayoutInfo (add new child seq in candidate)
//////////////////////////////////////////////////////////////////////////////////////
OverlapLayoutConsensus::ACTION OverlapLayoutConsensus::updateSeqLayoutInfo(const unsigned parentID, unsigned &num, const long nowDistance)
{

    ACTION action = CONTINUE;
    auto listIt = seqLink[parentID].tailLink.begin();
    auto listEnd = seqLink[parentID].tailLink.end();

    // check whether position has contradiction
    for (; listIt != listEnd; ++listIt) {
        POSITION_CHECK check = checkDistance(listIt->seqID, listIt->overlapLength, parentID, nowDistance);

        if (check == INCORRECT) {
            // decide which position is proper (parents or child)
            PROPER_LINK which = selectProperLink(parentID, listIt->seqID);
             if (checkLinkFar(seqLayoutInfo[listIt->seqID].startPosition, nowDistance + seqLink[parentID].seqLength - listIt->overlapLength, seqLink[parentID].seqLength)) which = UNKNOWN;
            switch (which) {
                case OWN:
                {
                    auto observedIt = observedSeqInfo.begin();
                    auto observedEnd = observedSeqInfo.end();
                    for (; observedIt != observedEnd; ++observedIt) {
                        deleteLink(observedIt->seqID, listIt->seqID);
                    }
                    action = REDO;
                    break;
                }
                case UNKNOWN:
                {
                    return STOP;
                    break;
                }
                case OTHERS:
                {
                    break;
                }
            }
        } else if (check == NOT_USED) {
            updateStartPosition(listIt->seqID, listIt->overlapLength, nowDistance, seqLink[parentID].seqLength);
            candidateSeqList.emplace_back(listIt->seqID);
            ++num;

        }
    }
    return action;
}



//////////////////////////////////////////////////////////////////////////////////////
// check child seq start position is appropriate
//////////////////////////////////////////////////////////////////////////////////////
OverlapLayoutConsensus::POSITION_CHECK OverlapLayoutConsensus::checkDistance(const unsigned childID, const unsigned overlapLength, const unsigned parentID, const long nowDistance) const
{
    if (seqLayoutInfo[childID].used == false)
        return NOT_USED;
    else if (isSameStartPositionAmongDifferentParents(childID, overlapLength, parentID, nowDistance))
        return CORRECT;
    else
        return INCORRECT;
}


//////////////////////////////////////////////////////////////////////////////////////
// judge child seq start position is same or not among different parents which has link the child.
//////////////////////////////////////////////////////////////////////////////////////
inline bool OverlapLayoutConsensus::isSameStartPositionAmongDifferentParents(const unsigned childID, const unsigned overlapLength, const unsigned parentID, const long nowDistance) const
{
    return seqLayoutInfo[childID].startPosition == static_cast<int>(nowDistance + seqLink[parentID].seqLength - overlapLength);
}




//////////////////////////////////////////////////////////////////////////////////////
// decide which parent link is more proper
// child                 -------------------    -|
// parent 1 ----------------------               |- which link better ?
// child                     -------------      -|
// parent 2     ----------------------
// the lowest missmathcs link is chosen
//////////////////////////////////////////////////////////////////////////////////////
OverlapLayoutConsensus::PROPER_LINK OverlapLayoutConsensus::selectProperLink(const unsigned parentID, const unsigned childID) const
{
    auto linkIt = seqLink[parentID].tailLink.begin();
    auto linkEnd = seqLink[parentID].tailLink.end();
    unsigned parentMiss = 0;
    bool win = false;

    for (; linkIt != linkEnd; ++linkIt) {
        if (linkIt->seqID == childID) {
            parentMiss = linkIt->miss;
            break;
        }
    }
    auto infoIt = observedSeqInfo.begin();
    auto infoEnd = observedSeqInfo.end();

    for (; infoIt != infoEnd; ++infoIt) {
        linkIt = seqLink[infoIt->seqID].tailLink.begin();
        linkEnd = seqLink[infoIt->seqID].tailLink.end();
        for (; linkIt != linkEnd; ++linkIt) {
            if (linkIt->seqID == childID) {
                if (linkIt->miss < parentMiss) {
                    return OTHERS;
                } else if (linkIt->miss == parentMiss) {
                    return UNKNOWN;
                } else {
                    win = true;
                    break;
                }
            }
        }
    }
    return win ? OWN : OTHERS;
}


//////////////////////////////////////////////////////////////////////////////////////
// check parent links are too far
// child                 -------------------
// parent 1 ---------------------
// child                                                             -------------
// parent 2                                         ----------------------
//////////////////////////////////////////////////////////////////////////////////////
bool OverlapLayoutConsensus::checkLinkFar(const int othersLinkStartPositionToChild, const int ownLinkStartPositionToChild, const unsigned ownLength) const
{
    return abs(othersLinkStartPositionToChild - ownLinkStartPositionToChild) >= ownLength;
}



//////////////////////////////////////////////////////////////////////////////////////
// pop new observed seq into candidates
//////////////////////////////////////////////////////////////////////////////////////
vector<unsigned> OverlapLayoutConsensus::popCandidate(const long seqPosition)
{
    vector<unsigned> newSeq;
    auto listIt = candidateSeqList.begin();
    auto listEnd = candidateSeqList.end();
    for (; listIt != listEnd; ) {
        if (seqPosition == seqLayoutInfo[*listIt].startPosition) {
            newSeq.emplace_back(*listIt);
            listIt = candidateSeqList.erase(listIt);
            continue;
        }
        ++listIt;
    }
    return newSeq;
}

//////////////////////////////////////////////////////////////////////////////////////
// update child start position
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::updateStartPosition(const unsigned childID, const int overlap, const int distance, const unsigned seqLength)
{
    seqLayoutInfo[childID].startPosition = distance + seqLength - overlap;
    seqLayoutInfo[childID].used = true;
}


//////////////////////////////////////////////////////////////////////////////////////
// move child seq form candidate to observed
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::moveFromCandidateToObserved(const unsigned newSeqID)
{
    ObservedSeq newObserved(newSeqID, gapInfo->seq[newSeqID], 0);
    observedSeqInfo.emplace_back(std::move(newObserved));
    auto listIt = candidateSeqList.begin();
    auto listEnd = candidateSeqList.end();
    for (; listIt != listEnd; ++listIt) {
        if (*listIt == newSeqID) {
            candidateSeqList.erase(listIt);
            break;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// delete link
//////////////////////////////////////////////////////////////////////////////////////
void OverlapLayoutConsensus::deleteLink(const unsigned parentID, const unsigned childID)
{
    auto linkIt = seqLink[parentID].tailLink.begin();
    auto linkEnd = seqLink[parentID].tailLink.end();
    for (; linkIt != linkEnd; ++linkIt) {
        if (linkIt->seqID == childID) {
            seqLink[parentID].tailLink.erase(linkIt);
            break;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// check newSeq is tail
//////////////////////////////////////////////////////////////////////////////////////
inline bool OverlapLayoutConsensus::isTailSeq(const unsigned newSeq) const
{
    return newSeq == TAILID;

}


//////////////////////////////////////////////////////////////////////////////////////
// check newSeq is overlapeed ownself
//////////////////////////////////////////////////////////////////////////////////////
inline bool OverlapLayoutConsensus::isOwnOverlap(const unsigned newSeq) const
{
    return seqLink[newSeq].ownOverlap;

}

//////////////////////////////////////////////////////////////////////////////////////
// check newSeq is overlapeed ownself
//////////////////////////////////////////////////////////////////////////////////////
inline bool OverlapLayoutConsensus::notConsensus(const char base) const
{
    return base == 4;

}


//////////////////////////////////////////////////////////////////////////////////////
// check observed seq is empty
//////////////////////////////////////////////////////////////////////////////////////
inline bool OverlapLayoutConsensus::notExtension(void) const
{
    return observedSeqInfo.empty();

}


