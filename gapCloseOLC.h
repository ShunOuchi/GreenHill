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

#ifndef GAPCLOSE_OLC_H
#define GAPCLOSE_OLC_H

#include "common.h"
#include "seqlib.h"
#include "gapClose.h"
#include "kmer.h"
#include <list>
#include <unordered_map>



//////////////////////////////////////////////////////////////////////////////////////
// OLC class
//////////////////////////////////////////////////////////////////////////////////////
class OverlapLayoutConsensus
{
private:


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// LinkInformation class
// this class has child seq (called "link" in this programs) information
    struct LinkInformation
    {
        unsigned seqID;
        unsigned overlapLength;
        unsigned miss;
        unsigned seqLength;

        LinkInformation(): seqID(0), overlapLength(0), miss(0), seqLength(0) {}
        LinkInformation(const unsigned diff, const unsigned m): seqID(0), overlapLength(diff), miss(m), seqLength(0) {}
        LinkInformation(const LinkInformation &) = default;
        LinkInformation &operator=(const LinkInformation &) = default;
        ~LinkInformation() = default;

        bool operator<(const LinkInformation &link) const
        {
            short tmp = overlapLength - link.overlapLength;
            if (tmp != 0)
                return tmp > 0;
            else
                return seqLength > link.seqLength;
        }

    };
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// SeqLink class
// this class has parent seq link information
    struct SeqLink
    {
        unsigned seqID;
        bool ownOverlap;
        std::list<LinkInformation> tailLink;
        unsigned seqLength;

        SeqLink(): seqID(0), ownOverlap(false), tailLink(), seqLength(0) {}
        SeqLink(const SeqLink &) = default;
        SeqLink &operator=(const SeqLink &) = default;
        ~SeqLink() = default;

    };
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Score class
// this class has pairwise alignment score whether between seqs have ovelap
    struct Score
    {
        unsigned score;
        unsigned miss;

        Score():score(0), miss(0) {}
        Score(const Score &) = default;
        Score &operator=(const Score &) = default;
        ~Score() = default;
    };
//////////////////////////////////////////////////////////////////////////////////////





    // layout phase

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// ObserbedSeq class
// this class has seq list which is being used base consensus just now
    struct ObservedSeq
    {
        unsigned seqID;
        std::vector<char> seq;
        unsigned nowPosition;

        ObservedSeq() = delete;
        ObservedSeq(const unsigned id, const std::vector<char> &s, const unsigned now): seqID(id), seq(s), nowPosition(now) {}
        ObservedSeq(const ObservedSeq &) = default;
        ObservedSeq &operator=(const ObservedSeq &) = default;
        ~ObservedSeq() = default;
    };
//////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// LayoutInfo class
// this class has seq layout information where seq is placed in gap position
    struct LayoutInfo
    {
        int startPosition; // dist
        bool used;

        LayoutInfo(): startPosition(0), used(false) {}
        LayoutInfo(const LayoutInfo &) = default;
        LayoutInfo &operator=(const LayoutInfo &) = default;
        ~LayoutInfo() = default;
    };
//////////////////////////////////////////////////////////////////////////////////////

    enum ALIGN_STATUS
    { YET = 0, ALIGNED, NOT_ALIGNED};

    enum OVERLAP
    { TO1 = 0, TO2, SAME};

    enum POSITION_CHECK
    { NOT_USED = 0, CORRECT, INCORRECT};

    enum ACTION
    { CONTINUE = 0, REDO, STOP};

    enum PROPER_LINK
    { OWN = 0, OTHERS, UNKNOWN};

    enum REASON
    { CLOSED = 0, NO_SEQ, UNKNOWN_BASE, DUPLICATION, OWN_OVERLAP};




    //////////////////////////////////////////////////////////////////////////////////////
    // constant
    //////////////////////////////////////////////////////////////////////////////////////
    static const unsigned HEADID;
    static const unsigned TAILID;



    const GAP * const gapInfo;
    const unsigned numSeq;
    std::vector<SeqLink> seqLink;

    // use only overlap phase
    const unsigned kmerLength;
    std::unordered_map<std::string, std::list<unsigned long long> > kmerList;


    // use only layout phase
    int seqLength;
    std::vector<char> seq;
    std::list<unsigned> candidateSeqList;
    std::list<ObservedSeq> observedSeqInfo;
    std::vector<LayoutInfo> seqLayoutInfo;



    inline void updateCheckTable(std::vector<std::vector<ALIGN_STATUS> > &checkTable, const unsigned long long readID1, const unsigned long long readID2, ALIGN_STATUS align);
    void initSeqLink(void);
    inline ALIGN_STATUS pairwiseAlignment(const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const double maxMissRate);
    inline void pairwiseAlignmentOwnself(const unsigned long long readID, const unsigned minOverlap, const double maxMissRate);
    ALIGN_STATUS pairwiseAlignmentOLC(const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const double maxMissRate, const bool ownself);
    ALIGN_STATUS checkAlignment(std::vector<std::vector<Score> > &table, const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const unsigned maxEditDistance, const bool ownself);
    ALIGN_STATUS checkAlignmentMissRate(std::vector<std::vector<Score> > &table, const unsigned long long readID1, const unsigned long long readID2, const unsigned minOverlap, const double maxMissRate, const bool ownself);
    void updateSeqLink(const unsigned long long readID1, const unsigned long long readID2, const std::pair<unsigned, unsigned> &position, const unsigned miss, const OVERLAP pattern);

    char consensusBase(const double threshold) const;
    void updateObservedSeqInfo(void);
    void greedyFirst(const unsigned topID, unsigned &numSeqUsed);
    void addCandidateSeq(const LinkInformation &link, const unsigned seqLength, const unsigned startPosition);
    ACTION updateSeqLayoutInfo(const unsigned parentID, unsigned &num, const long nowDistance);
    POSITION_CHECK checkDistance(const unsigned childID, const unsigned overlapLength, const unsigned parentID, const long nowDistance) const;
    bool isSameStartPositionAmongDifferentParents(const unsigned childID, const unsigned overlapLength, const unsigned parentID, const long nowDistance) const;
    PROPER_LINK selectProperLink(const unsigned parentID, const unsigned childID) const;
    bool checkLinkFar(const int othersLinkStartPositionToChild, const int ownLinkStartPositionToChild, const unsigned ownLength) const;
    std::vector<unsigned> popCandidate(const long seqPosition);
    void updateStartPosition(const unsigned childID, const int overlap, const int distance, const unsigned seqLength);
    void moveFromCandidateToObserved(const unsigned newSeqID);
    void deleteLink(const unsigned parentID, const unsigned childID);
    bool isTailSeq(const unsigned newSeq) const;
    bool isOwnOverlap(const unsigned newSeq) const;
    bool notConsensus(const char base) const;
    bool notExtension(void) const;

    bool checkAlreadyAlign(const std::vector<std::vector<ALIGN_STATUS> > &table, unsigned long long read1, unsigned long long read2) const
    { return table[read1][read2 - read1 - 1] == YET; }

    unsigned getMostOverlapHeadSeqID(void) const { return seqLink[HEADID].tailLink.front().seqID; }
    void deleteKmerList(void)
    { kmerList.clear(); }

    void clearObject(void)
    {
        candidateSeqList.clear();
        observedSeqInfo.clear();
        seqLayoutInfo.clear();
        seqLayoutInfo.resize(numSeq);
    }
public:

    OverlapLayoutConsensus(const unsigned kmer, GAP *gap): gapInfo(gap), numSeq(gapInfo->seq.size()), seqLink(numSeq), kmerLength(kmer), kmerList(), seqLength(gap->length), seq(), candidateSeqList(), observedSeqInfo(), seqLayoutInfo(numSeq) {}
    OverlapLayoutConsensus(const OverlapLayoutConsensus &olc) = delete;
    OverlapLayoutConsensus &operator=(const OverlapLayoutConsensus &olc) = delete;
    ~OverlapLayoutConsensus() = default;

    void makeKmerTable(void);
    void makeOverlapGraph(const unsigned minOverlap, const double maxMissRate);
    bool greedyExtension(const double threshold);

    // getter and setter
    unsigned getNumSeq(void) const {return numSeq; }
    std::vector<char> getLayoutSeq(void) const {return seq; }
    int getLayoutLength(void) const {return seqLength; }

    bool hasHeadSeqLink(void) const { return seqLink[HEADID].tailLink.size() != 0; }
};


#endif
