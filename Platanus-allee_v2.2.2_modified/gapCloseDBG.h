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

#ifndef GAPCLOSE_DBG_H
#define GAPCLOSE_DBG_H

#include "common.h"
#include "gapClose.h"
#include "counter.h"
#include "graph.h"


// using tag dispatch

template <int v>
struct Int2Type
{
    enum {value = v};
};


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// de bruijn graph gap closing class
// low kmer length and high kmer length are already defined in compile, using template
template <unsigned M, unsigned N>
class GapCloseDBG
{
private:
    static const double BRANCH_THRESHOLD;
    static const double BUBBLE_THRESHOLD;
    BruijnGraph<Kmer31> kmer31Graph;
    BruijnGraph<KmerN<Binstr63> > kmer63Graph;
    BruijnGraph<KmerN<Binstr95> > kmer95Graph;
    BruijnGraph<KmerN<Binstr127> > kmer127Graph;
    BruijnGraph<KmerN<Binstr159> > kmer159Graph;
    BruijnGraph<KmerN<binstr_t> > kmerNGraph;
    Counter<Kmer31> kmer31Counter;
    Counter<KmerN<Binstr63> > kmer63Counter;
    Counter<KmerN<Binstr95> > kmer95Counter;
    Counter<KmerN<Binstr127> > kmer127Counter;
    Counter<KmerN<Binstr159> > kmer159Counter;
    Counter<KmerN<binstr_t> > kmerNCounter;
    const unsigned minKmer;
    const unsigned maxKmer;
    const long minCoverage;
    const unsigned minOverlap;
    const double maxMissRate;
    const GAP * const gapInfo;
    int seqLength;
    std::vector<char> seq;
	int remainedGap;
    enum KMERTYPE
    { KMER31, KMER63, KMER95, KMER127, KMER159, KMERN};
    static const int minType = (M <=  32) ? KMER31  :
                               (M <=  64) ? KMER63  :
                               (M <=  96) ? KMER95  :
                               (M <= 128) ? KMER127 :
                               (M <= 160) ? KMER159 : KMERN;
    static const int maxType = (N <=  32) ? KMER31  :
                               (N <=  64) ? KMER63  :
                               (N <=  96) ? KMER95  :
                               (N <= 128) ? KMER127 :
                               (N <= 160) ? KMER159 : KMERN;

    void lowKmerCount(void)
    {
        this->kmerCount(minKmer, Int2Type<minType>());
    }
    void highKmerCount(void)
    {
        this->kmerCount(maxKmer, Int2Type<maxType>());
    }
    void kmerCount(const unsigned kmer, Int2Type<KMER31>)
    {
        kmer31Counter.setKmerLength(kmer);
        kmer31Counter.countKmerForGapClose(gapInfo->seq);
    }
    void kmerCount(const unsigned kmer, Int2Type<KMER63>)
    {
        kmer63Counter.setKmerLength(kmer);
        kmer63Counter.countKmerForGapClose(gapInfo->seq);
    }
    void kmerCount(const unsigned kmer, Int2Type<KMER95>)
    {
        kmer95Counter.setKmerLength(kmer);
        kmer95Counter.countKmerForGapClose(gapInfo->seq);
    }
    void kmerCount(const unsigned kmer, Int2Type<KMER127>)
    {
        kmer127Counter.setKmerLength(kmer);
        kmer127Counter.countKmerForGapClose(gapInfo->seq);
    }
    void kmerCount(const unsigned kmer, Int2Type<KMER159>)
    {
        kmer159Counter.setKmerLength(kmer);
        kmer159Counter.countKmerForGapClose(gapInfo->seq);
    }
    void kmerCount(const unsigned kmer, Int2Type<KMERN>)
    {
        kmerNCounter.setKmerLength(kmer);
        kmerNCounter.countKmerForGapClose(gapInfo->seq);
    }



    void lowKmerMakeBruijnGraph(void)
    {
        this->makeBruijnGraph(minKmer, Int2Type<minType>());
    }
    void highKmerMakeBruijnGraph(void)
    {
        this->makeBruijnGraph(maxKmer, Int2Type<maxType>());
    }
    void makeBruijnGraph(const unsigned kmerLength, Int2Type<KMER31>)
    {
        kmer31Graph.setKmerLength(kmerLength);
        kmer31Graph.makeBruijnGraphForGapClose(kmer31Counter, minCoverage);
    }
    void makeBruijnGraph(const unsigned kmerLength, Int2Type<KMER63>)
    {
        kmer63Graph.setKmerLength(kmerLength);
        kmer63Graph.makeBruijnGraphForGapClose(kmer63Counter, minCoverage);
    }
    void makeBruijnGraph(const unsigned kmerLength, Int2Type<KMER95>)
    {
        kmer95Graph.setKmerLength(kmerLength);
        kmer95Graph.makeBruijnGraphForGapClose(kmer95Counter, minCoverage);
    }
    void makeBruijnGraph(const unsigned kmerLength, Int2Type<KMER127>)
    {
        kmer127Graph.setKmerLength(kmerLength);
        kmer127Graph.makeBruijnGraphForGapClose(kmer127Counter, minCoverage);
    }
    void makeBruijnGraph(const unsigned kmerLength, Int2Type<KMER159>)
    {
        kmer159Graph.setKmerLength(kmerLength);
        kmer159Graph.makeBruijnGraphForGapClose(kmer159Counter, minCoverage);
    }
    void makeBruijnGraph(const unsigned kmerLength, Int2Type<KMERN>)
    {
        kmerNGraph.setKmerLength(kmerLength);
        kmerNGraph.makeBruijnGraphForGapClose(kmerNCounter, minCoverage);
    }




    void mergeGraph(void)
    {
        this->mergeGraph(Int2Type<minType>(), Int2Type<maxType>());
    }
    void mergeGraph(Int2Type<KMER31>, Int2Type<KMER63>)
    {
        kmer63Counter.clearOccurrenceTableGap();
        kmer31Graph.saveLargeKmerForGapClose(kmer63Counter);
        kmer63Graph.saveLargeKmerForGapClose(kmer63Counter);
    }
    void mergeGraph(Int2Type<KMER31>, Int2Type<KMER95>)
    {
        kmer95Counter.clearOccurrenceTableGap();
        kmer31Graph.saveLargeKmerForGapClose(kmer95Counter);
        kmer95Graph.saveLargeKmerForGapClose(kmer95Counter);
    }
    void mergeGraph(Int2Type<KMER31>, Int2Type<KMER127>)
    {
        kmer127Counter.clearOccurrenceTableGap();
        kmer31Graph.saveLargeKmerForGapClose(kmer127Counter);
        kmer127Graph.saveLargeKmerForGapClose(kmer127Counter);
    }
    void mergeGraph(Int2Type<KMER31>, Int2Type<KMER159>)
    {
        kmer159Counter.clearOccurrenceTableGap();
        kmer31Graph.saveLargeKmerForGapClose(kmer159Counter);
        kmer159Graph.saveLargeKmerForGapClose(kmer159Counter);
    }
    void mergeGraph(Int2Type<KMER31>, Int2Type<KMERN>)
    {
        kmerNCounter.clearOccurrenceTableGap();
        kmer31Graph.saveLargeKmerForGapClose(kmerNCounter);
        kmerNGraph.saveLargeKmerForGapClose(kmerNCounter);
    }
    void mergeGraph(Int2Type<KMER63>, Int2Type<KMER95>)
    {
        kmer95Counter.clearOccurrenceTableGap();
        kmer63Graph.saveLargeKmerForGapClose(kmer95Counter);
        kmer95Graph.saveLargeKmerForGapClose(kmer95Counter);
    }
    void mergeGraph(Int2Type<KMER63>, Int2Type<KMER127>)
    {
        kmer127Counter.clearOccurrenceTableGap();
        kmer63Graph.saveLargeKmerForGapClose(kmer127Counter);
        kmer127Graph.saveLargeKmerForGapClose(kmer127Counter);
    }
    void mergeGraph(Int2Type<KMER63>, Int2Type<KMER159>)
    {
        kmer159Counter.clearOccurrenceTableGap();
        kmer63Graph.saveLargeKmerForGapClose(kmer159Counter);
        kmer159Graph.saveLargeKmerForGapClose(kmer159Counter);
    }
    void mergeGraph(Int2Type<KMER63>, Int2Type<KMERN>)
    {
        kmerNCounter.clearOccurrenceTableGap();
        kmer63Graph.saveLargeKmerForGapClose(kmerNCounter);
        kmerNGraph.saveLargeKmerForGapClose(kmerNCounter);
    }
    void mergeGraph(Int2Type<KMER95>, Int2Type<KMER127>)
    {
        kmer127Counter.clearOccurrenceTableGap();
        kmer95Graph.saveLargeKmerForGapClose(kmer127Counter);
        kmer127Graph.saveLargeKmerForGapClose(kmer127Counter);
    }
    void mergeGraph(Int2Type<KMER95>, Int2Type<KMER159>)
    {
        kmer159Counter.clearOccurrenceTableGap();
        kmer95Graph.saveLargeKmerForGapClose(kmer159Counter);
        kmer159Graph.saveLargeKmerForGapClose(kmer159Counter);
    }
    void mergeGraph(Int2Type<KMER95>, Int2Type<KMERN>)
    {
        kmerNCounter.clearOccurrenceTableGap();
        kmer95Graph.saveLargeKmerForGapClose(kmerNCounter);
        kmerNGraph.saveLargeKmerForGapClose(kmerNCounter);
    }
    void mergeGraph(Int2Type<KMER127>, Int2Type<KMER159>)
    {
        kmer159Counter.clearOccurrenceTableGap();
        kmer127Graph.saveLargeKmerForGapClose(kmer159Counter);
        kmer159Graph.saveLargeKmerForGapClose(kmer159Counter);
    }
    void mergeGraph(Int2Type<KMER127>, Int2Type<KMERN>)
    {
        kmerNCounter.clearOccurrenceTableGap();
        kmer127Graph.saveLargeKmerForGapClose(kmerNCounter);
        kmerNGraph.saveLargeKmerForGapClose(kmerNCounter);
    }
    void mergeGraph(Int2Type<KMER159>, Int2Type<KMERN>)
    {
        kmerNCounter.clearOccurrenceTableGap();
        kmer159Graph.saveLargeKmerForGapClose(kmerNCounter);
        kmerNGraph.saveLargeKmerForGapClose(kmerNCounter);
    }


    unsigned long long lowKmerCutBranch(void)
    {
        return this->cutBranch(Int2Type<minType>());
    }
    unsigned long long highKmerCutBranch(void)
    {
        return this->cutBranch(Int2Type<maxType>());
    }
    unsigned long long cutBranch(Int2Type<KMER31>)  {return kmer31Graph.cutBranch(); }
    unsigned long long cutBranch(Int2Type<KMER63>)  {return kmer63Graph.cutBranch(); }
    unsigned long long cutBranch(Int2Type<KMER95>)  {return kmer95Graph.cutBranch(); }
    unsigned long long cutBranch(Int2Type<KMER127>) {return kmer127Graph.cutBranch(); }
    unsigned long long cutBranch(Int2Type<KMER159>) {return kmer159Graph.cutBranch(); }
    unsigned long long cutBranch(Int2Type<KMERN>)   {return kmerNGraph.cutBranch(); }



    unsigned long long lowKmerCrushBubble(void)
    {
        return this->crushBubble(Int2Type<minType>());
    }
    unsigned long long highKmerCrushBubble(void)
    {
        return this->crushBubble(Int2Type<maxType>());
    }
    unsigned long long crushBubble(Int2Type<KMER31>)  {return kmer31Graph.crushBubble(UINT16_MAX); }
    unsigned long long crushBubble(Int2Type<KMER63>)  {return kmer63Graph.crushBubble(UINT16_MAX); }
    unsigned long long crushBubble(Int2Type<KMER95>)  {return kmer95Graph.crushBubble(UINT16_MAX); }
    unsigned long long crushBubble(Int2Type<KMER127>) {return kmer127Graph.crushBubble(UINT16_MAX); }
    unsigned long long crushBubble(Int2Type<KMER159>) {return kmer159Graph.crushBubble(UINT16_MAX); }
    unsigned long long crushBubble(Int2Type<KMERN>)   {return kmerNGraph.crushBubble(UINT16_MAX); }



    void lowKmerConcatinateNodes(void)
    {
        this->concatinateNodes(Int2Type<minType>());
    }
    void highKmerConcatinateNodes(void)
    {
        this->concatinateNodes(Int2Type<maxType>());
    }
    void concatinateNodes(Int2Type<KMER31>) {kmer31Graph.concatinateNodes(); }
    void concatinateNodes(Int2Type<KMER63>) {kmer63Graph.concatinateNodes(); }
    void concatinateNodes(Int2Type<KMER95>) {kmer95Graph.concatinateNodes(); }
    void concatinateNodes(Int2Type<KMER127>) {kmer127Graph.concatinateNodes(); }
    void concatinateNodes(Int2Type<KMER159>) {kmer159Graph.concatinateNodes(); }
    void concatinateNodes(Int2Type<KMERN>) {kmerNGraph.concatinateNodes(); }

    bool isClosedGap(Int2Type<KMER31>) {return closedGap(kmer31Graph); }
    bool isClosedGap(Int2Type<KMER63>) {return closedGap(kmer63Graph); }
    bool isClosedGap(Int2Type<KMER95>) {return closedGap(kmer95Graph); }
    bool isClosedGap(Int2Type<KMER127>) {return closedGap(kmer127Graph); }
    bool isClosedGap(Int2Type<KMER159>) {return closedGap(kmer159Graph); }
    bool isClosedGap(Int2Type<KMERN>) {return closedGap(kmerNGraph); }

    bool isClosedGapPartial(Int2Type<KMER31>) {return closedGapPartial(kmer31Graph); }
    bool isClosedGapPartial(Int2Type<KMER63>) {return closedGapPartial(kmer63Graph); }
    bool isClosedGapPartial(Int2Type<KMER95>) {return closedGapPartial(kmer95Graph); }
    bool isClosedGapPartial(Int2Type<KMER127>) {return closedGapPartial(kmer127Graph); }
    bool isClosedGapPartial(Int2Type<KMER159>) {return closedGapPartial(kmer159Graph); }
    bool isClosedGapPartial(Int2Type<KMERN>) {return closedGapPartial(kmerNGraph); }


    void save(FILE *contigFP, Int2Type<KMER31>) const { kmer31Graph.saveContigSimple(contigFP, 1.0); }
    void save(FILE *contigFP, Int2Type<KMER63>) const { kmer63Graph.saveContigSimple(contigFP, 1.0); }
    void save(FILE *contigFP, Int2Type<KMER95>) const { kmer95Graph.saveContigSimple(contigFP, 1.0); }
    void save(FILE *contigFP, Int2Type<KMER127>) const { kmer127Graph.saveContigSimple(contigFP, 1.0); }
    void save(FILE *contigFP, Int2Type<KMER159>) const { kmer159Graph.saveContigSimple(contigFP, 1.0); }
    void save(FILE *contigFP, Int2Type<KMERN>) const { kmerNGraph.saveContigSimple(contigFP, 1.0); }


    inline void assembleLowKmer(void);
    inline void assembleHighKmer(void);

    std::pair<long, double> calcMissmatchLeftEdgeSeq(const Straight &straight);
    std::pair<long, double> calcMissmatchRightEdgeSeq(const Straight &straight);
    template <typename KMER> bool closedGap(BruijnGraph<KMER> &graph);
    template <typename KMER> bool closedGapPartial(BruijnGraph<KMER> &graph);


    bool extendEdge(Int2Type<KMER31>) {return this->extendEdgeEntity(kmer31Graph); }
    bool extendEdge(Int2Type<KMER63>) {return this->extendEdgeEntity(kmer63Graph); }
    bool extendEdge(Int2Type<KMER95>) {return this->extendEdgeEntity(kmer95Graph); }
    bool extendEdge(Int2Type<KMER127>) {return this->extendEdgeEntity(kmer127Graph); }
    bool extendEdge(Int2Type<KMER159>) {return this->extendEdgeEntity(kmer159Graph); }
    bool extendEdge(Int2Type<KMERN>) {return this->extendEdgeEntity(kmerNGraph); }

    template <typename KMER> bool extendEdgeEntity(BruijnGraph<KMER> &graph);

public:

    GapCloseDBG() = delete;
    GapCloseDBG(const unsigned min, const unsigned max, const long cov, const unsigned overlap, const double miss, const GAP *gap)
    : kmer31Graph(), kmer63Graph(), kmer95Graph(), kmer127Graph(), kmer159Graph(), kmerNGraph(), kmer31Counter(), kmer63Counter(), kmer95Counter(), kmer127Counter(), kmer159Counter(), kmerNCounter(), minKmer(min), maxKmer(max), minCoverage(cov), minOverlap(overlap), maxMissRate(miss), gapInfo(gap), seqLength(0), seq(), remainedGap(0)
    {
        kmer31Graph.setBubbleAndBranch(BUBBLE_THRESHOLD, BRANCH_THRESHOLD);
        kmer63Graph.setBubbleAndBranch(BUBBLE_THRESHOLD, BRANCH_THRESHOLD);
        kmer95Graph.setBubbleAndBranch(BUBBLE_THRESHOLD, BRANCH_THRESHOLD);
        kmer127Graph.setBubbleAndBranch(BUBBLE_THRESHOLD, BRANCH_THRESHOLD);
        kmer159Graph.setBubbleAndBranch(BUBBLE_THRESHOLD, BRANCH_THRESHOLD);
        kmerNGraph.setBubbleAndBranch(BUBBLE_THRESHOLD, BRANCH_THRESHOLD);
    }

    GapCloseDBG(const GapCloseDBG &) = delete;
    GapCloseDBG &operator=(const GapCloseDBG &) = delete;
    ~GapCloseDBG() = default;

    bool enableGapClose(void);
    void gapAssemble(void);
    bool isClosedGap(void)
    {
        return this->isClosedGap(Int2Type<maxType>());
    }

    bool isClosedGapPartial(void)
    {
        return this->isClosedGapPartial(Int2Type<maxType>());
    }

    std::vector<char> getSeq(void) const {return seq; }
    int getSeqLength(void) const {return seqLength; }
	int getRemainedGap(void) const {return remainedGap; }
    void saveContig(FILE *contigFP) const;
    bool extendEdge(void)
    {
        return extendEdge(Int2Type<maxType>());
    }

};
//////////////////////////////////////////////////////////////////////////////////////

template <unsigned M, unsigned N>
const double GapCloseDBG<M, N>::BUBBLE_THRESHOLD = 0.1;

template <unsigned M, unsigned N>
const double GapCloseDBG<M, N>::BRANCH_THRESHOLD = 0.5;


//////////////////////////////////////////////////////////////////////////////////////
// check enable gap closed
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
bool GapCloseDBG<M, N>::enableGapClose(void)
{
    return (gapInfo->headSeq.size() >= minOverlap) && (gapInfo->tailSeq.size() >= minOverlap) && gapInfo->seq.size() >= 2;
}



//////////////////////////////////////////////////////////////////////////////////////
// exec De Bruijn assemble around gaps
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
void GapCloseDBG<M, N>::gapAssemble(void)
{
    unsigned long long numDelete = 0;


    this->assembleLowKmer();
    this->assembleHighKmer();

    this->mergeGraph();
    this->highKmerMakeBruijnGraph();
    do {
        numDelete = this->highKmerCutBranch();
        this->highKmerConcatinateNodes();
    } while (numDelete > 0);

    do {
        numDelete = this->highKmerCrushBubble();
        this->highKmerConcatinateNodes();
    } while (numDelete > 0);

}



//////////////////////////////////////////////////////////////////////////////////////
// De Bruijn assemble with low value kmer
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
inline void GapCloseDBG<M, N>::assembleLowKmer(void)
{
    unsigned long long numDelete = 0;

    this->lowKmerCount();
    this->lowKmerMakeBruijnGraph();
    do {
        numDelete = this->lowKmerCutBranch();
        this->lowKmerConcatinateNodes();
    } while (numDelete > 0);
}


//////////////////////////////////////////////////////////////////////////////////////
// De Bruijn assemble with high value kmer
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
inline void GapCloseDBG<M, N>::assembleHighKmer(void)
{
    unsigned long long numDelete = 0;

    this->highKmerCount();
    this->highKmerMakeBruijnGraph();
    do {
        numDelete = this->highKmerCutBranch();
        this->highKmerConcatinateNodes();
    } while (numDelete > 0);
}



//////////////////////////////////////////////////////////////////////////////////////
// check created contig with De Bruijn assemble overlapped with gap edge seq
// count the number of missmatch with gap edge seq (this program means headSeq and tailSeq)
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
template <typename KMER> 
bool GapCloseDBG<M, N>::closedGap(BruijnGraph<KMER> &graph)
{
    typedef std::pair<long, double> LENGTHandRATIO;
    long LENGTHandRATIO::* overlapLength = &LENGTHandRATIO::first;
    double LENGTHandRATIO::* missRate = &LENGTHandRATIO::second;
    long maxLeftOverlap = 0, maxRightOverlap = 0;
    long contigLength = 0;
    double minRate = 1.0;
    const Straight *maxStraight = NULL;


    LENGTHandRATIO left, right;
    for (auto it = graph.getLeftStraightTableBegin(); it != graph.getLeftStraightTableEnd(); ++it) {
        if (it->second.coverage == UINT16_MAX) continue;

        const Straight &straight = it->second;
        contigLength = straight.length + maxKmer - 1;
        left = calcMissmatchLeftEdgeSeq(straight);
        if (left.*overlapLength == 0) continue;
        right = calcMissmatchRightEdgeSeq(straight);
        if (right.*overlapLength == 0) continue;

		double bothMissRate = (right.*missRate * right.*overlapLength + left.*missRate * left.*overlapLength) / (right.*overlapLength + left.*overlapLength);

        if (bothMissRate < minRate || (bothMissRate == minRate && (right.*overlapLength + left.*overlapLength) > (maxRightOverlap + maxLeftOverlap))) {
            maxLeftOverlap = left.*overlapLength;
            maxRightOverlap = right.*overlapLength;
            minRate = bothMissRate;
            maxStraight = &straight;
        }
    }

    if (maxLeftOverlap == 0) return false;
    
    contigLength = maxStraight->length + maxKmer - 1;
    this->seqLength = contigLength - maxLeftOverlap - maxRightOverlap;
    if (this->seqLength > 0) {
        for (int i = 0; i < this->seqLength; ++i) {
            seq.emplace_back(maxStraight->get(contigLength - 1 - maxLeftOverlap - i));
        }
    } else {
        unsigned long long tmpLength = -seqLength;
        if (tmpLength > gapInfo->headSeq.size() || tmpLength > gapInfo->tailSeq.size()) return false;
        for (int i = 0; i < -seqLength; ++i) {
            if (gapInfo->headSeq[gapInfo->headSeq.size() - 1 - i] != gapInfo->tailSeq[-seqLength - 1 - i]) return false;
        }
        --seqLength;
    }
    return true;

}


template <unsigned M, unsigned N>
template <typename KMER> 
bool GapCloseDBG<M, N>::closedGapPartial(BruijnGraph<KMER> &graph)
{
	const int MIN_GAP_LENGTH = 10;

    typedef std::pair<long, double> LENGTHandRATIO;
    long LENGTHandRATIO::* overlapLength = &LENGTHandRATIO::first;
    double LENGTHandRATIO::* missRate = &LENGTHandRATIO::second;
    long maxLeftOverlap = 0, maxRightOverlap = 0;
    double minLeftRate = 1.0, minRightRate = 1.0;
    const Straight *maxLeftStraight = NULL;
    const Straight *maxRightStraight = NULL;


    LENGTHandRATIO left, right;
    for (auto it = graph.getLeftStraightTableBegin(); it != graph.getLeftStraightTableEnd(); ++it) {
        if (it->second.coverage == UINT16_MAX) continue;

        const Straight &straight = it->second;

        left = calcMissmatchLeftEdgeSeq(straight);
		if (left.*overlapLength > 0 && (left.*missRate < minLeftRate || (left.*missRate < minLeftRate && left.*overlapLength > maxLeftOverlap))) {
			maxLeftOverlap = left.*overlapLength;
			minLeftRate = left.*missRate;
			maxLeftStraight = &straight;
			continue;
		}

        right = calcMissmatchRightEdgeSeq(straight);
		if (right.*overlapLength > 0 && (right.*missRate < minRightRate || (right.*missRate < minRightRate && right.*overlapLength > maxRightOverlap))) {
			maxRightOverlap = right.*overlapLength;
			minRightRate = right.*missRate;
			maxRightStraight = &straight;
		}
    }

    if (maxLeftOverlap == 0 && maxRightOverlap == 0) return false;
    
    int leftContigLength = maxLeftStraight != NULL ? maxLeftStraight->length + maxKmer - 1 : 0;
    int rightContigLength = maxRightStraight != NULL ? maxRightStraight->length + maxKmer - 1 : 0;
    this->seqLength = (leftContigLength - maxLeftOverlap) + (rightContigLength - maxRightOverlap);

	int gapLength = MIN_GAP_LENGTH;
	if ((gapInfo->end - gapInfo->start) - this->seqLength > MIN_GAP_LENGTH)
		gapLength = (gapInfo->end - gapInfo->start) - this->seqLength;
	this->seqLength += gapLength;
	this->remainedGap = gapLength;

	for (int i = 0; i < leftContigLength - maxLeftOverlap; ++i)
		seq.emplace_back(maxLeftStraight->get(leftContigLength - 1 - maxLeftOverlap - i));
	for (int i = 0; i < gapLength; ++i)
		seq.emplace_back(4);
	for (int i = 0; i < rightContigLength - maxRightOverlap; ++i)
		seq.emplace_back(maxRightStraight->get(i));

    return true;
}



//////////////////////////////////////////////////////////////////////////////////////
// count the number of missmatch with Leftside edge seq (headSeq)
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
std::pair<long, double> GapCloseDBG<M, N>::calcMissmatchLeftEdgeSeq(const Straight &straight)
{
    long overlapLength = 0;
    double missRate = 1.0;
    const long contigLength = straight.length + maxKmer - 1;

    for (long length = contigLength; length >= minOverlap; --length) {
        double missTolerence = length * maxMissRate + 0.5;
        unsigned long long numMiss = 0;
        long i = 0;
        unsigned long long headSeqSize = gapInfo->headSeq.size();
        for (; i < minOverlap; ++i) {
            if (gapInfo->headSeq[headSeqSize - 1 - i] != straight.get(contigLength - length + i)) {
                ++numMiss;
                if (numMiss > missTolerence) break;
            }
        }
        if (i < minOverlap) continue;

        const long maxOverlap = std::min(static_cast<unsigned long long>(length), headSeqSize);
        missTolerence = maxOverlap * maxMissRate + 0.5;
        for (; i < maxOverlap; ++i) {
            if (gapInfo->headSeq[headSeqSize - 1 - i] != straight.get(contigLength - length + i)) {
                ++numMiss;
                if (numMiss > missTolerence) break;
            }
        }
        if (numMiss <= missTolerence && static_cast<double>(numMiss) / i < missRate) {
            missRate = static_cast<double>(numMiss) / i;
            overlapLength = length;
        }
    }
    return std::make_pair(overlapLength, missRate);
}



//////////////////////////////////////////////////////////////////////////////////////
// count the number of missmatch with Righttside edge seq (tailSeq)
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
std::pair<long, double> GapCloseDBG<M, N>::calcMissmatchRightEdgeSeq(const Straight &straight)
{
    long overlapLength = 0;
    double missRate = 1.0;
    const long contigLength = straight.length + maxKmer - 1;

    for (long length = contigLength; length >= minOverlap; --length) {
        double missTolerence = length * maxMissRate + 0.5;
        unsigned long long numMiss = 0;
        long i = 0;
        for (; i < minOverlap; ++i) {
            if (gapInfo->tailSeq[i] != straight.get(length - 1 - i)) {
                ++numMiss;
                if (numMiss > missTolerence) break;
            }
        }
        if (i < minOverlap) continue;

        const long maxOverlap = std::min(static_cast<unsigned long long>(length), static_cast<unsigned long long>(gapInfo->tailSeq.size()));
        missTolerence = maxOverlap * maxMissRate + 0.5;
        for (; i < maxOverlap; ++i) {
            if (gapInfo->tailSeq[i] != straight.get(length - 1 - i)) {
                ++numMiss;
                if (numMiss > missTolerence) break;
            }
        }
        if (numMiss <= missTolerence && static_cast<double>(numMiss) / i < missRate) {
            missRate = static_cast<double>(numMiss) / i;
            overlapLength = length;
        }
    }
    return std::make_pair(overlapLength, missRate);
}


//////////////////////////////////////////////////////////////////////////////////////
// save contig
//////////////////////////////////////////////////////////////////////////////////////
template <unsigned M, unsigned N>
void GapCloseDBG<M, N>::saveContig(FILE *contigFP) const
{
    this->save(contigFP, Int2Type<maxType>());
}



template <unsigned M, unsigned N>
template <typename KMER>
bool GapCloseDBG<M, N>::extendEdgeEntity(BruijnGraph<KMER> &graph)
{
    typedef std::pair<long, double> LENGTHandRATIO;
    long LENGTHandRATIO::* overlapLength = &LENGTHandRATIO::first;
    double LENGTHandRATIO::* missRate = &LENGTHandRATIO::second;
    long contigLength = 0;
    double minRate = 1.0;
    long minRateOverlap = 0;
    LENGTHandRATIO edge;
    const Straight *maxStraight = NULL;


    for (auto itr = graph.getLeftStraightTableBegin(); itr != graph.getLeftStraightTableEnd(); ++itr) {
        if (itr->second.coverage == UINT16_MAX) continue;

        const Straight &straight = itr->second;
        contigLength = straight.length + maxKmer - 1;

        if (gapInfo->start != 0) {
            edge = calcMissmatchLeftEdgeSeq(straight);
        } else {
            edge = calcMissmatchRightEdgeSeq(straight);
        }

        if (edge.*missRate < minRate) {
            minRateOverlap = edge.*overlapLength;
            minRate = edge.*missRate;
            maxStraight = &straight;
        }
    }

    if (minRateOverlap == 0) return false;

    contigLength = maxStraight->length + maxKmer - 1;
    this->seqLength = contigLength - minRateOverlap;

    if (this->seqLength <= 0) return false;

    if (gapInfo->start != 0) {
        for (long i = 0; i < this->seqLength; ++i) {
            seq.emplace_back(maxStraight->get(contigLength - 1 - minRateOverlap - i));
        }
    } else {
        for (long i = 0; i < this->seqLength; ++i) {
            seq.emplace_back(maxStraight->get(contigLength - 1 - i));
        }
    }
    return true;
}





#endif


