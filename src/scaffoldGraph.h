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

#ifndef SCAFFOLD_GRAPH_H
#define SCAFFOLD_GRAPH_H

#include "seqlib.h"
#include "mapper.h"
#include <unordered_map>
#include <vector>
#include <array>


class ScaffoldGraph
{
protected:
    /*
     * define inner classes using Scaffold class menber function
     */
    struct ScaffoldPart
    {
        long id;
        long start;
        long end;

        ScaffoldPart(): id(0), start(0), end(0) {}
        ScaffoldPart(long a, long b, long c): id(a), start(b), end(c) {}
        ScaffoldPart(const ScaffoldPart &) = default;
        ScaffoldPart &operator=(const ScaffoldPart &) = default;
        ~ScaffoldPart() = default;

        //added by ouchi
        bool operator==(const ScaffoldPart &a) const
        {
            if (id == a.id) {
                return start == a.start && end == a.end;
            } else {
                return false;
            }
        }
        bool operator!=(const ScaffoldPart &a) const
        {
            return !(*this == a);
        }
        //

    };

    struct GraphEdge
    {
        char direction;
        char state;
        long end;
        long length;
        long numLink;
        long score;
        std::vector<long> breakdown;

        GraphEdge(): direction(0), state(0), end(0), length(0), numLink(0), score(0), breakdown() {}
        GraphEdge(const GraphEdge &) = default;
        GraphEdge &operator=(const GraphEdge &) = default;
        ~GraphEdge() = default;

        bool isForward(void) const {return direction > 0; }

        bool operator<(const GraphEdge &a) const
        {
            if (direction != a.direction)
                return direction < a.direction;
            else
                return end < a.end;
        }
    };

//added by ouchi
    struct MPLink
    {
        char direction;
        long end;
        long offset1;
        long offset2;

        MPLink(): direction(0), end(0), offset1(0), offset2(0) {}
        MPLink(const MPLink &) = default;
        MPLink &operator=(const MPLink &) = default;
        ~MPLink() = default;

        bool operator<(const MPLink &a) const
        {
            if (end == a.end) {
                if (offset1 == a.offset1)
                    return offset2 < a.offset2;
                else
                    return offset1 < a.offset1;
            } else {
                return end < a.end;
            }
        }
    };
//

//added by ouchi
    struct LongLink
    {
        char direction;
        long end;
        long offset1;
        long offset2;
        long gap;

        LongLink(): direction(0), end(0), offset1(0), offset2(0), gap(0) {}
        LongLink(const LongLink &) = default;
        LongLink &operator=(const LongLink &) = default;
        ~LongLink() = default;

        bool operator<(const LongLink &a) const
        {
            if (end == a.end) {
                if (offset1 == a.offset1)
                    return offset2 < a.offset2;
                else
                    return offset1 < a.offset1;
            } else {
                return end < a.end;
            }
        }
    };
//

//added by ouchi
    struct HiCLink
    {
        long end;
        long offset1;
        long offset2;

        HiCLink(): end(0), offset1(0), offset2(0) {}
        HiCLink(const HiCLink &) = default;
        HiCLink &operator=(const HiCLink &) = default;
        ~HiCLink() = default;

        bool operator<(const HiCLink &a) const
        {
            if (end == a.end) {
                if (offset1 == a.offset1)
                    return offset2 < a.offset2;
                else
                    return offset1 < a.offset1;
            } else {
                return end < a.end;
            }
        }
    };
//

    struct GraphNode
    {
        bool isHomo;
        char state;
        long length;
        long numEdge;
        std::vector<GraphEdge> edge;
        std::vector<std::vector<MPLink> > MPLinks; //added by ouchi
        //std::vector<std::vector<MPLink> > MPINLinks; //added by ouchi
        std::vector<LongLink> LongLinks; //added by ouchi
        std::vector<HiCLink> HiCLinks; //added by ouchi
        long numContig;
        std::vector<ScaffoldPart> contig;
        std::unordered_map<int, unsigned> numMappedTag;
		long oppositeBubbleNodeID;

        GraphNode(): isHomo(false), state(0), length(0), numEdge(0), edge(), MPLinks(), LongLinks(), HiCLinks(), numContig(0), contig(), numMappedTag(), oppositeBubbleNodeID(0) {} //edited by ouchi (MPLinks(), HiCLinks())
        GraphNode(const GraphNode &) = default;
        GraphNode &operator=(const GraphNode &) = default;
        ~GraphNode() = default;

        //added by ouchi
        bool operator==(const GraphNode &a) const
        {
            if (numContig != a.numContig) {
                return false;
            } else {
                long i;
                char direction = 0;
                for (i = 0; i < numContig; ++i) {
                    if (contig[i].id != a.contig[i].id) {
                        if (direction == 1) break;
                        direction = -1;
                    }
                    if (contig[i].id != -(a.contig[numContig-1-i].id)) {
                        if (direction == -1) break;
                        direction = 1;
                    }
                }
                return i == numContig;
            }
        }
        bool operator!=(const GraphNode &a) const
        {
            return !(*this == a);
        }
        //
    };

    struct GraphLink
    {
        long id1;
        long id2;
        long offset1;
        long offset2;
        long gap;

        GraphLink(): id1(0), id2(0), offset1(0), offset2(0), gap(0) {}
        ~GraphLink() = default;
        void clearValue(void)
        {
            id1 = 0;
            id2 = 0;
            offset1 = 0;
            offset2 = 0;
            gap = 0;
        }
        bool operator<(const GraphLink &a) const
        {
            if (id1 != a.id1)
                return (id1 - a.id1) < 0;
            else if (id2 != a.id2)
                return (id2 - a.id2) < 0;
            else
                return (gap - a.gap) < 0;
        }
    };

    struct GraphLinkPoolIndex
    {
        unsigned long index;
        long numLink;

        GraphLinkPoolIndex(): index(0), numLink(0) {}
        GraphLinkPoolIndex(unsigned long idx): index(idx), numLink(0) {}
        ~GraphLinkPoolIndex() = default;
    };

    struct GraphLinkPoolIndexGreater
    {
        bool operator() (const GraphLinkPoolIndex& link1, const GraphLinkPoolIndex& link2) const
        { return link1.numLink > link2.numLink; }
    };

    struct GraphLayout
    {
        long id;
        long start;
        long end;
        long distance;
        long numLink;
        char h; //added by ouchi
        long score; //added by ouchi for Platanus-allee v2.3.2

        GraphLayout(): id(0), start(0), end(0), distance(0), numLink(0), h(0), score(0) {} //edited by ouchi(h(0), score(0))
        ~GraphLayout() = default;

        bool operator<(const GraphLayout &a) const
        {
            return (start - a.end - (a.start - end)) < 0;
        }
    };

    struct Overlap
    {
        int id1;
        int id2;
        int length;

        Overlap(): id1(0), id2(0), length(0) {}
        ~Overlap() = default;
    };

    struct HaplotypeBlockInfo
    {
        bool isHetero;
		std::vector<char> altSeq;
		long start;
		long end;

        HaplotypeBlockInfo(): isHetero(false), altSeq(0), start(0), end(0) {}
        HaplotypeBlockInfo &operator=(const HaplotypeBlockInfo &) = default;
        ~HaplotypeBlockInfo() = default;
    };

    struct NodeIDWithGap
	{
		long id;
		long gap;

        NodeIDWithGap(): id(0), gap(0) {}
        NodeIDWithGap(long i, long n): id(i), gap(n) {}
	};

    struct GraphPathGapped
    {
		long selfID;
		std::vector<NodeIDWithGap> node;
		unsigned sumLink;

        GraphPathGapped(): selfID(), node(), sumLink(0) {}
        GraphPathGapped(long ID, unsigned long size): selfID(ID), node(size) {}
        GraphPathGapped(long ID, unsigned long size, unsigned sum): selfID(ID), node(size), sumLink(sum) {}
	};

    struct GraphPathGappedSelfIDLess
    {
        bool operator() (const GraphPathGapped& path1, const GraphPathGapped& path2) const
        { return path1.selfID < path2.selfID; }
    };

	struct ResultSeq
	{
		std::string seq;
		std::string name;
		std::string component;
		bool bubbleFlag;
		bool redundantFlag;

        ResultSeq(): seq(), name(), component(), bubbleFlag(false), redundantFlag(false) {}
	};


    // end definition inner classes

    static const unsigned TABLE_DIVID;
    static const double MAX_DIFF_RATE;
    static const double EDGE_EXPECTED_RATE_TH;
    static const double EDGE_EXPECTED_RATE_UPPER_TH;
    static const double CHECK_USING_LONGER_LIB_TH;
    static const unsigned SC_REP;
    static const unsigned SC_INC;
    static const unsigned SC_DEL;
    static const unsigned SC_RR; //added by ouchi
    static const unsigned SC_RL; //added by ouchi
    static const double MAX_HOMO_RATE;
    static const double MAX_HETERO_RATE;
    static const double MAX_OVERLAP_IDENTITY_DIFF;
	static const long MIN_NUM_MAPPED_TAG;

    long seedLength;
    long keyLength;
    long minOverlap;
    long hashOverlap;
    long indexLength;
    long minLink;
    long tolerence;
    long minTolerenceFactor;
    long cutoffLength;
    long genomeSize;
    long numContig;
    long numNode;
    long maxFragmentLengthOfTag;
    double averageCoverage;
    double bubbleThreshold;
    FILE* contigFP;
    FILE* bubbleFP;
    FILE* bubbleOpositeFP;
    FILE* deletedHeteroFP;
    FILE* remainingHeteroFP;
    FILE* overlapFP;
    std::vector<FILE*> graphLinkFP;
    std::vector<platanus::SEQ> contig;
    std::vector<std::string> contigName;
	std::unordered_map<std::string, unsigned> contigNameIndex;
    std::unordered_set<std::pair<int, int>, platanus::PairHash, platanus::PairEqual> contigUnlinkSet;
    std::vector<std::vector<SeqLib> > *allLibraryMT;
    std::vector<SeqLib> *longReadLibraryMT;
    std::vector<SeqLib> *tagLibraryMT;
    std::vector<SeqLib> *HiCLibraryMT; //added by ouchi
   	long targetLibraryIndex;
    std::vector<unsigned short> coverage;
	std::vector<char> contigState;
    std::vector<char> seqPool;
    std::vector<GraphNode> node;
    std::vector<long> numBubble;
    std::vector<std::vector<platanus::Position> > contigPositionInScaffold; //edited by ouchi
    std::vector<std::unordered_map<std::pair<int, int>, Overlap, platanus::PairHash, platanus::PairEqual> > overlapTable;
	std::unordered_map<long, long> bubblePairNodeIDMap;
	std::vector<ResultSeq> resultSeq;
	std::vector<std::unordered_map<std::pair<int, int>, unsigned char, platanus::PairHash, platanus::PairEqual> > contigTagCounterUchar;
	std::vector<std::unordered_map<std::pair<int, int>, unsigned int, platanus::PairHash, platanus::PairEqual> > contigTagCounterUint;

    void destroyGraph(void);
    long getOverlap(long id1, long id2);
    long getShortOverlap(long id1, long id2) const;
    long getScaffoldOverlap(long id1, long id2);
    bool checkDeleteEdge(const GraphEdge &edge1, const GraphEdge &edge2, const GraphNode &node1, const GraphNode &node2);
    long deleteErroneousEdge(const long numThread);
	long deleteErroneousEdgeNumLinkRate(const long numThread);
    void deleteEdges(std::vector<long> &ids);
    void remake(const long numNewNode, const long numContigPoolSize, FILE *scaffoldFP);
    double calcExpectedLink(const double coverage, const double link1, const double link2, const double g) const;
    virtual double calcNodeCoverage(const GraphNode &node);
    unsigned long long crushBubble(const double bubbleThreshold, const double averageCoverage, const long numThread);
    void layoutNodes(GraphNode *newNode, std::vector<GraphLayout> &ret, std::vector<GraphLayout> &work);
    void layout2seq(const std::vector<GraphLayout> &lerfOverlap, const long startPoint, const long numLeftOverlap, std::vector<char> &ret);
    long alignScaffold(const std::vector<char> &scaffold1, const std::vector<char> &scaffold2, std::vector<long> &work, const long scoreThreshold) const;
    double layoutAverageCoverage(const std::vector<GraphLayout> &leftOverlap, const long startPoint, const long leftOverlapSize) const;
    long getSimilarOverlap(const long id1, const long id2);
    void calcLinkAndWriteGraphLinkFile(const std::vector<GraphLink>& links, const GraphLinkPoolIndex& index, const long libraryIndex, const bool multiThreadFlag);
    long estimateGapSize(const std::vector<GraphLink>& links, const unsigned index, const unsigned size) const;
    long estimateGapSizeAverage(const std::vector<GraphLink>& links, const unsigned index, const unsigned size) const;
    long calcNumPossiblePosition(const long length1, const long length2, const long distance, const long insSize) const;
    long calcNumPossiblePositionNode(const GraphNode& node1, const GraphNode& node2, const long distance, const long insSize) const;
    long calcNumPossiblePositionNodeTemp(const long node1Start, const long node1End, const long node1Length, const GraphNode& node2, const long distance, const long insSize) const;
    double calcExpectedLinkNode(const GraphNode& node1, const GraphNode& node2, const long distance) const;
    double calcExpectedLinkNodeTemp(const long node1Start, const long node1End, const long node1Length, const GraphNode& node2, const long distance) const;
	void storeGraphLinkFromMappedPair(std::vector<GraphLink> &graphLinkPool, long numThread);
	long getCommonTagBetweenNodePair(const long nodeID1, const long nodeID2);
	long getNumLinkFromIDPair(long leftNodeID, long rightNodeID);
	long remakeGraphAccordingToGappedPath(std::vector<GraphPathGapped> &pathBuffer);
	long remakeGraphAccordingToGappedPathPair(std::vector<GraphPathGapped> &pathBuffer);
	long writeAndMarkGappedNodes(const std::vector<NodeIDWithGap> &gappedPath, FILE *storeFP);
	void writeSingletonNode(long &numNewNode, long &newContigPoolSize, FILE *storeFP);

    template <typename T>
    inline T id2Index(T id) { return std::abs(id) - 1; }

    unsigned decideTableID(const unsigned long long key)
    {
        static const unsigned long long ander = TABLE_DIVID - 1;
        return key & ander;
    }

	long calcGapSizeFromPositions(const long start1, const long end1, const long start2, const long end2)
	{
		return abs(std::min(end1, end2) - std::min(start1, start2)) - ((end1 - start1) + (end2 - start2));
	}


public:
    ScaffoldGraph();
    ~ScaffoldGraph();
    ScaffoldGraph(const ScaffoldGraph&) = delete;
    ScaffoldGraph& operator=(const ScaffoldGraph&) = delete;

    void saveOverlap(const Mapper &map, const long hashOverlapValue, const long cutoffLengthi, const long numThread);
    virtual void makeGraph(const long numThread);
    virtual void calcLink(const long libraryIndex, const long linkThreshold, const long numThread);
    virtual void cutAndPrintSeq(const long minSeqLength, const unsigned long long readLength, const std::string &outFilename, const std::string &componentFilename);
    virtual void detectRepeat(const double averageCoverage);
    virtual void deleteRepeatEdge(void);
    void deleteThinEdge(void);
	void deleteThinEdgeConstant(const long linkThreshold);
    void deleteErroneousEdgeIterative(const long numThread);
	void deleteErroneousEdgeNumLinkRateIterative(const long numThread);
    void split(void);
    virtual void makeScaffold(void);
    long estimateLink(void);
    void crushBubbleIterative(const double bubleThreshold, const double averageCoverage, const long numThread);
    virtual unsigned long long crushHeteroBubble(const double averageCoverage);
    void initScaffolding(std::vector<unsigned short> &cov, Mapper &mapper, const double ave, const double bubble);
    void countBubble(const platanus::Contig &bubble, const std::vector<platanus::Position> &bubblePosition);
    void classifyNode(void);
    virtual long deleteHeteroEdge(void);
    void removeHeteroOverlap(void);
    void printScaffoldBubble(const std::string &outFilename);
    virtual void insertSizeDistribution(std::vector<SeqLib>& library, std::vector<long>& distribution, const long numThread);
    void scaffoldLengthList(std::vector<long>& list);
	void clearEdges();
	void node2seq(const GraphNode &node, std::vector<char> &ret);
	void writeNodeSeq(const GraphNode &targetNode, FILE *fp);
	long readNodeSeq(std::vector<char> &nodeSeq, std::vector<long> &contigID, FILE *fp);
	void setContigName(const std::vector<std::string> &name);
	void setContigNameIndex(const std::unordered_map<std::string, unsigned> nameIndex);
	void setContigUnlinkFlags(const std::string fileName);
	void dumpAllEdges(const std::string &outputFilename);
	void outputFastg(const std::string &outputFilename);
	void countMappedTagForEachContig(long numThread);
	void countMappedTagForEachScaffold(long numThread);
	virtual void loadResultSeq(const long minSeqLength, const unsigned long long readLength);
	void printBinSeqConstLineLength(const std::vector<char> &seq, std::ofstream &out);
	void printResultSeq(const std::string &outFilename);
	void mapBubbleToResultSeq();

    long getTolerence(void) const { return tolerence; }
    long getCutoffLength(void) const { return cutoffLength; }
    long getNumNode(void) const { return numNode; }
    long getMinLink(void) const { return minLink; }
    void setSeedLength(const long len) { seedLength = len; }
    void setMaxFragmentLengthOfTag(const long len) { maxFragmentLengthOfTag = len; }
    void setMinOverlap(const long olp) { minOverlap = olp; }
    void setAllLibraryMT(std::vector<std::vector<SeqLib> > *input) { allLibraryMT = input; }
    void setLongReadLibraryMT(std::vector<SeqLib> *input) { longReadLibraryMT = input; }
    void setTagLibraryMT(std::vector<SeqLib> *input) { tagLibraryMT = input; }
//    void setTargetLibraryIndex(const long index) { targetLibraryIndex = index; averageCoverage = (*allLibraryMT)[index][0].getAverageCoverage(); }
    void setTargetLibraryIndex(const long index) { targetLibraryIndex = index; if (allLibraryMT != NULL && index < static_cast<long>((*allLibraryMT).size())) averageCoverage = (*allLibraryMT)[index][0].getAverageCoverage(); }
    void setHiCLibraryMT(std::vector<SeqLib> *input) {HiCLibraryMT = input; } //added by ouchi
    void clearSeqLib() { (*allLibraryMT)[targetLibraryIndex].clear(); averageCoverage = 0.0; }
    void clearResultSeq() { resultSeq.clear(); }
    void setTolerence(const long tol) { tolerence = tol; }
    void setMinTolerenceFactor(const long fac) { minTolerenceFactor = fac; }
    void setCutoffLength(const long cut) { cutoffLength = cut; }
    void setMinLink(const unsigned long num) { minLink = num; }
    void setAverageCoverage(const double cov) { averageCoverage = cov; }
	void getLinkedNode(const long sourceNodeIndex, const char targetDirection, std::vector<NodeIDWithGap> &nodeIDBuffer);
	void getLinkedUniqueNode(const long sourceNodeIndex, const char targetDirection, std::vector<NodeIDWithGap> &nodeIDBuffer);
	virtual void updateInsertLengthFP(std::vector<SeqLib>& lib, const long numThread);
	virtual long joinUnambiguousNodePairGapped(const long numThread);
	void joinUnambiguousNodePairGappedIterative(const long numThread);
	void deleteEdgeFromShortNode(const long lengthThreshold);
	void deleteEdgeFromShortAndAbnormalCoverageNode(const long minLength, const double minCoverage, const double maxCoverage);

    template <typename T>
    inline T sign(T val) { return (T(0) < val) - (T(0) > val); }

    template <typename T>
    std::pair<long, long> trimCommonEdgePart(T first1, T last1, T first2, T last2)
	{
		std::pair<long, long> trimmedLength;

		long halfSize = std::min(last1 - first1, last2 - first2) / 2;
		for (trimmedLength.first = 0; trimmedLength.first < halfSize; ++trimmedLength.first) {
			if (*(first1 + trimmedLength.first) != *(first2 + trimmedLength.first))
				break;
		}

		halfSize = std::min(last1 - first1, last2 - first2) / 2;
		for (trimmedLength.second = 0; trimmedLength.second < halfSize; ++trimmedLength.second) {
			if (*(last1 - trimmedLength.second - 1) != *(last2 - trimmedLength.second - 1))
				break;
		}

		return trimmedLength;
	}

    template <typename T>
	void mergeAndClearMultiThreadedVector(std::vector<std::vector<T> > &threadVector, std::vector<T> &mergedVector)
	{
		mergedVector.clear();
		for (auto itr = threadVector.begin() ; itr != threadVector.end(); ++itr) {
			std::copy(itr->begin(), itr->end(), std::back_inserter(mergedVector));
			itr->clear();
		}
	}

	template <typename T>
	void reverseComplement(T &seq)
	{
		std::reverse(seq.begin(), seq.end());
		for (auto itr = seq.begin(); itr != seq.end(); ++itr) {
			if (*itr < 4)
				*itr ^= 0x3;
		}
	}

	template <typename T>
	void printBinSeqConstLineLength(const T &seq, std::ofstream &out)
	{
		unsigned long j;
		for (j = 0; j < seq.size(); ++j) {
			out << platanus::Bin2Char(seq[j]);
			if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
				out.put('\n');
		}
		if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
			out.put('\n');
	}
};


//////////////////////////////////////////////////////////////////////////////////////
// calc expected number of links
//////////////////////////////////////////////////////////////////////////////////////
inline double ScaffoldGraph::calcExpectedLink(const double coverage, const double link1, const double link2, const double g) const
{
    const double averageIns = (double)(*allLibraryMT)[targetLibraryIndex][0].getAverageInsSize();
    const double sdIns = (double)(*allLibraryMT)[targetLibraryIndex][0].getSDInsSize();
    const double average = (double)(*allLibraryMT)[targetLibraryIndex][0].getAverageLength();
    double numLink = 0;

    numLink += (link1 + g - averageIns + link2) * erf((link1 + g - averageIns + link2) / (M_SQRT2 * sdIns));
    numLink += (M_SQRT2 * sdIns / sqrt(M_PI)) * exp(-pow(((link1 + g - averageIns + link2) / (M_SQRT2 * sdIns)), 2.0));
    numLink -= (average + g - averageIns + link2) * erf((average + g - averageIns + link2) / (M_SQRT2 * sdIns));
    numLink -= (M_SQRT2 * sdIns / sqrt(M_PI)) * exp(-pow(((average + g - averageIns + link2) / (M_SQRT2 * sdIns)), 2.0));

    numLink -= (link1 + g - averageIns + average) * erf((link1 + g - averageIns + average) / (M_SQRT2 * sdIns));
    numLink -= (M_SQRT2 * sdIns / sqrt(M_PI)) * exp(-pow(((link1 + g - averageIns + average) / (M_SQRT2 * sdIns)), 2.0));
    numLink += (average + g - averageIns + average) * erf((average + g - averageIns + average) / (M_SQRT2 * sdIns));
    numLink += (M_SQRT2 * sdIns / sqrt(M_PI)) * exp(-pow(((average + g - averageIns + average) / (M_SQRT2 * sdIns)), 2.0));

    numLink *= coverage / (4.0 * average);

    return numLink;

}


//////////////////////////////////////////////////////////////////////////////////////
// calc expected number of links between two nodes
//////////////////////////////////////////////////////////////////////////////////////
inline double ScaffoldGraph::calcExpectedLinkNode(const GraphNode& node1, const GraphNode& node2, const long distance) const
{
    if (node1.contig.empty() || node2.contig.empty()) return 1.0;

    double expected = 0;
    long node1Start = node1.contig[0].start;
    long node1End   = node1.contig[0].end;
    for (unsigned idx1 = 1; idx1 < node1.contig.size(); ++idx1) {
        if (node1.contig[idx1].start <= node1End) {
            node1End = node1.contig[idx1].end;
            continue;
        }
        expected += calcExpectedLinkNodeTemp(node1Start, node1End, node1.length, node2, distance);
        node1Start = node1.contig[idx1].start;
        node1End   = node1.contig[idx1].end;
    }
    return expected + calcExpectedLinkNodeTemp(node1Start, node1End, node1.length, node2, distance);
}


//////////////////////////////////////////////////////////////////////////////////////
// calc expected number of links between two nodes (utility)
//////////////////////////////////////////////////////////////////////////////////////
inline double ScaffoldGraph::calcExpectedLinkNodeTemp(const long node1Start, const long node1End, const long node1Length, const GraphNode& node2, const long distance) const
{
    double expected = 0;
    long node2Start = node2.contig[0].start;
    long node2End   = node2.contig[0].end;
    for (unsigned idx2 = 1; idx2 < node2.contig.size(); ++idx2) {
        if (node2.contig[idx2].start < node2End) {
            node2End = node2.contig[idx2].end;
            continue;
        }
        expected += calcExpectedLink(averageCoverage,
                node1End - node1Start + 1, node2End - node1Start + 1,
                distance + node1Length - node1End + node2Start);
        node2Start = node2.contig[idx2].start;
        node2End   = node2.contig[idx2].end;
    }
    expected += calcExpectedLink(averageCoverage,
            node1End - node1Start + 1, node2End - node1Start + 1,
            distance + node1Length - node1End + node2Start);
    return expected;
}


//////////////////////////////////////////////////////////////////////////////////////
// calc possible number of positions that reads can be mapped
//////////////////////////////////////////////////////////////////////////////////////
inline long ScaffoldGraph::calcNumPossiblePosition(const long length1, const long length2, const long distance, const long insSize) const
{
    const long minNodeLength = std::min(length1, length2);
    const long totalNodeLength = length1 + length2;
    const long readLength = (*allLibraryMT)[targetLibraryIndex][0].getAverageLength() / 2;

    long way = std::max(0l, insSize - (distance > 0 ? distance : 0) - readLength * 2 + 1);
    way = std::min(way, std::max(0l, minNodeLength + readLength + 1));
    way = std::min(way, std::max(0l, totalNodeLength + distance - insSize + 1));
    return way;
}


//////////////////////////////////////////////////////////////////////////////////////
// calc possible number of positions that reads can be mapped between two nodes
//////////////////////////////////////////////////////////////////////////////////////
inline long ScaffoldGraph::calcNumPossiblePositionNode(const GraphNode& node1, const GraphNode& node2, const long distance, const long insSize) const
{
    if (node1.contig.empty() || node2.contig.empty()) return 0;

    long way = 0;
    long node1Start = node1.contig[0].start;
    long node1End   = node1.contig[0].end;
    for (unsigned idx1 = 1; idx1 < node1.contig.size(); ++idx1) {
        if (node1.contig[idx1].start <= node1End) {
            node1End = node1.contig[idx1].end;
            continue;
        }
        way += calcNumPossiblePositionNodeTemp(node1Start, node1End, node1.length, node2, distance, insSize);
        node1Start = node1.contig[idx1].start;
        node1End   = node1.contig[idx1].end;
    }
    return way + calcNumPossiblePositionNodeTemp(node1Start, node1End, node1.length, node2, distance, insSize);
}


//////////////////////////////////////////////////////////////////////////////////////
// calc possible number of positions that reads can be mapped between two nodes (utility)
//////////////////////////////////////////////////////////////////////////////////////
inline long ScaffoldGraph::calcNumPossiblePositionNodeTemp(const long node1Start, const long node1End, const long node1Length, const GraphNode& node2, const long distance, const long insSize) const
{
    long way = 0;
    long node2Start = node2.contig[0].start;
    long node2End   = node2.contig[0].end;
    for (unsigned idx2 = 1; idx2 < node2.contig.size(); ++idx2) {
        if (node2.contig[idx2].start < node2End) {
            node2End = node2.contig[idx2].end;
            continue;
        }
        way += calcNumPossiblePosition(
                node1End - node1Start + 1, node2End - node2Start + 1,
                distance + node1Length - node1End + node2Start, insSize);
        node2Start = node2.contig[idx2].start;
        node2End   = node2.contig[idx2].end;
    }
    way += calcNumPossiblePosition(
            node1End - node1Start + 1, node2End - node2Start + 1,
            distance + node1Length - node1End + node2Start, insSize);
    return way;
}




#endif
