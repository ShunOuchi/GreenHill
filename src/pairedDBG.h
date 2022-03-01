/*
Copyright (C) 2022 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus-3D.

Platanus-3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus-3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus-3D; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef PAIRED_DBG_H
#define PAIRED_DBG_H
#define STATIC_VERSION

#include "seqlib.h"
#include "mapper.h"
#include "scaffoldGraph.h"
#include <unordered_map>
#include <vector>
#include <array>


class PairedDBG : public ScaffoldGraph
{
public:
    enum CROSS_RESOLUTION_MODE {LINK, TAG, SCORE, HIC}; //edited by ouchi
    enum DIVISION_MODE {SWITCH, GAP, MIS};

protected:

	struct MapInfoForGraph
	{
		platanus::Position position;
		long contigIndex;
		long score;
		long tStart; //added by ouchi
		long tEnd; //added by ouchi
		long qStart; //added by ouchi
		long qEnd; //added by ouchi

		struct PositionIDLessScoreGreater
		{
			bool operator()(const MapInfoForGraph &a, const MapInfoForGraph &b)
			{
				return (a.position.id != b.position.id) ? (a.position.id < b.position.id) : ((a.score != b.score) ? (a.score > b.score) : (a.contigIndex < b.contigIndex));
			}
		};

		struct PositionIDLess
		{
			bool operator()(const MapInfoForGraph &a, const MapInfoForGraph &b)
			{
				return (a.position.id < b.position.id);
			}
		};

//		bool operator <(const MapInfoForGraph &a) const { return (position < a.position); }
		MapInfoForGraph(): position(0, 0), contigIndex(0), score(0), tStart(0), tEnd(0), qStart(0), qEnd(0) {} //edited by ouchi
		MapInfoForGraph(const platanus::Position p, const long c, const long s): position(p), contigIndex(c), score(s),tStart(0), tEnd(0), qStart(0), qEnd(0) {} //edited by ouchi
		MapInfoForGraph(const platanus::Position p, const long c, const long s, const long ts, const long te, const long qs, const long qe): position(p), contigIndex(c), score(s),tStart(ts), tEnd(te), qStart(qs), qEnd(qe) {} //added by ouchi
		~MapInfoForGraph() {}
	};

    struct GraphLinkWithFlag : public GraphLink
    {
		bool overlapFlag;
		long score;

        GraphLinkWithFlag(): GraphLink(), overlapFlag(false), score(1) {}
        ~GraphLinkWithFlag() = default;
	};

    struct GraphLinkWithFlagPoolIndex : public GraphLinkPoolIndex
    {
		bool overlapFlag;
		long gap;

        GraphLinkWithFlagPoolIndex(): GraphLinkPoolIndex(), overlapFlag(false), gap(0) {}
        GraphLinkWithFlagPoolIndex(unsigned long idx): GraphLinkPoolIndex(idx), overlapFlag(false), gap(0) {}
        ~GraphLinkWithFlagPoolIndex() = default;
    };

    struct GraphLinkWithFlagPoolIndexGreater
    {
        bool operator() (const GraphLinkWithFlagPoolIndex& link1, const GraphLinkWithFlagPoolIndex& link2) const
        { return link1.numLink > link2.numLink; }
    };

    //added by ouchi
    struct ConsensusLink
    {
        long id1;
        long id2;
        char h1;
        char h2;

        ConsensusLink(): id1(0), id2(0), h1(0), h2(0) {}
        ~ConsensusLink() = default;
        void clearValue(void)
        {
            id1 = 0;
            id2 = 0;
            h1 = 0;
            h2 = 0;
        }
        bool operator<(const ConsensusLink &a) const
        {
            if (id1 != a.id1)
                return (id1 - a.id1) < 0;
            else
                return (id2 - a.id2) < 0;
        }
    };
    //

    //added by ouchi
    struct ConsensusLinkPoolIndex
    {
        unsigned long index;
        std::array<long, 2> numLinks;
        long numLink;

        ConsensusLinkPoolIndex(): index(0), numLinks(), numLink(0) {}
        ConsensusLinkPoolIndex(unsigned long idx): index(idx), numLinks(), numLink(0) {}
        ~ConsensusLinkPoolIndex() = default;
    };
    //

    //added by ouchi
    struct ConsensusLinkPoolIndexGreater
    {
        bool operator() (const ConsensusLinkPoolIndex& link1, const ConsensusLinkPoolIndex& link2) const
        { return link1.numLink > link2.numLink; }
    };
    //

    struct GraphPath
    {
		long selfID;
		std::vector<long> nodeID;
		unsigned sumLink;

        GraphPath(): selfID(), nodeID(), sumLink(0) {}
        GraphPath(long ID, unsigned long size): selfID(ID), nodeID(size) {}
        GraphPath(long ID, unsigned long size, unsigned sum): selfID(ID), nodeID(size), sumLink(sum) {}
	};

    struct GraphPathSelfIDLess
    {
        bool operator() (const GraphPath& path1, const GraphPath& path2) const
        { return path1.selfID < path2.selfID; }
    };

    struct NodeInfoForPathSearch
    {
		unsigned numVisit;
		long distance;
		long preNodeID;

        NodeInfoForPathSearch(): numVisit(0), distance(0) {}
        NodeInfoForPathSearch(unsigned n, long d, long p):  numVisit(n), distance(d), preNodeID(p) {}
	};

	struct ContigBubbleInfo
	{
		std::array<long, 2> joinedContigID;
		std::array<long, 2> oppositeContigID;

        ContigBubbleInfo() { joinedContigID.fill(0); oppositeContigID.fill(0); }
	};

    //added by ouchi
    struct ConsensusIDWithGap
    {
        long id;
        long gap;
        char h;
        long numLink;

        ConsensusIDWithGap(): id(0), gap(0), h(0), numLink(0) {}
        ConsensusIDWithGap(long i, long n, char c, long l): id(i), gap(h), h(c), numLink(l) {}
    };
    //

    //added by ouchi
    struct ConsensusPathGapped
    {
        long selfID;
        std::vector<ConsensusIDWithGap> consensusnode;

        ConsensusPathGapped(): selfID(), consensusnode() {}
        ConsensusPathGapped(long ID, unsigned long size): selfID(ID), consensusnode(size) {}
    };
    //

    //added by ouchi
    struct ConsensusEdge
    {
        char direction;
        long end;
        long length;
        std::array<long, 2> numLink;

        ConsensusEdge(): direction(0), end(0), length(0), numLink() {}
        ConsensusEdge(const ConsensusEdge &) = default;
        ConsensusEdge &operator=(const ConsensusEdge &) = default;
        ~ConsensusEdge() = default;

        bool operator<(const ConsensusEdge &a) const
        {
            if (direction != a.direction)
                return direction < a.direction;
            else
                return end < a.end;
        }
    };
    //

/*    //added by ouchi
    struct ConsensusPathEdge
    {
        char direction;
        long end;
        long length;
        std::array<long, 2> numLink;
        std::vector<std::vector<long> > path;
        std::vector<std::vector<char>> h;

        ConsensusPathEdge(): direction(0), end(0), length(0), numLink(), path(), h() {}
        ConsensusPathEdge(const ConsensusPathEdge &) = default;
        ConsensusPathEdge &operator=(const ConsensusPathEdge &) = default;
        ~ConsensusPathEdge() = default;
    };
    //
*/

    //added by ouchi
    struct HiCConsensusEdge
    {
        long id1;
        long id2;
        char h;
        long numLink;

        HiCConsensusEdge(): id1(0), id2(0), h(0), numLink(0) {}
        HiCConsensusEdge(const HiCConsensusEdge &) = default;
        HiCConsensusEdge &operator=(const HiCConsensusEdge &) = default;
        ~HiCConsensusEdge() = default;
    };
    //

    //added by ouchi
    struct ConsensusNode
    {
        char state;
        //char searchState;
        long length;
        long numEdge;
        //long numPathEdge;
        std::vector<ConsensusEdge> edge;
        //std::vector<ConsensusPathEdge> pathEdge;
        //std::vector<HiCLink> HiCLinks;
        std::array<long, 2> numNode;
        std::array<std::vector<ScaffoldPart>, 2> node;

        ConsensusNode(): state(0), length(0), numEdge(0), edge(), numNode(), node() {}
        ConsensusNode(const ConsensusNode &) = default;
        ConsensusNode &operator=(const ConsensusNode &) = default;
        ~ConsensusNode() = default;
    };
    //

    //added by ouchi
    struct ConsensusPart
    {
        long id;
        long start;
        long end;
        char h;

        ConsensusPart(): id(0), start(0), end(0), h(0) {}
        ConsensusPart(long a, long b, long c, char d): id(a), start(b), end(c), h(d) {}
        ConsensusPart(const ConsensusPart &) = default;
        ConsensusPart &operator=(const ConsensusPart &) = default;
        ~ConsensusPart() = default;

    };
    //

    //added by ouchi
    struct HiCNode
    {
        char state;
        long length;
        long numScaffold;
        std::vector<ConsensusPart> scaffold;

        HiCNode(): state(0), length(0), numScaffold(0), scaffold() {}
        HiCNode(const HiCNode &) = default;
        HiCNode &operator=(const HiCNode &) = default;
        ~HiCNode() = default;
    };
    //

    //added by ouchi
    /*struct HiCNodeGraph
    {
        std::vector<long> par;
        std::vector<long> sizes;
        std::vector<long> degrees;
        std::vector<std::vector<long> > scaffolds;

        HiCNodeGraph(): par(), sizes(), degrees(), scaffolds() {}
        HiCNodeGraph(const HiCNodeGraph &) = default;
        HiCNodeGraph &operator=(const HiCNodeGraph &) = default;
        ~HiCNodeGraph() = default;

        void initGraph(const long numScaffolds)
        {
            scaffolds.resize(numScaffolds);
            par.resize(numScaffolds);
            sizes.resize(numScaffolds, 2);
            degrees.resize(numScaffolds, 1);
            for (long i = 0; i < numScaffolds; ++i) {
                scaffolds[i].push_back(i+1);
                scaffolds[i].push_back(-(i+1));
                par[i] = i;
            }

        }
        long findParent(long id)
        {
            if (par[id] == id) {
                return id;
            } else {
                return par[id] = findParent(par[id]);
            }
        }
        std::vector<long> getScaffold(long id, char mode)
        {
            long p = findParent(id);
            std::vector<long> scaffold = scaffolds[p];
            if (mode == -1) return scaffold;
            if (degrees[id] == 1) {
                if (scaffold[0] == id) {
                    if (mode == 1) {
                        std::reverse(scaffold.begin(), scaffold.end());
                    }
                    return scaffold;
                } else if (scaffold.back() == id) {
                    if (mode == 0) {
                        std::reverse(scaffold.begin(), scaffold.end());
                    }
                    return scaffold;
                } else {
                    std::cerr << "failed to get Scaffold" << std::endl;
                }
            } else {
                long iditr;
                for (iditr = 0; iditr < scaffold.size(); ++iditr) {
                    if (scaffold[iditr] == id) break;
                }
                if (iditr == scaffold.size())    std::cerr << "failed to get Scaffold" << std::endl;
                std::vector<long> scaffold1, scaffold2;
                if (iditr % 2 == 0) {
                    scaffold2 = std::vector<long>(scaffold.begin(), scaffold.begin()+iditr);
                    scaffold1 = std::vector<long>(scaffold.begin()+iditr, scaffold.end());
                    if (mode == 0) {
                        return scaffold1;
                    } else if (mode == 1) {
                        std::reverse(scaffold1.begin(), scaffold1.end());
                        return scaffold1;
                    } else if (mode == 2) {
                        std::reverse(scaffold2.begin(), scaffold2.end());
                        return scaffold2;
                    } else if (mode == 3) {
                        return scaffold2;
                    }
                } else {
                    scaffold1 = std::vector<long>(scaffold.begin(), scaffold.begin()+iditr+1);
                    scaffold2 = std::vector<long>(scaffold.begin()+iditr+1, scaffold.end());
                    if (mode == 0) {
                        std::reverse(scaffold1.begin(), scaffold1.end());
                        return scaffold1;
                    } else if (mode == 1) {
                        return scaffold1;
                    } else if (mode == 2) {
                        return scaffold2;
                    } else if (mode == 3) {
                        std::reverse(scaffold2.begin(), scaffold2.end());
                        return scaffold2;
                    }
                }
            }
            return scaffold;
        }
        void uniteScaffolds(long id1, long id2)
        {
            long p1 = findParent(id1);
            long p2 = findParent(id2);
            if (p1 == p2) return;
            if (degrees[id1] != 1 || degrees[id2] != 1) {
                std::cerr << "failed to unite scaffolds" << std::endl;
                return;
            }
            if (sizes[p1] < sizes[p2])    std::swap(p1, p2);
            std::vector<long> scaffold1 = getScaffold(id1, 1);
            std::vector<long> scaffold2 = getScaffold(id2, 0);
            std::copy(scaffold2.begin(), scaffold2.end(), back_inserter(scaffold1));
            scaffolds[p1] = scaffold1;
            par[p2] = p1;
            sizes[p1] += sizes[p2];
            ++degrees[id1];
            ++degrees[id2];
        }
        bool judgeSameScaffold(long id1, long id2) { return (findParent(id1) == findParent(id2));}
    };
    // */

    //added by ouchi
    struct Alignment
    {
        long qStart;
        long qEnd;
        long tID;
        long tStart;
        long tEnd;
        long match;

        bool operator <(const Alignment &a) const {
            if (tID == a.tID)
                return (qStart == a.qStart) ? (qEnd < a.qEnd) : (qStart < a.qStart);
            else
                return (tID < a.tID);
        }

        Alignment(): qStart(0), qEnd(0), tID(0), tStart(0), tEnd(0), match(0) {}
        Alignment(long a, long b, long c, long d, long e, long f): qStart(a), qEnd(b), tID(c), tStart(d), tEnd(e), match(f) {}
        Alignment(const Alignment &) = default;
        Alignment &operator=(const Alignment &) = default;
        ~Alignment() = default;
    };
    //


    long contigMaxK;
	double heteroCoverage;
	double tolerenceFactor;
	double cutoffLengthFactor;
	std::vector<long> contigPreviousParentNodeID;
	std::vector<ContigBubbleInfo> contigBubbleInfo;
	unsigned mode;
	unsigned long numInputBubbleContig;
	std::vector<std::vector<long> > nodeBreakpointPosition;

    //added by ouchi
    long numConsensusNode;
    std::vector<ConsensusNode> consensusnode;
    std::vector<std::vector<std::pair<long, platanus::Position> > > nodePositionInConsensus;
    //long numHiCNode;
    //std::vector<HiCNode> hicnode;
    //std::vector<platanus::Position> consensusPositionInHiCNode;
    //

	// for node-state
    static const unsigned DBG_HETERO;
	static const unsigned DBG_PRIMARY_BUBBLE;
    static const unsigned DBG_SECONDARY_BUBBLE;

	// for edge-state
    static const char DBG_OVERLAP;

	// for contig-state
    static const char DBG_CONTIG_BUBBLE_JUNCTION;
	static const char DBG_CONTIG_PRIMARY_BUBBLE;
	static const char DBG_CONTIG_SECONDARY_BUBBLE;


	void storeGraphLinkFromOverlap(std::vector<GraphLinkWithFlag> &graphLinkPool);
	void storeGraphLinkFromMappedPair(std::vector<GraphLinkWithFlag> &graphLinkPool, long numThread);
	void storeGraphLinkFromMappedLongRead(std::vector<GraphLinkWithFlag> &graphLinkPool, long numThread);
    virtual void calcLink(const long libraryIndex, const long linkThreshold, const long numThread);
    virtual double calcNodeCoverage(const GraphNode &node);
	void calcLinkAndWriteGraphLinkWithFlagFile(const std::vector<GraphLinkWithFlag>& graphLinkPool, const GraphLinkWithFlagPoolIndex& index, const long libraryIndex, const bool multiThreadFlag);
	void makeHiCLink(const long numThread, const long HiCLinkMode = 0); //added by ouchi
	void getOverlappedBubbleNodeIndex(std::vector<long> &bubbleNodeIndex);
	void getOverlappedBubbleNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID);
	void getOverlappedForkNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID, std::vector<char> &directionBuffer);
	void getGappedBubbleNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID);
	void getGappedConflictingNodePairID(std::vector<std::pair<long, long> > &bubbleNodePairID, double heteroNodeCoverage);
	void markBubbleHeteroNode(const std::vector<long> &candidateNodeIndex, const double maxHeteroCoverageFactor);
	void calculateHeteroCoverage(const std::vector<long> &bubbleNodeIndex);
	long writeAndMarkOverlappedNodes(const std::vector<long> &nodeID, FILE *storeFP);
	long remakeGraphAccordingToPath(std::vector<GraphPath> &pathBufferForEachThread);
	long remakeGraphAccordingToPathPair(std::vector<GraphPath> &pathBufferForEachThread);
	void setOppositeBubbleNodeID(std::vector<std::array<long, 2> > &nodeIDVector, const std::vector<ScaffoldPart> &partVector);
	void flipOppositeBubbleNodeID(std::vector<std::array<long, 2> > &nodeIDVector, const std::vector<ScaffoldPart> &partVector);
	void setOppositeBubbleNodeIDStrandAware(std::vector<std::array<long, 2> > &nodeIDVector, const std::vector<ScaffoldPart> &partVector);
	long maxLengthContigID(std::vector<std::array<long, 2> > &IDs, const long start, const long end);
	std::pair<long, long> fillMajorityIDRunConvergenceAware(std::vector<std::array<long, 2> > &IDs, const GraphNode &targetNode, const std::pair<long, long> &ends, const double scoreFactor);
	long getNonGapContigLengthOfNode(const GraphNode &targetNode);
	long getNumEdgeDirectionOfNode(const GraphNode &targetNode);
	long deleteDifferentBubbleEdge(const long numThread);
	double calcNodeCoveragePartial(const GraphNode &node, const long start, const long end);
	void setCorrespondingNodePosition(std::vector<platanus::Position> &positionVector, const std::vector<ScaffoldPart> &partVector);
	void smoothNodeIDVector(std::vector<std::array<long, 2> > &nodeIDVector, const GraphNode &targetNode, const double scoreFactor);
	long deleteErroneousEdgeNumTagRate(const long numThread);
    void divideErroneousLink(const std::vector<std::vector<unsigned> >& numErroneousPair, const std::vector<std::vector<unsigned> >& numSpanningPair, const std::vector<std::vector<double> >& sumExpectedLink, std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink, const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize=0);
    void countPairsSpanningGap(std::vector<std::vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, const long numThread);
	void countLongReadSpanningGap(std::vector<std::vector<unsigned> >& numSpanningPair, const DIVISION_MODE divisionMode, long numThread);
    void countLinksInsideContigs(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void countLongReadLinksInsideContigs(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void countSwitchErrorLinks(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void countLongReadSwitchErrorLinks(std::vector<std::vector<unsigned> >& numPair, const long numThread);
	void calculatePhysicalCoverage(std::vector<std::vector<unsigned> >& physicalCoverage, const long insertTolerence, const long numThread);
	void calculateLongReadPhysicalCoverage(std::vector<std::vector<unsigned> >& physicalCoverage, const long insertTolerence, const long numThread);
	void compensatePhysicalCoverageBasedOnGapRate(std::vector<std::vector<unsigned> >& physicalCoverage, const long windowSize, const long numThread);
	void calculateDiffCoverage(std::vector<std::vector<unsigned> >& diffCoverage, const long lengthThreshold, const long numThread);
	void calculateLongReadDiffCoverage(std::vector<std::vector<unsigned> >& diffCoverage, const long lengthThreshold, const long numThread);
	void markJunctionContigJoinedToBubble();
	long getScoreFromIDPair(long leftNodeID, long rightNodeID);
	void detectBreakpointBasedOnCoverage(const std::vector<unsigned>& physicalCoverage, const std::vector<unsigned>& diffCoverage, const long edgeLength, const double minCoverageRate, const double maxDiffCoverageRate, const long minCoverage, std::vector<char>& breakpoint);
	void deleteShortRunOfBreakpoints(const long minRunLength, std::vector<char>& breakpoint);
	bool detectContigBoundaryBreakpoints(const long edgeLength, const GraphNode &targetNode, const std::vector<char>& baseBreakpoint, std::vector<char>& contigBreakpoint);
	void node2GapFlagsUnmappableContig(const GraphNode &node, std::vector<char> &ret);
	void storeBreakpointBasedOnCoverage(const std::vector<std::vector<unsigned> >& physicalCoverage, const std::vector<std::vector<unsigned> >& diffCoverage, const long numThread);
	bool checkNodeConvergenceBothSide(const long nodeID1, const long nodeID2);
	bool checkSingleNodeConvergenceBothSide(const long nodeID);


public:
	PairedDBG(): ScaffoldGraph(), contigMaxK(0), heteroCoverage(0.0), tolerenceFactor(1.0), cutoffLengthFactor(0.0), contigPreviousParentNodeID(0), contigBubbleInfo(0), mode(0x1 | 0x2), numInputBubbleContig(0), nodeBreakpointPosition(), numConsensusNode(0), consensusnode(), nodePositionInConsensus() {} //edited by ouchi (contigMak(0), numConsensusNode(0), consensusnode(), nodePositionInConsensus(), numHiCNode(0), hicnode(), consensusPositionInHiCNode())

    void setContigMaxK(const long len) { contigMaxK = len; }
	void setTolerenceFactor(const double factor) { tolerenceFactor = factor; }
	void setCutoffLengthFactor(const double factor) { cutoffLengthFactor = factor; }
	void setMode(const unsigned bits) { mode = bits; }
	void unsetMode(unsigned bits) { mode &= ~bits; }
	void setHeteroCoverage(const double cov) { heteroCoverage = cov; }
	void setNumInputBubbleContig(const unsigned long num) {numInputBubbleContig = num; }

	unsigned getMode() { return mode; }

    virtual void makeGraph(const long numThread);
	void makeGraphAllLibraries(const long numThread);
	void resetGraph();
	void getOverlappedNode(const long sourceNodeIndex, const char targetDirection, std::vector<long> &nodeIDBuffer);
	long getNumEdgeOneDirection(const GraphNode &targetNode, const char targetDirection);
	void extractDBGBubbleInformation();
	void calculateHeteroAndAverageCoverageUnphase();
	void markHeteroNode(const double maxHeteroCoverageFactor);
	void deleteNonOverlapHomoEdge();
	long joinUnambiguousNodePair(const long numThread);
	long joinUnambiguousNodePath(const long numThread);
	long solveSimpleCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	long solveUniquePathBetweenLinkedNodePair(const long numThread);
	void solveUniquePathBetweenLinkedNodePairAllLibrariesIterative(const long numThread);
	void searchUniqueOverlapPathGuidedByEdge(const long startNodeID, const GraphEdge &guideEdge, std::vector<GraphPath> &pathBuffer);
    virtual void cutAndPrintSeq(long minSeqLength, const unsigned long long readLength, const std::string &outFilename, const std::string &componentFilename);
	long crushSimpleDBGBubble();
    virtual unsigned long long crushHeteroBubble(const double averageCoverage);
    virtual long deleteHeteroEdge(void);
	virtual void loadResultSeq(const long minSeqLength, const unsigned long long readLength, const std::string prefix);
	void loadDividedContigResultSeq(const long minSeqLength, const unsigned long long readLength);
	void markRedundantResultSeq(const long numThread);
	void solveSimpleCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	void solveSimpleCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	void joinUnambiguousNodePairIterative(const long numThread);
	void joinUnambiguousNodePairGappedIterativeAllLibraries(const long numThread);
	long solveSimpleGappedCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	void solveSimpleGappedCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	void solveSimpleGappedCrossStructureAllLibrariesIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread);
	long solveDoubleGappedCrossStructure(const double linkRateThreshold, const CROSS_RESOLUTION_MODE resolutionMode, const long numThread); //added by ouchi
	void solveDoubleGappedCrossStructureIterative(const CROSS_RESOLUTION_MODE resolutionMode, const long numThread); //added by ouchi
	void solveUniquePathBetweenLinkedNodePairIterative(const long numThread);
	void etOppositeContigIDGappedConflicting(const long numThread);
	void setOppositeBubbleContigIDOverlapped(const long numThread);
	void setOppositeBubbleContigIDGapped(const long numThread);
	void setOppositeForkContigIDOverlapped(const long numThread);
	long divideNodeUsingBubbleContigPair(const long numThread);
	long divideNodeUsingBubbleContigPairStrandAware(const long numThread);
	void setOppositeBubbleNodeIDForEachNode(const long numThread);
	void setOppositeBubbleNodeIDAndStateForEachNode();
	void outputResultSeqWithBubble(const std::string filePrefix, const std::string &primaryBubbleSuffix, const std::string &secondaryBubbleSuffix, const std::string &primaryForkSuffix, const std::string &secondaryForkSuffix, const std::string &nestedBubbleSuffix, const std::string &nonBubbleOtherSuffix, const std::string &pairSuffix, const long minSeqLength, const long readLength);
	void outputResultSeqComponent(const std::string filePrefix, const std::string &fileSuffix);
	void deleteDifferentBubbleEdgeIterative(const long numThread);
	void deleteThinEdgeCostantKeepingOverlap(const long linkThreshold);
	long deleteConflictingBubbleEdge(const long numThread);
	long deleteSecondaryBubbleNodeAndEdge(const long numThread);
	long deleteShortAndLowCoverageBranch(const long lengthThreshold, const double coverageThreshold, const long numThread);
	void deleteShortAndLowCoverageBranchIterative(const long lengthThreshold, const double coverageThreshold, const long numThread);
	void setBubbleJunctionContigIDOverlapped();
	void setForkJunctionContigIDOverlapped();
	long divideBubbleJunctionNode(const bool gapDivideFlag);
	long divideBubbleContigInNonHeteroNode();
	void divideGappedNode(const long minGapSize);
	void copyAllNodes(PairedDBG &targetGraph);
	long getQuartileLengthOfBubble(const unsigned long quartileNumber);
    virtual void detectRepeat(const double averageCoverage);
    virtual void deleteRepeatEdge(void);
	void getUniqueConflictingNode(const long sourceNodeIndex, const char targetDirection, std::vector<NodeIDWithGap> &nodeIDBuffer);
	void divideNestedBubbleNode(const long numThread);
	void deleteLongEdge(const long maxEdgeLength);
	void deleteErroneousEdgeNumTagRateIterative(const long numThread);
	void deleteEdgeFromShortNodeKeepingBubble(const long lengthThreshold);
	void trimSparseEnd();
	void trimRepeatEnd(); //added by ouchi for Platanus-allee v2.3.2
	long divideInconsistentBubbleEnd();
	void adjustOppositeBubbleNodeIDDirection();
	void saveBubble();
	void deleteEdgeFromSecondaryBubble();
	void divideNodeBasedOnBubblesIterative(const bool strandFlag, const long numThread);
	void divideNodeByBubbleBoundary(const long numThread);
	void remakeGraphRecoveringSecondaryBubble(PairedDBG &bubbleGraph);
	void divideNodeByPrimaryBubbleBoundary(const PairedDBG &bubbleGraph);
	void divideErroneousNode(const long minLink, const DIVISION_MODE divisionMode, const long numThread, const long maxGapSize=0);
	void divideErroneousNodeBaseLevel(const long pairLibraryBegin, const long pairLibraryEnd, const long readLength, const bool bubbleFlag, const bool longLibraryFlag, const bool storeOnlyFlag, const long numThread);
    virtual void makeScaffold(void);
	void deleteEdgeFromDifferentPreviousParent(const long numThread);
	void clearContigPreviousParentNodeID();
	void makeScaffoldCombine(void); //added by ouchi for debug
	void setOppositeBubbleContigIDByEndMatch();
	void setOppositeBubbleContigIDByOneEndMatch();
    virtual void insertSizeDistribution(std::vector<SeqLib>& library, std::vector<long>& distribution, const long numThread);
	virtual void updateInsertLengthFP(std::vector<SeqLib>& lib, const long numThread);
	void dumpLongestNodeCoverage(std::vector<std::vector<unsigned> >& coverage, const long numOutputNode, std::string &outputFilename);
	void divideNodeBasedOnCoverage(const std::vector<std::vector<unsigned> >& physicalCoverage, const std::vector<std::vector<unsigned> >& diffCoverage, const bool bubbleFlag, const long numThread);

//added by ouchi
	void setOppositeBubbleContigIDByIndex();
	//long binsize;
	std::vector<double> idealcontact;
	//double percentilvalue;
	//void setBinsize(const long size) {binsize = size;}
	long deleteErroneousEdgebyHiC(const long numThread); //added by ouchi
	void deleteErroneousEdgebyHiCIterative(const long numThread); //added by ouchi
	bool checkHiCLinkBetweenNodePair(const long nodeID1, const char direction, const long nodeID2); //added by ouchi
	void calcIdealContact(const long numThread); //added by ouchi
	void calcPercentilValue(const long percent, const long numThread);
	void deleteGapRegionFromContactmap(std::vector<std::vector<int> > &contactmap, long &boundary, std::vector<long> &gapRegion);
	void makeContactmap(const long nodeID, std::vector<std::vector<int> > &contactmap); //added by ouchi
	long makeContactmap(const long nodeID1, const char direction, const long nodeID2, std::vector<std::vector<int> > &contactmap); //added by ouchi
	void calcSeparation(const std::vector<std::vector<int> > &contactmap, const long n1, const long n2, const long n3, std::array<double, 6> &result, bool divided); //added by ouchi
	void detectContactChange(const std::vector<std::vector<int> > &contactmap, const long start, const long end, std::vector<long> &locate); //added by ouchi
	template <typename T>
	void detectPeak(const std::vector<T> &v, const long delta, const long bin, const long maxh, std::vector<long> &locate); //added by ouchi
	void calculateHiCReadPhysicalCoverage(const long nodeID, std::vector<unsigned>& physicalCoverage);
	double getHiCLinkScoreBetweenNodePair(const long nodeID1, const char targetDirection, const long nodeID2, const long leftLimit, const long rightLimit);
	long getHiCLinkFromIDPair(long leftNodeID, long rightNodeID);
	void getConflictRegionBetweenNodePair(const GraphEdge &edge1, const GraphEdge &edge2, const GraphNode &node1, const GraphNode &node2, std::vector<char> &conflictpoint1, std::vector<char> &conflictpoint2);
	void getHeteroHiCLinkScore(long centerNodeIndex, const std::array<std::array<NodeIDWithGap ,2>, 2> &externalNodeID, std::array<double, 2> &sumLinkForHaplotype, std::array<double, 4> &HiCScores);
	void makeMPLink(const long numThread, const long pairLibraryBegin, const long pairLibraryEnd, const long LinkMode = 0);
	void makeMPINLink(const long numThread, const long pairLibraryBegin, const long pairLibraryEnd, const long LinkMode = 0);
	void makeLongLink(const long numThread);
	void makeLongINLink(const long numThread);
	void getHeteroMPLinkNum(long centerNodeIndex, const std::array<std::array<NodeIDWithGap ,2>, 2> &externalNodeID, std::array<double, 2> &sumLinkForHaplotype, std::array<double, 4> &MPScores);


	bool searchPathBetweenNodePair(const long nodeID1, const char targetDirection, const long nodeID2, const long limitLength);
	void detectRepeatSeverely(const double averageCoverage);
	bool checkNodeConflict(const GraphEdge &edge1, const GraphEdge &edge2, const GraphNode &node1, const GraphNode &node2);
	long joinUnambiguousUniqueNodePairGapped(const long numThread);
	void joinUnambiguousUniqueNodePairGappedIterative(const long numThread);
	long joinNoBranchNodePairGapped(const long numThread);
	void joinNoBranchNodePairGappedIterative(const long numThread);
	void joinNoBranchNodePairGappedIterativeAllLibraries(const long numThread);
	virtual long joinUnambiguousNodePairGapped(const long numThread); //added by ouchi

	void initConsensusNode(const long numThread);
	void remakeConsensus(const long numNewConsensusNode, const long newNodePoolsize, FILE *scaffoldFP);
	void writeSingletonConsensusNode(long &numNewConsensusNode, long &newNodePoolSize, FILE *storeFP);
	long writeAndMarkGappedConsensusNodes(const std::vector<ConsensusIDWithGap> &gappedPath, FILE *storeFP);
	long remakeConsensusGraphAccordingToGappedPath(std::vector<ConsensusPathGapped> &pathBuffer);
	bool judgeConflictConsensusNode(const long offset1, const long offset2, const long consensusID1, const long consensusID2, const char h1, const char h2, long tol);
	void setOppositeBubbleConsensusNode(const long consensusnodeIndex, std::array<long, 4>  &oppositeConsensusInfo);
	void mergeHeteroNode(const long numThread);
	void makeConsensusGraph(const long numThread, bool allLib = true);
	void getLinkedConsensusNode(const long sourceConsensusIndex, const char targetDirection, std::vector<ConsensusIDWithGap> &consensusIDBuffer);
	long joinUnambiguousConsensusPair(const long numThread);
	void joinUnambiguousConsensusPairIterative(const long numThread, bool allLib = true);
	void detectConflictConsensus(void);
	void consensusScaffolding(void);
	void deleteConsensusEdges(std::vector<long> &ids);
	void deleteErroneousConsensusEdgebyHiC(const long numThread);
	void divideErroneousConsensusNode(const long numThread);

	void divideNodeBasedOnHiC(long numThread);

	long makeConsensusContactmap(const long consensusID, std::vector<std::vector<int> > &contactmap);
	long makeConsensusContactmap(const long consensusID1, const char direction, const long consensusID2, std::vector<std::vector<int> > &contactmap);
	bool checkHiCLinkBetweenConsensusPair(const long consensusID1, const char direction, const long consensusID2);
	double calcConsensusCoverage(const ConsensusNode& consensusnode);
	void checkConsensusPath(long start, long end, long maxDistance, long &numLink, long &gap);
	void checkLinkBetweenHiCNodePair(std::vector<ConsensusPart> scaffold1, std::vector<ConsensusPart> scaffold2, std::vector<std::vector<long> > &physicalCoverage, std::vector<std::vector<long> > &diffCoverage, std::array<long, 5> &result);
	//void makeConsensusPathEdge(const long numThread);
	//void searchConsensusPath(long startID, std::vector<std::vector<long> > &path, std::vector<std::vector<char> > &h, std::vector<long> &lengths, std::vector<std::array<long, 2> > &numLinks);
	void HiC_Scaffolding(const long numThread, long lastFlag, std::string filePrefix);
	void calcHiCLinkScore(const long L, std::vector<HiCConsensusEdge> &es, long &numHiCNode, std::vector<HiCNode> &hicnode, std::vector<platanus::Position> &consensusPositionInHiCNode);
	void updateHiCLinkScore(std::vector<long> dividedIDs, std::vector<long> newIDs, long ei, const long L, std::vector<HiCConsensusEdge> &es, long &numHiCNode, std::vector<HiCNode> &hicnode, std::vector<platanus::Position> &consensusPositionInHiCNode);
	void compareHiCLinkScore(std::vector<ConsensusPart> scaffolda, std::vector<ConsensusPart> scaffoldb, std::vector<ConsensusPart> scaffolda2, std::vector<ConsensusPart> scaffoldb2);
	void getConsensusPart(long id, char h, char mode, std::vector<ConsensusPart>& scaffold, long &numHiCNode, std::vector<HiCNode> &hicnode, std::vector<platanus::Position> &consensusPositionInHiCNode);
	void uniteHiCNode(long id1, long id2, char h, long gap, long &numHiCNode, std::vector<HiCNode> &hicnode, std::vector<platanus::Position> &consensusPositionInHiCNode);
	std::array<long, 2> cutHiCNode(long id, char h, long offset, long &numHiCNode, std::vector<HiCNode> &hicnode, std::vector<platanus::Position> &consensusPositionInHiCNode, std::unordered_set<long>& dividedConsensusIDs, std::vector<std::vector<long> > &physicalCoverage, std::vector<std::vector<long> > &diffCoverage, long numThread);
	void calculateConsensusCoverage(std::vector<std::vector<long> > &physicalCoverage, std::vector<std::vector<long> > &diffCoverage, const long numThread);
	void detectBreakpointBasedOnConsensusCoverage(long consensusIndex, std::vector<long> physicalCoverage, std::vector<long> diffCoverage, char mode, std::vector<long> &breakpoint);
	long detectBreakpointBasedOnCoverage2(const std::vector<ConsensusPart> &scaffold, long start, long end);
    std::array<std::array<long, 2>, 2> divideConsensusNode(long consensusID, long offset, long numThread, char mode);
    std::array<std::array<long, 2>, 2> divideNode(long nodeID, long offset, long numThread, char mode);
    std::array<std::array<long, 2>, 2> divideContig(long contigID, long offset, long numThread);
	void phasing(const long numThread, const std::string filePrefix);
	void getHeteroMPLinkScoreBetweenScaffold(const std::array<std::array<std::vector<ScaffoldPart>, 2>, 2> &scaffolds, std::array<long, 2> &sumLinkForHaplotype);
	void getHeteroHiCLinkScoreBetweenScaffold(const std::array<std::array<std::vector<ScaffoldPart>, 2>, 2> &scaffolds, std::array<long, 2> &sumLinkForHaplotype);
	void getHeteroMPandHiCLinkScoreBetweenBlock(const std::vector<std::array<std::vector<ScaffoldPart>, 2> > &blocks, const std::vector<std::vector<long> > heteroBlocks, std::vector<std::vector<std::array<long, 2> > > &LinkScore);
	void consensus2seq(const ConsensusNode &consensusnode, std::vector<char> &ret, char mode = -1);
	bool checkHiCLinkBetweenHiCNodePair(std::vector<ConsensusPart> scaffold1, std::vector<ConsensusPart> scaffold2, std::array<std::vector<long>, 2> &breakpoint);
	long makeContactmapOfScaffold(std::vector<ConsensusPart> scaffold1, std::vector<ConsensusPart> scaffold2, std::vector<std::vector<int> >& contactmap);
	void outputConsensusFastg(const std::string &outputFilename);
	void outputHiCNode(const std::string &outFilename, const long numHiCNode, const std::vector<HiCNode> &hicnode);
	void outputConsensusCovTxt(const std::string &outputFilename, const long numThread);
	void checkHiCNode(const long numHiCNode, const std::vector<HiCNode> &hicnode);

	void mincingBubble(const std::string filePrefix, const long numThread);
	void calculateBaseCoverage(std::vector<std::vector<double> >& baseCoverage, const std::string PAFFilename, const long numThread);
	void mincingBubbleNodeBySelfAlignment(const std::string PAFFilename, const long numThread);
	bool setOppositeBubbleConsensusNodebySelfAlignment(const long consensusnodeIndex, const std::vector<std::vector<Alignment> > &selfAlignments, std::array<long, 4> &oppositeConsensusInfo);
//


	double getHeteroCoverage() { return heteroCoverage; }


	static const unsigned OVERLAP_MODE;
    static const unsigned PAIRED_END_LINK_MODE;
    static const unsigned LENGTH_CUTOFF_MODE;
    static const unsigned NON_DBG_OVERLAP_MODE;
	static const unsigned LONG_READ_LINK_MODE;
	static const unsigned BUBBLE_AWARE_MODE;
	static const unsigned SECONDARY_BUBBLE_REMOVAL_MODE;
	static const unsigned PREVIOUS_DIVISION_AWARE_MODE;

	static const long MAX_ITERATION_OF_CROSS_SOLUTION;

    static const double NON_PAIRED_END_TOLERENCE_FACTOR;
    static const double HETERO_COVERAGE_THRESHOLD_FACTOR;
    static const double HETERO_FORK_COVERAGE_THRESHOLD_FACTOR;
	static const double CROSS_LINK_RATE_THRESHOLD;
	static const double CROSS_SCORE_RATE_THRESHOLD;
	static const double MIN_BUBBLE_COUNT_FACTOR;

	//added by ouchi
	static const double HiC_CONTACTMAP_S_THRESHOLD;
	static const double HiC_CONTACTMAP_T_THRESHOLD;
	static const long HiC_CONTACTMAP_L_THRESHOLD;
	static const long HiC_CONTACTMAP_BINSIZE;
	//
};


#endif
