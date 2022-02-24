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

#ifndef GRAPH_H
#define GRAPH_H

#include "common.h"
#include "counter.h"
#include "kmer.h"
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <stack>
#include <set>
#include <bitset>
#include <utility>


struct NodeBase
{
    long id;
    NodeBase(): id(0) {}
    virtual ~NodeBase() {}
};




//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Straight class
// this class has contiguity kmer which doesn't have junction in de Bruijn graph
// out variable means terminal information (usually terminal has junction or no connection)
// length doesn't have actual straight seq length, (actual straight seq length) + (kmer length) - 1.
struct Straight : public NodeBase
{
    u64_t length;
    unsigned short coverage;
    unsigned char out;
    std::vector<unsigned long long> seq;

    Straight(): NodeBase(), length(0), coverage(0), out(0), seq() {}
    ~Straight() {}

    bool operator==(const Straight &stra) const { return seq == stra.seq; }


    // base setter and getter in seq
    unsigned char get(unsigned pos) const { return (seq[pos/32] >> ((pos%32)*2)) & 3; }
    void set(unsigned pos, unsigned char val)
    {
        seq[pos/32] = (seq[pos/32] & ~(0x3ull << ((pos%32)*2))) | (static_cast<u64_t>(val) << ((pos%32)*2));
    }

    // resize seq std::vector
    void resize(unsigned kmerLength)
    {
        seq.resize((length + kmerLength + 30) / 32);
    }

    // read and write seq in temporary file
    void readTemporaryFileSeq(FILE *fp)
    {
        for (unsigned long long i = 0, n = seq.size(); i < n; ++i)
            fread(&(seq[i]), sizeof(unsigned long long), 1, fp);
    }

    void writeTemporaryFileSeq(FILE *fp) const
    {
        for (unsigned long long i = 0, n = seq.size(); i < n; ++i)
            fwrite(&(seq[i]), sizeof(unsigned long long), 1, fp);
    }
};
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Junction class
// this class has junction in de Bruijn graph
// out variable means connected kmer information
struct Junction : public NodeBase
{
    unsigned short coverage;
    unsigned char out;

    Junction(): NodeBase(), coverage(0), out(0) {}
    ~Junction() {}
};
//////////////////////////////////////////////////////////////////////////////////////





template <typename KMER>
class BruijnGraph
{
    typedef unsigned long long u64_t;
public:
    enum DIRECTION {LTOR=0, RTOL, DUAL};

struct DuplicateJunction
{
    typename KMER::keyType leftStraightKey;
    typename KMER::keyType junctionKey;
    typename KMER::keyType rightStraightKey;
    DIRECTION direction;

    DuplicateJunction(const typename KMER::keyType &lPtr,  const typename KMER::keyType &jun, const typename KMER::keyType &rPtr, const DIRECTION dir): leftStraightKey(lPtr), junctionKey(jun), rightStraightKey(rPtr), direction(dir) {}
    bool operator==(const DuplicateJunction &dup) const { return (leftStraightKey == dup.leftStraightKey) && (junctionKey == dup.junctionKey) && (rightStraightKey == dup.rightStraightKey);}

  struct hasher
  {
    inline size_t operator()(DuplicateJunction j) const { return (size_t)(j.leftStraightKey & UINT64_MAX); }
  };
};

private:
    static const double BUBBLE_COVERAGE_RATE;
    unsigned OVERLAP_READ_LEN;
    unsigned kmerLength;
    u64_t numJunction;
    u64_t numStraight;
    std::unordered_map<typename KMER::keyType, Junction, typename KMER::hasher> junctionTable;
    std::unordered_map<typename KMER::keyType, Straight, typename KMER::hasher> leftStraightTable;
    std::unordered_map<typename KMER::keyType, Straight *, typename KMER::hasher> rightStraightTable;
    double bubbleThreshold;
    double branchThreshold;
    FILE *junctionFP;
    FILE *straightFP;
    FILE *bubbleFP;



    void searchBubbleStructure(Junction *rootJunction, std::vector<Junction*> &junction, std::vector<Straight*> &straight, const typename KMER::keyType &buffer);
    unsigned long long pairwiseAlignment(const Straight * const straight1, const Straight * const straight2) const;
    inline void deleteStraight(const Straight &straight);
    inline void deleteStraight(const Straight &straight, KMER &kmer);
    inline void deleteJunction(const typename KMER::keyType &key);
    inline void insertStraight(const Straight &sp);
    inline void insertStraight(const Straight &sp, KMER &kmer);
    inline Straight concatinateS_S(Straight *leftStraight, Straight *rightStraight);
    inline Straight concatinateS_J(const Straight *straight, const typename KMER::keyType &junctionKey);
    inline Straight concatinateJ_J(const typename KMER::keyType &leftJuncKey, const typename KMER::keyType &rightJuncKey);
    inline Straight changeJunctionToStraight(const typename KMER::keyType &key);
    unsigned long long getLeftMinimal(const std::vector<unsigned long long> &dist) const;

    void makeKmerFromNode(DoubleHash<typename KMER::keyType, platanus::Position> &kmerAssigner);

    void mapRead(const DoubleHash<typename KMER::keyType, platanus::Position> &kmerAssigner, std::vector<std::set<std::pair<long, long> > > &straightFillRead, FILE **readFP, const unsigned long long numThread);
    void divideNode(std::vector<std::set<std::pair<long, long> > > &straightFillRead);
    inline bool correctlyMapRead(DoubleHash<typename KMER::keyType, platanus::Position> &kmerAssigner, std::vector<platanus::Position> &nodeSetForward, std::vector<platanus::Position> &nodeSetReverse, platanus::SEQ &seq, KMER &kmer, const unsigned long long mask) const;
    inline std::unordered_map<long, std::pair<long, long> > gatherMappingPosition(std::vector<platanus::Position> &nodeSet);
    inline void updateFillRead(std::unordered_map<long, std::pair<long, long> >&result, std::vector<std::set<std::pair<long, long> > > &straightFillRead);
    inline void remakeStraightNode(const std::pair<long, long> &point, std::vector<Straight> &straightVector, const Straight &straight);

public:
    BruijnGraph(): OVERLAP_READ_LEN(0), kmerLength(0), numJunction(0), numStraight(0), junctionTable(), leftStraightTable(), rightStraightTable(), bubbleThreshold(), branchThreshold(), junctionFP(NULL), straightFP(NULL) {}
    ~BruijnGraph() {};
    BruijnGraph(const BruijnGraph &b) = delete;
    BruijnGraph *operator=(const BruijnGraph &b) = delete;

    // setter and finder

    void setBubbleAndBranch(const double bubble, const double branch)
    {
        bubbleThreshold = bubble;
        branchThreshold = branch;
    }

    void setJunctionTable(const Junction &junc, const typename KMER::keyType &key) {junctionTable[key] = junc; }

    void setKmerLength(const unsigned k) {kmerLength = k; }

    Straight *findStraightHead(const typename KMER::keyType &key)
    {
        auto it = leftStraightTable.find(key);
        if (it == leftStraightTable.end() || it->second.coverage == UINT16_MAX)
            return NULL;
        else
            return &(it->second);
    }

    Straight *findStraightTail(const typename KMER::keyType &key)
    {
        auto it = rightStraightTable.find(key);
        if (it == rightStraightTable.end() || it->second->coverage == UINT16_MAX)
            return NULL;
        else
            return it->second;
    }

    Junction *findJunction(const typename KMER::keyType &key)
    {
        auto it = junctionTable.find(key);
        if (it == junctionTable.end() || it->second.coverage == UINT16_MAX)
            return NULL;
        else
            return &(it->second);
    }

    // clear function
    void deleteAllTable(void)
    {
        junctionTable.clear();
        leftStraightTable.clear();
        rightStraightTable.clear();
    }

    std::unordered_map<typename KMER::keyType, Straight, typename KMER::hasher> moveStraightTable(void)
    {
        rightStraightTable.clear();
        return std::move(leftStraightTable);
    }

    typename std::unordered_map<typename KMER::keyType, Straight, typename KMER::hasher>::const_iterator getLeftStraightTableBegin(void)
    {
        return leftStraightTable.begin();
    }

    typename std::unordered_map<typename KMER::keyType, Straight, typename KMER::hasher>::const_iterator getLeftStraightTableEnd(void)
    {
        return leftStraightTable.end();
    }

    unsigned long long estimateNumKmerOnStraight(void) const
    {
        unsigned long long num = 0;

        for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
            if (it->second.coverage == UINT16_MAX) continue;
            num += it->second.length;
        }
        return num;
    }



    void makeInitialBruijnGraph(Counter<KMER> &counter, FILE *sortedKeyFP);
    unsigned long long crushBubble(const double averageCoverage);
    unsigned long long cutBranch(const unsigned long long numThread=1);
    void concatinateNodes(void);
    void printBubble(const std::string &outputFilename, const double occurrenceRatio) const;
    void printContig(const std::string &outputFilename, const double coverageRatio, const double averageLength) const;
    template <typename NEXT> void saveContig(const unsigned nextKmer, const double occurrenceRatio, DoubleHash<typename NEXT::keyType, unsigned short> &table);
    void saveContigSimple(FILE *contigFP, const double coverageRatio) const;
	void saveJunction(FILE *contigFP, const double coverageRatio) const;
    void cutBranchIterative(const unsigned long long numThread=1);
    void crushBubbleIterative(const double averageCoverage);
    void saveEdgeKmer(Counter<KMER> &counter, const unsigned nextKmerLength) const;
    void calcStraightStats(unsigned long long &lengthCut, unsigned long long &coverageCut, double &averageCoverage) const;
    unsigned long long getLeftMinimalCoverage() const;
    unsigned long long deleteErroneousStraightNode(const unsigned long long lengthCut, const unsigned long long coverageCut, const unsigned long long numThread);
    void deleteErroneousStraightNodeIterative(const unsigned long long lengthCut, const unsigned long long coverageCut, const unsigned long long numThread);
    void makeBruijnGraphForGapClose(Counter<KMER> &counter, const long minCoverage);
    template <typename NEXT> void saveLargeKmerForGapClose(Counter<NEXT> &counter);
    double getAverageCoverageExcludingBubble(void);


    void divideStraightNode(FILE **readFP, const unsigned long long coverageCutoff, const unsigned long long size, const unsigned long long numThread);
};

template<typename KMER>
const double BruijnGraph<KMER>::BUBBLE_COVERAGE_RATE = 1.5;


//////////////////////////////////////////////////////////////////////////////////////
// make initial brujin graph using kmers
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::makeInitialBruijnGraph(Counter<KMER> &counter, FILE *sortedKeyFP)
{
    unsigned char rightFlags, leftFlags;
    unsigned long long sum = 0;
    Junction leftJunction;
    Straight leftStraight;
    std::vector<char> leftSeq, rightSeq;
    std::unordered_map<unsigned long long, typename KMER::keyType> kmerStack;
    unsigned long long stackID = 0;
    KMER kmer(kmerLength);
    KMER startKmer(kmerLength);
    KMER leftKmer(kmerLength);
    KMER rightKmer(kmerLength);
    typename KMER::keyType key(kmerLength);
    typename KMER::keyType nextKey(kmerLength), tmpKey(kmerLength);
    KMER reservedKmer(kmerLength);
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));


    junctionTable.clear();
    leftStraightTable.clear();
    rightStraightTable.clear();

    numStraight = numJunction = 0;
    std::cerr << "connecting kmers..." << std::endl;

    leftSeq.reserve(100000);
    rightSeq.reserve(100000);

    if (junctionFP != NULL)
        fclose(junctionFP);
    junctionFP = platanus::makeTemporaryFile();


    if (straightFP != NULL)
        fclose(straightFP);
    straightFP = platanus::makeTemporaryFile();


    rewind(sortedKeyFP);
    while (startKmer.readTemporaryFileForward(sortedKeyFP)) {
        stackID = 0;
        kmerStack[stackID] = (startKmer.forward);
        ++stackID;
        while (stackID > 0) {
            startKmer.forward = kmerStack[stackID - 1];
            startKmer.reverseComplement();
            unsigned short kmerValue = counter.findValueAndOverwrite(std::min(startKmer.forward, startKmer.reverse), UINT16_MAX);
            if (kmerValue == UINT16_MAX) {
                --stackID;
                continue;
            }
            sum = kmerValue;

            // search parent kmer node 
            // example kmer "ATTGCAGC", find "*ATTGCAG"
            leftSeq.clear();
            leftKmer = startKmer;
            leftKmer.forward >>= 2;
            leftKmer.reverse <<= 2;
            leftKmer.maskReverse(mask);

            // search whether parent kmer node contains valueTable.
            leftFlags = 0;
            for (unsigned char base = 0; base < 4; ++base) {
                leftKmer.setForward(leftKmer.kmerLength - 1, base);
                leftKmer.setReverse(0, 0x3^base);
                if (counter.findValue(std::min(leftKmer.forward, leftKmer.reverse)) != 0) {
                    leftFlags |= 1 << base;
                }
            }
            if (leftFlags != 0)
                leftSeq.push_back(platanus::Flag2Base(leftFlags));
           // search parent kmer node
           // example kmer "ATTGCAGC", find "TTGCAGC*"
            rightKmer = startKmer;
            rightKmer.forward <<= 2;
            rightKmer.maskForward(mask);
            rightKmer.reverse >>= 2;
            rightSeq.clear();
            // search whether parent kmer node contains valueTable.
            rightFlags = 0;
            for (unsigned char base = 0; base < 4; ++base) {
                rightKmer.setForward(0, base);
                rightKmer.setReverse(rightKmer.kmerLength - 1, 0x3^base);
                if (counter.findValue(std::min(rightKmer.forward, rightKmer.reverse)) != 0) {
                    rightFlags |= 1 << base;
                }
            }
            if (rightFlags != 0)
                rightSeq.push_back(platanus::Flag2Base(rightFlags));

            // if left/fight seq value is 4, stacked and and write junction file.
            if ((leftSeq.size() > 0 && leftSeq[0] == 4) || (rightSeq.size() > 0 && rightSeq[0] == 4)) {
                --stackID;
                // kmerStack.pop_back();
                for (unsigned char base = 0; base < 4; ++base) {
                    if ((rightFlags & (1 << base)) != 0) {
                        rightKmer.setForward(0, base);
                        kmerStack[stackID] = (rightKmer.forward);
                        ++stackID;
                    }
                    if ((leftFlags & (1 << base)) != 0) {
                        leftKmer.setForward(leftKmer.kmerLength - 1, base);
                        kmerStack[stackID] = (leftKmer.forward);
                        ++stackID;
                    }
                }

                leftJunction.coverage = sum;
                leftJunction.out = (leftFlags << 4) | rightFlags;
                setJunctionTable(leftJunction, startKmer.forward);
                ++numJunction;
                continue;
            }

            // extend de brujin graph, which has single edge, in right side.
            if (rightSeq.size() > 0) {
                // set kmer key
                kmer = startKmer;
                kmer.forward <<= 2;
                kmer.maskForward(mask);
                kmer.reverse >>= 2;
                kmer.setForward(0, rightSeq[0]);
                kmer.setReverse(kmer.kmerLength - 1, 0x3^rightSeq[0]);
                key = std::min(kmer.forward, kmer.reverse);
                auto keyIt = counter.getOccurrenceIterator(key);
                auto nextIt = keyIt;
                while (1) {
                    // check whether only one node is found in reverse relation.
                    // i.e. we check one son "ATGGCAG" in before phase
                    // ATGGCAG -> AATGGCA, check AATGGCA <- only ATGGCAG, not exist ATGGCAA, ATGGCAC, ATGGCAT
                    rightKmer = kmer;
                    kmer.forward >>= 2;
                    kmer.reverse <<= 2;
                    kmer.maskReverse(mask);
                    unsigned char tmpFlags = 0;
                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(kmer.kmerLength - 1, base);
                        kmer.setReverse(0, 0x3^base);
                        if (counter.findValue(std::min(kmer.forward, kmer.reverse)) != 0) {
                            tmpFlags |= 1 << base;
                        }
                    }
                    if (platanus::FlagCount(tmpFlags) > 1 || keyIt->second == UINT16_MAX) {
                        rightFlags = 0 | (1 << rightSeq.back());
                        rightSeq.pop_back();
                        break;
                    }

                    // search parent kmer node "iterative"
                    // example kmer "ATTGCAGC", find "*ATTGCAG"
                    kmer.forward <<= 4;
                    kmer.maskForward(mask);
                    kmer.reverse >>= 4;
                    kmer.setForward(1, rightSeq.back());
                    kmer.setReverse(kmer.kmerLength - 2, 0x3^rightSeq.back());
                    rightFlags = 0;

                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(0, base);
                        kmer.setReverse(kmer.kmerLength - 1, 0x3^base);
                        tmpKey = std::min(kmer.forward, kmer.reverse);
                        auto tmpIt = counter.getOccurrenceIterator(tmpKey);
                        if (tmpIt != counter.getOccurrenceEnd() && tmpIt->second != 0) {
                            rightFlags |= 1 << base;
                            nextKey = tmpKey;
                            nextIt = tmpIt;
                        }
                    }

                    if (rightFlags == 0) {
                        sum += keyIt->second;
                        keyIt->second = static_cast<unsigned short>(UINT16_MAX);
                        break;
                    } else if (platanus::FlagCount(rightFlags) > 1) {
                        rightFlags = 0 | (1 << rightSeq.back());
                        rightSeq.pop_back();
                        break;
                    }

                    rightSeq.push_back(platanus::Flag2Base(rightFlags));
                    sum += keyIt->second;
                    keyIt->second = static_cast<unsigned short>(UINT16_MAX);
                    key = nextKey;
                    keyIt = nextIt;
                }
            }

            // same process in if (rightSeq.size() > 0)
            if (leftSeq.size() > 0) {
                kmer = startKmer;
                kmer.forward >>= 2;
                kmer.reverse <<= 2;
                kmer.maskReverse(mask);
                kmer.setForward(kmer.kmerLength - 1, leftSeq[0]);
                kmer.setReverse(0, 0x3^leftSeq[0]);
                key = std::min(kmer.forward, kmer.reverse);
                auto keyIt = counter.getOccurrenceIterator(key);
                auto nextIt = keyIt;
                while (1) {
                    leftKmer = kmer;
                    kmer.forward <<= 2;
                    kmer.maskForward(mask);
                    kmer.reverse >>= 2;
                    unsigned char tmpFlags = 0;
                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(0, base);
                        kmer.setReverse(kmer.kmerLength - 1, 0x3^base);
                        if (counter.findValue(std::min(kmer.forward, kmer.reverse)) != 0) {
                            tmpFlags |= 1 << base;
                        }
                    }
                    if (platanus::FlagCount(tmpFlags) > 1 || keyIt->second == UINT16_MAX) {
                        leftFlags = 0 | (1 << leftSeq.back());
                        leftSeq.pop_back();
                        break;
                    }
                    kmer.forward >>= 4;
                    kmer.reverse <<= 4;
                    kmer.maskReverse(mask);
                    kmer.setForward(kmer.kmerLength - 2, leftSeq.back());
                    kmer.setReverse(1, 0x3^leftSeq.back());
                    leftFlags = 0;

                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(kmer.kmerLength - 1, base);
                        kmer.setReverse(0, 0x3^base);
                        tmpKey = std::min(kmer.forward, kmer.reverse);
                        auto tmpIt = counter.getOccurrenceIterator(tmpKey);
                        if (tmpIt != counter.getOccurrenceEnd() && tmpIt->second != 0) {
                            leftFlags |= 1 << base;
                            nextKey = tmpKey;
                            nextIt = tmpIt;
                        }
                    }
                    if (leftFlags == 0) {
                        sum += keyIt->second;
                        keyIt->second = static_cast<unsigned short>(UINT16_MAX);
                        break;
                    } else if (platanus::FlagCount(leftFlags) > 1) {
                        leftFlags = 0 | (1 << leftSeq.back());
                        leftSeq.pop_back();
                        break;
                    }

                    leftSeq.push_back(platanus::Flag2Base(leftFlags));
                    sum += keyIt->second;
                    keyIt->second = static_cast<unsigned short>(UINT16_MAX);
                    key = nextKey;
                    keyIt = nextIt;
                }
            }
            --stackID;

            if (leftFlags != 0) {
                leftKmer.setForward(leftKmer.kmerLength - 1, platanus::Flag2Base(leftFlags));
                kmerStack[stackID] = (leftKmer.forward);
                ++stackID;
            }

            if (rightFlags != 0) {
                rightKmer.setForward(0, platanus::Flag2Base(rightFlags));
                kmerStack[stackID] = (rightKmer.forward);
                ++stackID;
            }

            // make straight struct
            leftStraight.length = leftSeq.size() + rightSeq.size() + 1;
            binstr_t seq(leftStraight.length + kmerLength);
            for (unsigned j = 0, n = leftSeq.size(); j < n; ++j) {
                seq.set(leftSeq.size() - 1 - j, leftSeq[leftSeq.size() - 1 - j]);
            }
            seq <<= 2 * kmerLength;
            for (unsigned j = 0; j < kmerLength; ++j) {
                seq.set(kmerLength - 1 - j, (startKmer.forward >> (2 * (kmerLength - 1 - j))) & 0x3);
            }
            seq <<= 2 * rightSeq.size();
            for (unsigned j = 0, n = rightSeq.size(); j < n; ++j) {
                seq.set(rightSeq.size() - 1 - j, rightSeq[j]);
            }
            leftStraight.coverage = (static_cast<double>(sum) / static_cast<double>(leftStraight.length) + 0.5);
            leftStraight.out = (leftFlags << 4) | rightFlags;
            seq.convertToVector(leftStraight.seq);
            this->insertStraight(leftStraight, reservedKmer);
            ++numStraight;
        }
    }

}


//////////////////////////////////////////////////////////////////////////////////////
// crush cubble
// bubble means like below structure
//      |---------------straight-------|
// rootJunciton                      junction
//      |---------------straight-------|
// bubbleClusterID means what straight regars bubble.
// This variable value contains minimum base number in bubble cluster (the set of straight which regard bubble each other)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long BruijnGraph<KMER>::crushBubble(const double averageCoverage)
{
    typedef unsigned long long u64_t;
    u64_t numCrush = 0;
    u64_t alignmentThreshold = 0;
    u64_t bubbleClusterID[4];
    std::vector<Junction*> junction(4);
    std::vector<Straight*> straight(4);
    std::vector<Straight*> maxStraight(4);
    std::vector<KMER> key(4, KMER(kmerLength));
    Junction *rootJunction;
    KMER reservedKmer(kmerLength);

    unsigned short coverageThreshold = static_cast<unsigned short>(std::min(static_cast<u64_t>(averageCoverage * BUBBLE_COVERAGE_RATE + 0.5), static_cast<u64_t>(UINT16_MAX)));
	if (averageCoverage >= UINT16_MAX)
		coverageThreshold = UINT16_MAX;

    for (auto it = junctionTable.begin(), end = junctionTable.end(); it != end; ++it) {
        rootJunction = &(it->second);
        if (rootJunction->coverage == UINT16_MAX || platanus::FlagCount(rootJunction->out & 0x0f) < 2) continue;

        for (unsigned char base = 0; base < 4; ++base) {
            bubbleClusterID[base] = base;
            straight[base] = NULL;
            junction[base] = NULL;
        }

        searchBubbleStructure(rootJunction, junction, straight, it->first);

        for (unsigned char base1 = 0; base1 < 3; ++base1) {
            if (straight[base1] == NULL) continue;
            for (unsigned char base2 = base1 + 1; base2 < 4; ++base2) {
                if (straight[base2] == NULL || junction[base1] == NULL || junction[base1] != junction[base2] || straight[base1]->coverage + straight[base2]->coverage > coverageThreshold) continue;
                alignmentThreshold = static_cast<u64_t>((std::max(straight[base1]->length, straight[base2]->length) + kmerLength - 1) * bubbleThreshold + 0.5);

                // if straight length is too short, regard bubble
                if (std::min(straight[base1]->length, straight[base2]->length) + 1 <= kmerLength) {
                    if (abs(static_cast<long>(straight[base1]->length) - static_cast<long>(straight[base2]->length)) > alignmentThreshold) continue;
                    if (bubbleClusterID[base2] == base2)
                        bubbleClusterID[base2] = bubbleClusterID[base1];
                    else
                        bubbleClusterID[base1] = bubbleClusterID[base2];
                    continue;
                }

                // pairwise alignment
                if (pairwiseAlignment(straight[base1], straight[base2]) > alignmentThreshold) continue;
                if (bubbleClusterID[base2] == base2)
                    bubbleClusterID[base2] = bubbleClusterID[base1];
                else
                    bubbleClusterID[base1] = bubbleClusterID[base2];
            }
        }

        // decide maximum coverage in cluster without itself for each straight
        for (unsigned char base1 = 0; base1 < 3; ++base1) {
            maxStraight[base1] = NULL;
            for (unsigned char base2 = base1; base2 < 4; ++base2) {
                if (bubbleClusterID[base2] != base1 || straight[base2] == NULL) continue;
                if (maxStraight[base1] == NULL || straight[base2]->coverage > maxStraight[base1]->coverage)
                    maxStraight[base1] = straight[base2];
            }
        }

        for (unsigned char base1 = 0; base1 < 3; ++base1) {
            if (maxStraight[base1] == NULL) continue;
            for (unsigned char base2 = base1; base2 < 4; ++base2) {
                if (bubbleClusterID[base2] != base1 || straight[base2] == NULL || straight[base2] == maxStraight[base1]) continue;

                u64_t tmp = maxStraight[base1]->coverage * maxStraight[base1]->length + straight[base2]->coverage * straight[base2]->length;
                tmp = static_cast<double>(tmp) / maxStraight[base1]->length + 0.5;
                maxStraight[base1]->coverage = std::min(tmp, static_cast<u64_t>(UINT16_MAX - 1));
                rootJunction->out ^= 1 << base2;
                junction[base2]->out ^= 1 << (straight[base2]->get(kmerLength - 1) + 4);

                if (bubbleFP != NULL) {
                    tmp = straight[base2]->length + kmerLength - 1;
                    fwrite(&tmp, sizeof(u64_t), 1, bubbleFP);
                    fwrite(&(straight[base2]->coverage), sizeof(unsigned short), 1, bubbleFP);
                    binstr_t seq(tmp);
                    seq.convertFromVector(straight[base2]->seq);
                    seq.writeTemporaryFile(bubbleFP);
                }
                deleteStraight(*straight[base2], reservedKmer);
                ++numCrush;
            }
        }
    }
    return numCrush;
}


//////////////////////////////////////////////////////////////////////////////////////
// search bubble structure crush cubble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::searchBubbleStructure(Junction *rootJunction, std::vector<Junction*> &junction, std::vector<Straight*> &straight, const typename KMER::keyType &buffer)
{
    KMER key(kmerLength);
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    key.forward = buffer;
    key.forward <<= 2;
    key.maskForward(mask);
    for (unsigned char base = 0; base < 4; ++base) {
        if (!(rootJunction->out & (1 << base))) continue;
        key.setForward(0, base);
        straight[base] = findStraightHead(key.forward);
        if (straight[base] != NULL && (straight[base]->out & 0x0f)) {
            for (unsigned j = 0; j < kmerLength - 1; ++j)
                key.setReverse(j + 1, straight[base]->get(j));
            key.setReverse(0, platanus::Flag2Base(straight[base]->out & 0x0f));
            junction[base] = findJunction(key.reverse);
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// pairwise alignment between bubbles
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long BruijnGraph<KMER>::pairwiseAlignment(const Straight * const straight1, const Straight * const straight2) const
{
    unsigned long long alignmentSize = kmerLength;
    std::vector<u64_t> alignmentScore(alignmentSize * 2);

    unsigned long long alignmentLength1 = straight1->length - kmerLength + 1;
    unsigned long long alignmentLength2 = straight2->length - kmerLength + 1;
    if (alignmentLength2 + 1 > alignmentSize) {
        alignmentSize = alignmentLength2 + 1;
        alignmentScore.resize(alignmentSize * 2);
    }


    u64_t m, n;
    for (n = 0; n < alignmentLength2 + 1; ++n)
        alignmentScore[n] = n;
    for (m = 0; m < alignmentLength1; ++m) {
        alignmentScore[(m & 1) * alignmentSize] = m;
        alignmentScore[(~m & 1) * alignmentSize] = m + 1;
        for (n = 0; n < alignmentLength2; ++n) {
            if (straight1->get(straight1->length - 1 - m) == straight2->get(straight2->length - 1 - n)) {
                alignmentScore[(~m & 1) * alignmentSize + n + 1] = alignmentScore[(m & 1) * alignmentSize + n];
                continue;
            }
            alignmentScore[(~m&1)*alignmentSize + n+1] = alignmentScore[(m&1)*alignmentSize + n] + 1;
            if (alignmentScore[(~m&1)*alignmentSize + n+1] > alignmentScore[(m&1)*alignmentSize + n+1] + 1)
                alignmentScore[(~m&1)*alignmentSize + n+1] = alignmentScore[(m&1)*alignmentSize + n+1] + 1;
            if (alignmentScore[(~m&1)*alignmentSize + n+1] > alignmentScore[(~m&1)*alignmentSize + n] + 1)
                alignmentScore[(~m&1)*alignmentSize + n+1] = alignmentScore[(~m&1)*alignmentSize + n] + 1;
        }
    }

    return alignmentScore[(m&1)*alignmentSize + n];
}



//////////////////////////////////////////////////////////////////////////////////////
// delete straight in table
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::deleteStraight(const Straight &straight)
{
    KMER kmer(kmerLength);
    deleteStraight(straight, kmer);
}



//////////////////////////////////////////////////////////////////////////////////////
// delete straight in table
// this function add KMER argument compared with above function
// it spends to make KMER constructor too many time,
// so this function diminishs KMER constructor to be given KMER argument
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::deleteStraight(const Straight &straight, KMER &kmer)
{
    for (unsigned i = 0; i < kmerLength; ++i) {
        // kmer.setReverse(i, straight.get(i));
        kmer.setForward(i, straight.get(straight.length + i - 1));
    }

    leftStraightTable[kmer.forward].coverage = UINT16_MAX;

}



//////////////////////////////////////////////////////////////////////////////////////
// delete junction in table
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::deleteJunction(const typename KMER::keyType &key)
{
    junctionTable[key].coverage = UINT16_MAX;
}


//////////////////////////////////////////////////////////////////////////////////////
// insert straight in table
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::insertStraight(const Straight &sp)
{
    KMER kmer(kmerLength);
    insertStraight(sp, kmer);
}


//////////////////////////////////////////////////////////////////////////////////////
// insert straight in table
// this function add KMER argument compared with above function
// it spends to make KMER constructor too many time,
// so this function diminishs KMER constructor to be given KMER argument
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::insertStraight(const Straight &sp, KMER &kmer)
{
    for (unsigned i = 0; i < kmerLength; ++i) {
        kmer.setReverse(i, sp.get(i));
        kmer.setForward(i, sp.get(sp.length + i - 1));
    }


    leftStraightTable[kmer.forward] = sp;
    rightStraightTable[kmer.reverse] = &(leftStraightTable[kmer.forward]);

}



//////////////////////////////////////////////////////////////////////////////////////
// cut graph branch
// branch means the short straight node which has little coverage
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long BruijnGraph<KMER>::cutBranch(const unsigned long long numThread)
{
    unsigned long long numDelete = 0;
    Straight *straight;
    KMER key(kmerLength);
    omp_lock_t lock;
    omp_init_lock(&lock);
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));

    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1) private(straight) firstprivate(key) reduction(+: numDelete)
    for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
        unsigned long long i = 0;
        for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it, ++i) {
            if (i % numThread != threadID) continue;
            straight = &(it->second);
            if (straight->coverage == UINT16_MAX || straight->length > kmerLength) continue;


            // try to cut branch left side (expand forward)
            // the highest coverage is rested, others are deleted.
            if (straight->out >> 4) {
                if (straight->out & 0xf) continue;

                key.forward = it->first;
                key.forward >>= 2;
                key.setForward(kmerLength - 1, platanus::Flag2Base(straight->out >> 4));
                Junction *junction = findJunction(key.forward);
                if (junction == NULL) continue;

                key.forward <<= 2;
                key.maskForward(mask);
                unsigned long long maxCoverage = 0;

                //search branch node and decide the highest coverage node
                for (unsigned char base = 0; base < 4; ++base) {
                    if (!(junction->out & (1 << base))) continue;

                    key.setForward(0, base);
                    Straight *tmpStraight = findStraightHead(key.forward);
                    if (tmpStraight != straight && tmpStraight != NULL && tmpStraight->coverage > maxCoverage) {
                        maxCoverage = tmpStraight->coverage;
                        continue;
                    }

                    Junction *tmpJunction = findJunction(key.forward);
                    if (tmpJunction != NULL && tmpJunction->coverage > maxCoverage)
                        maxCoverage = tmpJunction->coverage;
                }
                if (straight->coverage > maxCoverage * branchThreshold) continue;

                // overwritte junction node information (cut branch actually)
                omp_set_lock(&lock);
                junction->out ^= 1 << straight->get(straight->length - 1);
                omp_unset_lock(&lock);

            } else if (straight->out & 0xf) {
                for (unsigned j = 1; j < kmerLength; ++j)
                    key.setForward(j, straight->get(j - 1));
                key.setForward(0, platanus::Flag2Base(straight->out & 0xf));
                Junction *junction = findJunction(key.forward);
                if (junction == NULL) continue;

                key.forward >>= 2;
                unsigned long long maxCoverage = 0;

                //search branch node and decide the highest coverage node
                for (unsigned char base = 0; base < 4; ++base) {
                    if(!(junction->out & (1 << (base + 4)))) continue;

                    key.setForward(kmerLength - 1, base);
                    Straight *tmpStraight = findStraightTail(key.forward);
                    if (tmpStraight != straight && tmpStraight != NULL && tmpStraight->coverage > maxCoverage) {
                        maxCoverage = tmpStraight->coverage;
                        continue;
                    }

                    Junction *tmpJunction = findJunction(key.forward);
                    if (tmpJunction != NULL && tmpJunction->coverage > maxCoverage)
                        maxCoverage = tmpJunction->coverage;
                }
                if (straight->coverage > maxCoverage * branchThreshold) continue;

                // overwritte junction node information (cut branch actually)
                omp_set_lock(&lock);
                junction->out ^= 1 << (4 + straight->get(kmerLength - 1));
                omp_unset_lock(&lock);
            }
            else continue;

            deleteStraight(*straight);
            ++numDelete;
        }
    }
    omp_destroy_lock(&lock);
    return numDelete;
}


//////////////////////////////////////////////////////////////////////////////////////
// when delete node (straight or junction)
// retry connecting nodes
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::concatinateNodes(void)
{
    Straight newStraight;
    KMER key(kmerLength);

    // retry concatinate nodes in junction side
    for (auto it = junctionTable.begin(), end = junctionTable.end(); it != end; ++it) {
        Junction *junction = &(it->second);
        // check this junction connects one node at most for ecah side
        if (junction->coverage == UINT16_MAX || platanus::FlagCount(junction->out >> 4) > 1 || platanus::FlagCount(junction->out & 0xf) > 1) continue;
        Straight *leftStraight, *rightStraight;
        Junction *leftJunction, *rightJunction;

        if (platanus::FlagCount(junction->out >> 4)) {
            key.forward = it->first;
            key.forward >>= 2;
            key.setForward(kmerLength - 1, platanus::Flag2Base(junction->out >> 4));

            // check connected node ecists straight node
            leftStraight = findStraightTail(key.forward);

            // find straight for the first time
            if (leftStraight == NULL) {
                leftJunction = findJunction(key.forward);

                if (leftJunction == NULL || platanus::FlagCount(leftJunction->out >> 4) > 1 || platanus::FlagCount(leftJunction->out & 0xf) > 1) {
                    leftJunction = NULL;
                }

            }
        } else {
            leftStraight = NULL;
            leftJunction = NULL;
        }

        if (leftStraight != NULL) {
            newStraight = concatinateS_J(leftStraight, it->first);
            leftStraight = &newStraight;
        } else if (leftJunction != NULL) {
            newStraight = concatinateJ_J(key.forward, it->first);
            leftStraight = &newStraight;
        } else {
            newStraight = changeJunctionToStraight(it->first);
            leftStraight = &newStraight;
        }

        if (platanus::FlagCount(leftStraight->out & 0xf)) {
            for (unsigned j = 0; j < kmerLength - 1; ++j)
                key.setForward(j + 1, leftStraight->get(j));
            key.setForward(0, platanus::Flag2Base(leftStraight->out & 0xf));
            rightStraight = findStraightHead(key.forward);
            if (rightStraight == NULL) {
                rightJunction = findJunction(key.forward);

                if (rightJunction == NULL || platanus::FlagCount(rightJunction->out >> 4) > 1 || platanus::FlagCount(rightJunction->out & 0xf) > 1)
                    rightJunction = NULL;
            }
        } else {
            rightStraight = NULL;
            rightJunction = NULL;
        }

        if (rightStraight != NULL){
            concatinateS_S(leftStraight, rightStraight);
            }
        else if (rightJunction != NULL){
            concatinateS_J(leftStraight, key.forward);
        }
    }


    KMER key2(kmerLength);
    for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
        Straight *straight = &(it->second);
        Straight newStraight;

        if (straight->coverage == UINT16_MAX) continue;
        if (platanus::FlagCount(straight->out >> 4)) {
            key2.forward = it->first;
            key2.forward >>= 2;
            key2.setForward(kmerLength - 1, platanus::Flag2Base(straight->out >> 4));
            Straight *leftStraight = findStraightTail(key2.forward);
            if (leftStraight != NULL){
                newStraight = concatinateS_S(leftStraight, straight);
            }
                straight = &newStraight;

        }

        if (platanus::FlagCount(straight->out & 0xf)) {
            for (unsigned j = 1; j < kmerLength; ++j)
                key2.setForward(j, straight->get(j - 1));
            key2.setForward(0, platanus::Flag2Base(straight->out & 0xf));
            Straight *rightStraight = findStraightHead(key2.forward);
            if (rightStraight != NULL){
                concatinateS_S(straight, rightStraight);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// connect straight node and straight one actually
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
Straight BruijnGraph<KMER>::concatinateS_S(Straight *leftStraight, Straight *rightStraight)
{
    KMER reservedKmer(kmerLength);
    if (*leftStraight == *rightStraight)
        return *leftStraight;

    // make new node
    Straight newStraight;
    newStraight.length = leftStraight->length + rightStraight->length;
    newStraight.coverage = static_cast<double>(leftStraight->coverage * leftStraight->length + rightStraight->coverage * rightStraight->length)
                        / static_cast<double>(newStraight.length) + 0.5;
    newStraight.out = (leftStraight->out & 0xf0) | (rightStraight->out & 0x0f);
    newStraight.resize(kmerLength);
    for (unsigned i = 0; i < leftStraight->length + kmerLength - 1; ++i)
        newStraight.set(newStraight.length + kmerLength - 2 - i, leftStraight->get(leftStraight->length + kmerLength - 2 - i));
    for (unsigned i = 0; i < rightStraight->length; ++i)
        newStraight.set(rightStraight->length - 1 - i, rightStraight->get(rightStraight->length - 1 - i));
    // delete old node and insert new one
    deleteStraight(*leftStraight, reservedKmer);
    deleteStraight(*rightStraight, reservedKmer);
    insertStraight(newStraight, reservedKmer);
    return std::move(newStraight);
}

//////////////////////////////////////////////////////////////////////////////////////
// connect straight node and junction one actually
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
Straight BruijnGraph<KMER>::concatinateS_J(const Straight *straight, const typename KMER::keyType &junctionKey)
{
    KMER reservedKmer(kmerLength);

    // make new node
    Straight newStraight;
    Junction junction = junctionTable[junctionKey];
    newStraight.length = straight->length + 1;
    newStraight.coverage = static_cast<double>(straight->coverage * straight->length + junction.coverage) / static_cast<double>(newStraight.length) + 0.5;
    newStraight.out = (straight->out & 0xf0) | (junction.out & 0x0f);
    newStraight.resize(kmerLength);
    newStraight.set(0, static_cast<unsigned char>(junctionKey & 0x3));
    for (unsigned long long i = 0; i < straight->length + kmerLength - 1; ++i)
        newStraight.set(i + 1, straight->get(i));
    deleteStraight(*straight, reservedKmer);
    newStraight.id = straight->id;

    // delete old node and insert new one
    deleteJunction(junctionKey);
    insertStraight(newStraight, reservedKmer);
    return std::move(newStraight);
}

//////////////////////////////////////////////////////////////////////////////////////
// connect junction node and junction one actually
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
Straight BruijnGraph<KMER>::concatinateJ_J(const typename KMER::keyType &leftJuncKey, const typename KMER::keyType &rightJuncKey)
{
    if (leftJuncKey == rightJuncKey)
        return changeJunctionToStraight(leftJuncKey);

    // make new node
    Junction leftJunction = junctionTable[leftJuncKey];
    Junction rightJunction = junctionTable[rightJuncKey];
    Straight newStraight;
    newStraight.length = 2;
    newStraight.coverage = static_cast<double>(leftJunction.coverage + rightJunction.coverage) / 2.0 + 0.5;
    newStraight.out = (leftJunction.out & 0xf0) | (rightJunction.out & 0x0f);
    newStraight.resize(kmerLength);
    newStraight.set(kmerLength, static_cast<unsigned char>((leftJuncKey >> (2 * (kmerLength - 1))) & 0x3));
    for (unsigned i = 0; i < kmerLength; ++i) {
        newStraight.set(kmerLength - 1 - i, static_cast<unsigned char>((rightJuncKey >> (2 * (kmerLength - 1 - i))) & 0x3));
    }

    // delete old node and insert new one
    deleteJunction(leftJuncKey);
    deleteJunction(rightJuncKey);
    insertStraight(newStraight);
    return std::move(newStraight);
}


//////////////////////////////////////////////////////////////////////////////////////
// change junction node to straight one
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
Straight BruijnGraph<KMER>::changeJunctionToStraight(const typename KMER::keyType &key)
{
    // make new node
    Junction junction = junctionTable[key];
    Straight newStraight;
    newStraight.length = 1;
    newStraight.coverage = junction.coverage;
    newStraight.out = junction.out;
    newStraight.resize(kmerLength);

    for (unsigned i = 0; i < kmerLength; ++i)
        newStraight.set(i, static_cast<unsigned char>((key >> (2 * i)) & 0x3));

    // delete old node and insert new one
    deleteJunction(key);
    insertStraight(newStraight);
    return std::move(newStraight);
}


//////////////////////////////////////////////////////////////////////////////////////
// divide straight node
// all reads is mapped and merge mapping results consider kmer overlap
// if all seq not covered, divide this seq
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::divideStraightNode(FILE **readFP, const unsigned long long coverageCutoff, const unsigned long long size, const unsigned long long numThread)
{
    DoubleHash<typename KMER::keyType, platanus::Position> kmerAssigner(size);
    std::vector<std::set<std::pair<long, long> > > straightFillRead(leftStraightTable.size() + 1);
    this->makeKmerFromNode(kmerAssigner);
    mapRead(kmerAssigner, straightFillRead, readFP, numThread);
    DoubleHash<typename KMER::keyType, platanus::Position>().swap(kmerAssigner);
    divideNode(straightFillRead);
}



//////////////////////////////////////////////////////////////////////////////////////
// make kmer assinger table form de Bruijn graph nodes
// divide each node seq to kmer
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::makeKmerFromNode(DoubleHash<typename KMER::keyType, platanus::Position> &kmerAssigner)
{
    // make junction kmer table
    // junctionID must be negative value
    long junctionID = 0;
    platanus::Position pos;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    KMER kmer(this->kmerLength);
    for (auto it = junctionTable.begin(), end = junctionTable.end(); it != end; ++it) {
        if (it->second.coverage == UINT16_MAX) continue;
        --junctionID;
        kmer.forward = it->first;
        kmer.reverseComplement();
        pos.offset = 0;
        pos.id = junctionID;
        kmerAssigner[kmer.forward] = pos;
        kmerAssigner[kmer.reverse] = pos;
        it->second.id = junctionID;
    }

    // make straight kmer table
    // straightID must be positive value
    long straightID = 0;
    binstr_t str;
    for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
        if (it->second.coverage == UINT16_MAX) continue;
        str.resize(it->second.length + kmerLength - 1);
        str.convertFromVector(it->second.seq);
        ++straightID;
        it->second.id = straightID;
        kmer.forward = it->first;
        pos.offset = 0;
        pos.id = straightID;
        kmerAssigner[kmer.forward] = pos;
        for (unsigned long long i = 0; i < it->second.length - 1; ++i) {
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.setForward(0, str.get(it->second.length - 2 - i));
            pos.offset = i + 1;
            kmerAssigner[kmer.forward] = pos;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// map all reads and divide straight node
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::mapRead(const DoubleHash<typename KMER::keyType, platanus::Position> &kmerAssigner, std::vector<std::set<std::pair<long, long> > > &straightFillRead, FILE **readFP, const unsigned long long numThread)
{
    omp_lock_t updateLock;
    omp_init_lock(&updateLock);
    unsigned long long numMapped = 0;
    unsigned long long numUnmapped = 0;
    unsigned long long numShort = 0;
        const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));

    OVERLAP_READ_LEN = kmerLength;

    std::cerr << "mapping reads on de Bruijn Graph nodes..." << std::endl;

    omp_set_num_threads(numThread);

    // map all reads
    # pragma omp parallel for schedule(static, 1) reduction(+: numMapped, numUnmapped, numShort)
    for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
        platanus::SEQ seq;
        KMER kmer(kmerLength);
        rewind(readFP[threadID]);
        while(seq.readTemporaryFile(readFP[threadID])) {
            if (seq.length < kmerLength) {
                ++numShort;
                continue;
            }
            bool mapForward = true;
            bool mapReverse = true;
            std::vector<platanus::Position> nodeSetForward(seq.length - this->kmerLength + 1);
            std::vector<platanus::Position> nodeSetReverse(seq.length - this->kmerLength + 1);
            seq.positionUnknown.resize(seq.numUnknown + 1);
            seq.positionUnknown[seq.numUnknown] = platanus::ConstParam::MAX_READ_LEN + 1;
            seq.numUnknown = 0;
            for (unsigned j = 0; j < kmerLength - 1; ++j) {
                kmer.setForward(kmerLength - 2 - j, seq.base[j]);
                kmer.setReverse(j + 1, 0x3 ^ seq.base[j]);
            }
            for (unsigned j = 0; j < seq.length - kmerLength + 1; ++j) {
                kmer.forward <<= 2;
                kmer.maskForward(mask);
                kmer.setForward(0, seq.base[j + kmerLength - 1]);
                kmer.reverse >>= 2;
                kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[j + kmerLength - 1]);
                if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) < j + kmerLength) {
                    if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) <= j)
                        ++seq.numUnknown;
                    continue;
                }
                if (mapForward && mapReverse) {
                    auto itForward = kmerAssigner.find_any(kmer.forward);
                    auto itReverse = kmerAssigner.find_any(kmer.reverse);
                    if (itForward->second && itReverse->second) {
                        if (itForward->second.id > 0 || itReverse->second.id > 0) {
                            mapForward = false;
                            mapReverse = false;
                            break;
                        }
                    } else if (itForward->second) {
                        nodeSetForward[j] = itForward->second;
                        mapReverse = false;
                    } else if (itReverse->second) {
                        nodeSetReverse[j] = itReverse->second;
                        mapForward = false;
                    } else {
                        mapForward = false;
                        mapReverse = false;
                        break;
                    }
                } else if (mapForward) {
                    auto itForward = kmerAssigner.find_any(kmer.forward);
                    if (itForward->second) {
                        if (itForward->second.id < 0) {
                            mapReverse = true;
                        }
                        nodeSetForward[j] = itForward->second;
                    } else {
                        mapForward = false;
                        break;
                    }
                } else if (mapReverse) {
                    auto itForward = kmerAssigner.find_any(kmer.reverse);
                    if (itForward->second) {
                        if (itForward->second.id < 0) {
                            mapForward = true;
                        }
                        nodeSetReverse[j] = itForward->second;
                    } else {
                        mapReverse = false;
                        break;
                    }
                } else {
                    break;
                }
            }

            if (mapForward || mapReverse) {
                std::unordered_map<long, std::pair<long, long> > mappingResult = gatherMappingPosition(nodeSetForward);
                omp_set_lock(&updateLock);
                this->updateFillRead(mappingResult, straightFillRead);
                omp_unset_lock(&updateLock);

                mappingResult.clear();
                std::reverse(nodeSetReverse.begin(), nodeSetReverse.end());
                mappingResult = gatherMappingPosition(nodeSetReverse);
                omp_set_lock(&updateLock);
                this->updateFillRead(mappingResult, straightFillRead);
                omp_unset_lock(&updateLock);
                ++numMapped;
            } else {
                ++numUnmapped;
            }
        }
    }
    omp_destroy_lock(&updateLock);

    std::cerr << "TOTAL_MAPPED_READS=" << numMapped << std::endl;
    std::cerr << "TOTAL_UNMAPPED_READS=" << numUnmapped << std::endl;
    std::cerr << "TOTAL_SHORT_READS(<" << this->kmerLength << ")=" << numShort << std::endl;
}

template <typename KMER>
void BruijnGraph<KMER>::divideNode(std::vector<std::set<std::pair<long, long> > > &straightFillRead)
{
    // divide straight node
    std::vector<Straight> straightVector;
    unsigned long long numDelete = 0;
    unsigned long long numCut = 0;
    KMER kmer(kmerLength);
    // make new straight node (not add actually)
    for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
        if (it->second.coverage == UINT16_MAX) continue;
        if (straightFillRead[it->second.id].size() > 1) {
            Straight stra(it->second);
            deleteStraight(it->second, kmer);
            ++numCut;
            for (auto i = straightFillRead[it->second.id].begin(), n = straightFillRead[it->second.id].end(); i != n; ++i) {
                remakeStraightNode(*i, straightVector, stra);
            }
        } else if (straightFillRead[it->second.id].size() == 0) {
            ++numDelete;
            it->second.coverage = UINT16_MAX;
        } else {
            if (straightFillRead[it->second.id].begin()->first != 0 || it->second.length + kmerLength - 2 - straightFillRead[it->second.id].begin()->second != 0) {
                ++numCut;
                Straight stra(it->second);
                deleteStraight(it->second, kmer);
                remakeStraightNode(*(straightFillRead[it->second.id].begin()), straightVector, stra);
            }
        }
    }
    // add straight node actually
    for (auto it = straightVector.begin(), end = straightVector.end(); it != end; ++it) {
        insertStraight(*it);
    }

    std::cerr << "NUM_DELETE_NODE(reads are unmapped)=" << numDelete << std::endl;
    std::cerr << "NUM_CUT_NODE=" << numCut << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// check all kmer made from read mapped correctly
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline bool BruijnGraph<KMER>::correctlyMapRead(DoubleHash<typename KMER::keyType, platanus::Position> &kmerAssigner, std::vector<platanus::Position> &nodeSetForward, std::vector<platanus::Position> &nodeSetReverse, platanus::SEQ &seq, KMER &kmer, const unsigned long long mask) const
{
    bool mapForward = true;
    bool mapReverse = true;
    for (unsigned j = 0; j < kmerLength - 1; ++j) {
        kmer.setForward(kmerLength - 2 - j, seq.base[j]);
        kmer.setReverse(j + 1, 0x3 ^ seq.base[j]);
    }
    auto nodeSetFIt = nodeSetForward.begin();
    auto nodeSetRIt = nodeSetReverse.begin();
    for (unsigned j = 0; j < seq.length - this->kmerLength + 1; ++j, ++nodeSetFIt, ++nodeSetRIt) {
        kmer.forward <<= 2;
        kmer.maskForward(mask);
        kmer.setForward(0, seq.base[j + kmerLength - 1]);
        kmer.reverse >>= 2;
        kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[j + kmerLength - 1]);
        if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) < j + kmerLength) {
            if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) <= j)
                ++seq.numUnknown;
            continue;
        }
        // when first kmer or previous kmer mapped on junction
        if (mapForward && mapReverse) {
            auto itForward = kmerAssigner.find(kmer.forward);
            auto itReverse = kmerAssigner.find(kmer.reverse);
            auto end = kmerAssigner.end();
            if (itForward != end && itReverse != end) {
                if (itForward->second.id > 0 || itReverse->second.id > 0) {
                    return false;
                }
            } else if (itForward != end) {
                *nodeSetFIt = itForward->second;
                mapReverse = false;
            } else if (itReverse != end) {
                *nodeSetRIt = itReverse->second;
                mapForward = false;
            } else {
                return false;
            }
        // when previous kmer mapped on straight node forward direction
        } else if (mapForward) {
            auto itForward = kmerAssigner.find(kmer.forward);
            if (itForward != kmerAssigner.end()) {
                if (itForward->second.id < 0) {
                    mapReverse = true;
                }
                *nodeSetFIt = itForward->second;
            } else {
                return false;
            }
        // when previous kmer mapped on straight node reverse direction
        } else if (mapReverse) {
            auto itForward = kmerAssigner.find(kmer.reverse);
            if (itForward != kmerAssigner.end()) {
                if (itForward->second.id < 0) {
                mapForward = true;
                }
                *nodeSetRIt = itForward->second;
            } else {
                return false;
            }
        } else {
          return false;
        }
    }
    return true;
}


//////////////////////////////////////////////////////////////////////////////////////
// remake straight node and push this on straightVector
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::remakeStraightNode(const std::pair<long, long> &point, std::vector<Straight> &straightVector, const Straight &straight)
{
    if (point.first < 0 || point.second > static_cast<long>(straight.length) + kmerLength - 2) return;
    if (point.second - point.first + 2 - kmerLength > 0) {
        Straight newStraight(straight);
        binstr_t seqOld(straight.length + kmerLength - 1);
        seqOld.convertFromVector(newStraight.seq);
        binstr_t seqNew(point.second - point.first + 1);
        for (long i = point.first; i <= point.second; ++i) {
            seqNew.set(point.second - i, seqOld.get(straight.length + kmerLength - 2 - i));
        }
        newStraight.length = seqNew.len - kmerLength + 1;
        seqNew.convertToVector(newStraight.seq);
        newStraight.out = 0;
        if (point.first == 0) newStraight.out |= (straight.out & 0xf0);
        if (point.second == static_cast<long>(straight.length) + kmerLength - 2) newStraight.out |= (straight.out & 0x0f);
        straightVector.push_back(newStraight);
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// gather kmer mapping result and make read mapping result
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline std::unordered_map<long, std::pair<long, long> > BruijnGraph<KMER>::gatherMappingPosition(std::vector<platanus::Position> &nodeSet)
{
    std::unordered_map<long, std::pair<long, long> > result;
    auto it = nodeSet.begin();
    long nowID = it->id;
    long start = it->offset;
    long end = it->offset + kmerLength - 1;
    for (auto endIt = nodeSet.end(); it != endIt; ++it) {
        if (it->id == 0) continue;
        if (it->id == nowID) {
            if (it->offset < start) start = it->offset;
            if (it->offset + kmerLength - 1 > end) end = it->offset + kmerLength - 1;
        } else {
            if (nowID > 0) {
                result[nowID] = std::make_pair(start, end);
            }
            nowID = it->id;
            start = it->offset;
            end = it->offset + kmerLength - 1;
        }
    }
    if (nowID > 0) {
        result[nowID] = std::make_pair(start, end);
    }
    return result;
}





//////////////////////////////////////////////////////////////////////////////////////
// update fill read
// fill read show each straight node's position covered reads consider kmer overlap
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void BruijnGraph<KMER>::updateFillRead(std::unordered_map<long, std::pair<long, long> >&result, std::vector<std::set<std::pair<long, long> > > &straightFillRead)
{
    for (auto it = result.begin(), endIt = result.end(); it != endIt; ++it) {
        straightFillRead[it->first].insert(it->second);
        bool concatinateFront = false;
        bool concatinateBack = false;
        bool remove = false;
        bool removeBack = false;
        auto insertIt = straightFillRead[it->first].find(it->second);
        auto frontIt = insertIt;
        auto backIt = insertIt;
        auto end = straightFillRead[it->first].end();
        ++backIt;
        while (backIt != end) {
            if (backIt->second <= insertIt->second) {
                straightFillRead[it->first].erase(backIt);
                backIt = insertIt;
                ++backIt;
            } else {
                break;
            }
        }
        if (insertIt != straightFillRead[it->first].begin()) {
            --frontIt;
            if (frontIt->second >= insertIt->second) {
                remove = true;
            } else if (frontIt->second >= insertIt->first + OVERLAP_READ_LEN - 1) {
                concatinateFront = true;
            }
        }
        backIt = insertIt;
        ++backIt;
        if (backIt != straightFillRead[it->first].end()) {
            if (backIt->second <= insertIt->second) {
                removeBack = true;
            } else if (backIt->first + OVERLAP_READ_LEN - 1 <= insertIt->second) {
                concatinateBack = true;
            }
        }

        if (remove) {
            straightFillRead[it->first].erase(insertIt);
        } else if (removeBack) {
            straightFillRead[it->first].erase(backIt);
        } else if (concatinateFront && concatinateBack) {
            long start = frontIt->first;
            long finish = backIt->second;
            straightFillRead[it->first].erase(frontIt);
            straightFillRead[it->first].erase(insertIt);
            straightFillRead[it->first].erase(backIt);
            straightFillRead[it->first].insert(std::make_pair(start, finish));
        } else if (concatinateFront) {
            long start = frontIt->first;
            long finish = insertIt->second;
            straightFillRead[it->first].erase(frontIt);
            straightFillRead[it->first].erase(insertIt);
            straightFillRead[it->first].insert(std::make_pair(start, finish));
        } else if (concatinateBack) {
            long start = insertIt->first;
            long finish = backIt->second;
            straightFillRead[it->first].erase(backIt);
            straightFillRead[it->first].erase(insertIt);
            straightFillRead[it->first].insert(std::make_pair(start, finish));
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// print bubble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::printBubble(const std::string &outputFilename, const double occurrenceRatio) const
{
    unsigned long long length;
    unsigned long long numSeq = 0;
    std::ofstream out(outputFilename.c_str());

    rewind(bubbleFP);
    while(fread(&length, sizeof(unsigned long long), 1, bubbleFP)) {
        unsigned short occurrence = 0;
        unsigned long long i;
        binstr_t seq(length);
        fread(&occurrence, sizeof(unsigned short), 1, bubbleFP);
        seq.readTemporaryFile(bubbleFP);
        occurrence = occurrence * occurrenceRatio + 0.5;
        ++numSeq;
        out << ">seq" << numSeq << "_len" << length << "_cov" << occurrence << "\n";
        for (i = 0; i < length; ++i) {
            out.put(platanus::Bin2Char(seq.get(length - 1 - i)));
            if ((i + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
        }
        if (i % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }
    out.close();
}


//////////////////////////////////////////////////////////////////////////////////////
// print contig
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::printContig(const std::string &outputFilename, const double coverageRatio, const double averageLength) const
{
    std::ofstream out(outputFilename.c_str());
    std::vector<typename KMER::keyType> buffer;
    unsigned long long numSeq = 0;

    for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
        unsigned long long j;
        const Straight &straight = it->second;
        if (straight.coverage == UINT16_MAX) continue;
        ++numSeq;
        out << ">seq" << numSeq << "_len" << straight.length + kmerLength - 1 << "_cov" << static_cast<unsigned short>(straight.coverage * coverageRatio + 0.5) << "_read" << static_cast<long>(averageLength + 0.5) << "\n";
        for (j = 0; j < straight.length + kmerLength - 1; ++j) {
            out.put(platanus::Bin2Char(straight.get(straight.length + kmerLength - 2 - j)));
            if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                out.put('\n');
        }
        if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }
    out.close();
}


//////////////////////////////////////////////////////////////////////////////////////
// save contig
// when extend from one straight terminal, if it can be connected the straight, add extenstion to straight terminal
//
// image
//
// 1 ----------------- ---|
//                       junc  ------------------
// 2 ----------------- ---|
//
// it seems to see one straight from straight 1 seq (seq1 cannot see seq2 because of same direction connected in junction)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
template <typename NEXT>
void BruijnGraph<KMER>::saveContig(const unsigned nextKmer, const double occurrenceRatio, DoubleHash<typename NEXT::keyType, unsigned short> &table)
{
    unsigned long long diff = nextKmer - kmerLength;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    const unsigned long long nextMask = nextKmer >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * nextKmer));
    auto tableEnd = leftStraightTable.end();
    KMER key(kmerLength);
    NEXT next(nextKmer);
    typename NEXT::keyType selectedKey(nextKmer);
    for (auto tableIt = leftStraightTable.begin(); tableIt != tableEnd; ++tableIt) {
        const Straight *straight = &(tableIt->second);
        Straight writeStraight;
        if (straight->coverage == UINT16_MAX) {
            continue;
        }
        writeStraight.length = 0;
        writeStraight.coverage = straight->coverage;
        binstr_t str;


        // check left terminal
        if (straight->out >> 4) {
            key.forward = tableIt->first;
            key.forward >>= 2;
            key.setForward(kmerLength - 1, platanus::Flag2Base(straight->out >> 4));

            // check connective junction exists
            Junction *junction = findJunction(key.forward);
            if (junction != NULL) {
                if (platanus::FlagCount(junction->out >> 4) == 1) {
                    key.reverse = key.forward;
                    key.reverse >>= 2;
                    key.setReverse(kmerLength - 1, platanus::Flag2Base(junction->out >> 4));

                    // check connective straight exists
                    Straight *tmpStraight = findStraightTail(key.reverse);
                    if (tmpStraight != NULL) {
                        unsigned exLength = std::min(tmpStraight->length, diff);
                        str.resize(writeStraight.length + exLength);
                        // write tmpStraight seq terminal
                        str <<= 2 * exLength;
                        for (unsigned i = 0; i < exLength; ++i) {
                            str.set(exLength - 1- i, tmpStraight->get(kmerLength + exLength - 2 - i)); // tmpStraight terminal is same seq in straight
                        }
                        writeStraight.length += exLength;
                    }
                }

                // write junction seq
                str.resize(writeStraight.length + 1);
                str <<= 2;
                str.set(0, (key.forward >> (2 * (kmerLength - 1))) & 0x3);
                ++writeStraight.length;

            }
        }

        // write straight seq
        str.resize(writeStraight.length + straight->length + kmerLength - 1);
        str <<= 2 * (straight->length + kmerLength - 1);
        for (unsigned i = 0; i < straight->length + kmerLength - 1; ++i) {
            str.set(straight->length + kmerLength - 2 - i, straight->get(straight->length + kmerLength - 2 - i));
        }
        writeStraight.length += straight->length + kmerLength - 1;

        // check right terminal
        if (straight->out & 0xf) {
            for (unsigned i = 1; i < kmerLength; ++i)
                key.setForward(i, straight->get(i - 1));
            key.setForward(0, platanus::Flag2Base(straight->out & 0xf));
            Junction *junction = findJunction(key.forward);
            if (junction != NULL) {
                str.resize(writeStraight.length + 1);
                str <<= 2;
                str.set(0, key.forward & 0x3);
                ++writeStraight.length;
                if (platanus::FlagCount(junction->out & 0xf) == 1) {
                    key.forward <<= 2;
                    key.maskForward(mask);
                    key.setForward(0, platanus::Flag2Base(junction->out & 0xf));
                    Straight *tmpStraight = findStraightHead(key.forward);
                    if (tmpStraight != NULL) {
                        unsigned exLength = std::min(tmpStraight->length, diff);
                        str.resize(writeStraight.length + exLength);
                        str <<= 2 * exLength;
                        for (unsigned i = 0; i < exLength; ++i) {
                            str.set(exLength - 1 - i, tmpStraight->get(tmpStraight->length - 1 - i));
                        }
                        writeStraight.length += exLength;
                    }
                }
            }
        }
        if (writeStraight.length < nextKmer) continue;
        for (unsigned i = 0; i < nextKmer - 1; ++i) {
            next.setForward(nextKmer - i - 2, str.get(writeStraight.length - 1 - i));
            next.setReverse(i + 1, 0x3 ^ str.get(writeStraight.length - 1 - i));
        }
        for (unsigned long long i = 0; i < writeStraight.length - nextKmer + 1; ++i) {
            next.forward <<= 2;
            next.maskForward(nextMask);
            next.setForward(0, str.get(writeStraight.length - nextKmer - i));
            next.reverse >>= 2;
            next.setReverse(nextKmer - 1, 0x3 ^ str.get(writeStraight.length - nextKmer - i));
            selectedKey = std::min(next.forward,  next.reverse);
            unsigned short occurrence = writeStraight.coverage * occurrenceRatio + 0.5;
            auto it = table.find_any(selectedKey);
            if (it->second < occurrence) {
                it->second = occurrence;
                it->first = selectedKey;
            }
        }
    }

}


template <typename KMER>
void BruijnGraph<KMER>::saveContigSimple(FILE *contigFP, const double coverageRatio) const
{
    if (contigFP == NULL) {
        platanus::makeTemporaryFile();
    }
    for (auto tableIt = leftStraightTable.begin(), tableEnd = leftStraightTable.end(); tableIt != tableEnd; ++tableIt) {
        const Straight &straight = tableIt->second;
        Straight writeStraight;
        if (straight.coverage == UINT16_MAX) {
            continue;
        }

        unsigned long long length = straight.length + this->kmerLength - 1;
		unsigned short coverage = static_cast<unsigned short>(straight.coverage * coverageRatio + 0.5);
        binstr_t seq(length);
        seq.convertFromVector(straight.seq);
        fwrite(&length, sizeof(unsigned long long), 1, contigFP);
        fwrite(&(coverage), sizeof(unsigned short), 1, contigFP);
        seq.writeTemporaryFile(contigFP);
    }
}


template <typename KMER>
void BruijnGraph<KMER>::saveJunction(FILE *contigFP, const double coverageRatio) const
{
    if (contigFP == NULL) {
        platanus::makeTemporaryFile();
    }

	unsigned long long length = this->kmerLength;

    for (auto tableIt = junctionTable.begin(), tableEnd = junctionTable.end(); tableIt != tableEnd; ++tableIt) {
        const Junction &junction = tableIt->second;
        if (junction.coverage == UINT16_MAX) {
            continue;
        }

		unsigned short coverage = static_cast<unsigned short>(junction.coverage * coverageRatio + 0.5);
		KMER kmer(this->kmerLength);
		kmer.forward = tableIt->first;
        binstr_t seq(length);
        KMER(this->kmerLength).convertKeyToBinstr(seq, tableIt->first);

        fwrite(&length, sizeof(unsigned long long), 1, contigFP);
        fwrite(&(coverage), sizeof(unsigned short), 1, contigFP);
        seq.writeTemporaryFile(contigFP);
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// cut branch (less coverage contig)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::cutBranchIterative(const unsigned long long numThread)
{
    unsigned long long totalDelete = 0;
    
    std::cerr << "removing branches..." << std::endl;
    std::cerr << "BRANCH_DELETE_THRESHOLD=" << branchThreshold << std::endl;

    unsigned long long numDelete;
    do {
        numDelete = this->cutBranch(numThread);
        std::cerr << "NUM_CUT=" << numDelete << std::endl;
        concatinateNodes();
        totalDelete += numDelete;
    } while (numDelete > 0);
    std::cerr << "TOTAL_NUM_CUT=" << totalDelete << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// crush bubble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::crushBubbleIterative(const double averageCoverage)
{
    std::cerr << "removing bubbles..." << std::endl;

    bubbleFP = platanus::makeTemporaryFile();

    std::cerr << "BUBBLE_IDENTITY_THRESHOLD=" << bubbleThreshold << std::endl;
    unsigned long long totalCrush = 0;
    unsigned long long numCrush;

    do {
        numCrush = this->crushBubble(averageCoverage);
	std::cerr << "NUM_REMOVED_BUBBLES=" << numCrush << std::endl;
        concatinateNodes();
        totalCrush += numCrush;
    } while (numCrush > 0);
    std::cerr << "TOTAL_NUM_REMOVED_BUBBLES=" << totalCrush << std::endl;

}



//////////////////////////////////////////////////////////////////////////////////////
// save only kmer terminal
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::saveEdgeKmer(Counter<KMER> &counter, const unsigned nextKmerLength) const
{

    const unsigned diff = nextKmerLength - kmerLength;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));

    if (counter.kmerFP != NULL)
        fclose(counter.kmerFP);
    counter.kmerFP = platanus::makeTemporaryFile();
    KMER leftKey(kmerLength);
    KMER rightKey(kmerLength);

    auto tableEnd = leftStraightTable.end();
    for (auto tableIt = leftStraightTable.begin(); tableIt != tableEnd; ++tableIt) {
        const Straight &straight = tableIt->second;
        if (straight.coverage == UINT16_MAX)
            continue;

        // if straight length is too short (all seq contains edge)
        if (straight.length < diff * 2) {
            for (unsigned i = 0; i < kmerLength - 1; ++i) {
                leftKey.setForward(kmerLength - 2 - i, straight.get(straight.length + kmerLength - 2 - i));
                leftKey.setReverse(i + 1, 0x3 ^ straight.get(straight.length + kmerLength - 2 - i));
            }
            // write all seq divided kmer in temporary files
            for (unsigned i = 0; i < straight.length; ++i) {
                leftKey.forward <<= 2;
                leftKey.maskForward(mask);
                leftKey.setForward(0, straight.get(straight.length - 1 - i));
                leftKey.reverse >>= 2;
                leftKey.setReverse(kmerLength - 1, 0x3 ^ straight.get(straight.length - 1 - i));
                if (straight.coverage != 0)
                    counter.setOccurrenceValue(std::min(leftKey.forward, leftKey.reverse), straight.coverage);
            }
        } else {
            for (unsigned i = 0; i < kmerLength - 1; ++i) {
                leftKey.setForward(kmerLength - 2 - i, straight.get(straight.length + kmerLength - 2 - i));
                leftKey.setReverse(i + 1, 0x3 ^ straight.get(straight.length + kmerLength - 2 - i));
                rightKey.setForward(kmerLength - 2 -i, straight.get(kmerLength + diff - 2 - i));
                rightKey.setReverse(i + 1, 0x3 ^ straight.get(kmerLength + diff - 2 - i));
            }
            for (unsigned i = 0; i < diff; ++i) {
                // write left seq divided kmer in temporary files
                leftKey.forward <<= 2;
                leftKey.maskForward(mask);
                leftKey.setForward(0, straight.get(straight.length - 1 - i));
                leftKey.reverse >>= 2;
                leftKey.setReverse(kmerLength - 1, 0x3 ^ straight.get(straight.length - 1 - i));
                if (straight.coverage != 0)
                    counter.setOccurrenceValue(std::min(leftKey.forward, leftKey.reverse), straight.coverage);
                // write right seq divided kmer in temporary files
                rightKey.forward <<= 2;
                rightKey.maskForward(mask);
                rightKey.setForward(0, straight.get(diff - 1 - i));
                rightKey.reverse >>= 2;
                rightKey.setReverse(kmerLength - 1, 0x3 ^ straight.get(diff - 1 - i));
                if (straight.coverage != 0)
                    counter.setOccurrenceValue(std::min(rightKey.forward, rightKey.reverse), straight.coverage);
            }
        }
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// calculate length cutoff, coverage cutoff and average coverage value from distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::calcStraightStats(unsigned long long &lengthCut, unsigned long long &coverageCut, double &averageCoverage) const
{
    unsigned long long sum = 0;
    unsigned long long num = 0;
    std::vector<unsigned long long> lengthDistribution;
    std::vector<unsigned long long> coverageDistribution(UINT16_MAX);
    auto tableEnd = leftStraightTable.end();
    for (auto tableIt = leftStraightTable.begin(); tableIt != tableEnd; ++tableIt) {
        const Straight &straight = tableIt->second;
        if (straight.coverage == UINT16_MAX)
            continue;

        if (straight.length + 1 > lengthDistribution.size()) {
            lengthDistribution.resize(straight.length + 1, 0);
        }
        ++lengthDistribution[straight.length];
        coverageDistribution[straight.coverage] += straight.length;
    }

    coverageCut = getLeftMinimal(coverageDistribution);
    lengthCut = getLeftMinimal(lengthDistribution);

    tableEnd = leftStraightTable.end();
    for (auto tableIt = leftStraightTable.begin(); tableIt != tableEnd; ++tableIt) {
        const Straight &straight = tableIt->second;
        if (straight.coverage == UINT16_MAX || (straight.coverage < coverageCut && straight.length < lengthCut))
            continue;
        sum += straight.coverage * straight.length;
        num += straight.length;
    }
    averageCoverage = static_cast<double>(sum) / num;
}


template <typename KMER>
unsigned long long BruijnGraph<KMER>::getLeftMinimalCoverage() const
{
    std::vector<unsigned long long> coverageDistribution(UINT16_MAX);
    auto tableEnd = leftStraightTable.end();
    for (auto tableIt = leftStraightTable.begin(); tableIt != tableEnd; ++tableIt) {
        const Straight &straight = tableIt->second;
        if (straight.coverage != UINT16_MAX)
			coverageDistribution[straight.coverage] += straight.length;
    }

    return getLeftMinimal(coverageDistribution);
}


//////////////////////////////////////////////////////////////////////////////////////
// get local minimal value in distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline unsigned long long BruijnGraph<KMER>::getLeftMinimal(const std::vector<unsigned long long> &dist) const
{
    unsigned long long i;
    unsigned long long pre;
    unsigned long long size = dist.size();

    if (size == 0) return 0;

    for (i = 0; i < size; ++i)
        if (dist[i] > 0) break;

    pre = dist[i];

    for (++i; i < size; ++i) {
        if (dist[i] >= pre)
            break;
        pre = dist[i];
    }
    if (i < size)
        return i - 1;
    else
        return 1;
}





//////////////////////////////////////////////////////////////////////////////////////
// delete too small coverage or length contig
// when delete, this seq terminal junction's out value change. (delete this seq in junction)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long BruijnGraph<KMER>::deleteErroneousStraightNode(const unsigned long long lengthCut, const unsigned long long coverageCut, const unsigned long long numThread)
{
    unsigned long long numDelete = 0;
    omp_lock_t lock;
    omp_init_lock(&lock);


    KMER key(kmerLength);
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1) firstprivate(key) reduction(+: numDelete)
    for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
        unsigned long long i = 0;
        for (auto tableIt = leftStraightTable.begin(), tableEnd = leftStraightTable.end(); tableIt != tableEnd; ++tableIt, ++i) {
            if (i % numThread != threadID) continue;
            const Straight &straight = tableIt->second;
            if (straight.coverage == UINT16_MAX || straight.coverage >= coverageCut || straight.length >= lengthCut) {
//            if (straight.coverage == UINT16_MAX || straight.coverage >= coverageCut || straight.length >= lengthCut || ((straight.out >> 4) && (straight.out & 0xf))) {
                continue;
            }
            if (straight.out >> 4) {
                key.forward = tableIt->first;
                key.forward >>= 2;
                key.setForward(kmerLength - 1, platanus::Flag2Base(straight.out >> 4));
                Junction *junction = findJunction(key.forward);
                if (junction != NULL) {
                    omp_set_lock(&lock);
                    junction->out ^= 1 << straight.get(straight.length - 1);
                    omp_unset_lock(&lock);
                }
            }

            if (straight.out & 0xf) {
                for (unsigned i = 0; i < kmerLength - 1; ++i)
                    key.setForward(i + 1, straight.get(i));
                key.setForward(0, platanus::Flag2Base(straight.out & 0xf));
                Junction *junction = findJunction(key.forward);
                if (junction != NULL) {
                    omp_set_lock(&lock);
                    junction->out ^= 1 << (4 + straight.get(kmerLength - 1));
                    omp_unset_lock(&lock);
                }
            }
            deleteStraight(straight, key);
            ++numDelete;
        }
    }
    omp_destroy_lock(&lock);

    return numDelete;
}




//////////////////////////////////////////////////////////////////////////////////////
// delete erroneousnode
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::deleteErroneousStraightNodeIterative(const unsigned long long lengthCut, const unsigned long long coverageCut, const unsigned long long numThread)
{
    unsigned long long total = 0;
    unsigned long long numDelete = 0;

    std::cerr << "removing erroneous nodes..." << std::endl;

    do {
        numDelete = cutBranch(numThread);
        numDelete += deleteErroneousStraightNode(lengthCut, coverageCut, numThread);
	std::cerr << "NUM_REMOVED_NODES=" << numDelete << std::endl;
        concatinateNodes();
        total += numDelete;
    } while (numDelete > 0);
    std::cerr << "TOTAL_NUM_REMOVED_NODES=" << total << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// make bruijn graph for gap close
// difference with this->makeInitBruijnGraph
// 1. exist occurence threshold(minCoverage)
// 2. think only one side direction (need not consider reverse strand)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void BruijnGraph<KMER>::makeBruijnGraphForGapClose(Counter<KMER> &counter, const long minCoverage)
{
    unsigned char leftFlags, rightFlags;
    unsigned long long sum = 0;
    std::vector<char> leftSeq, rightSeq;
    std::vector<typename KMER::keyType> sorter;
    std::stack<typename KMER::keyType> kmerStack;
    Junction leftJunction;
    Straight leftStraight;
    KMER kmer(this->kmerLength);
    KMER startKmer(this->kmerLength);
    KMER leftKmer(this->kmerLength);
    KMER rightKmer(this->kmerLength);
    typename KMER::keyType key(kmerLength);
    typename KMER::keyType nextKey(kmerLength), tmpKey(kmerLength);
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));

    this->deleteAllTable();
    int ss = 0;
    counter.assignOccurrenceKeysGap(sorter);

    for (auto it = sorter.begin(), end = sorter.end(); it != end; ++it) {
        if (counter.findValueGap(*it) < minCoverage) continue;
            ++ss;
        kmerStack.push(*it);

        if (kmerStack.size() > 0) {
            startKmer.forward = kmerStack.top();
            unsigned short kmerValue = counter.findValueAndOverwriteGap(startKmer.forward, UINT16_MAX);
            if (kmerValue == UINT16_MAX) {
                kmerStack.pop();
                continue;
            }
            sum = kmerValue;
            // search parent kmer node 
            // example kmer "ATTGCAGC", find "*ATTGCAG"
            leftKmer.forward = startKmer.forward;
            leftKmer.forward >>= 2;
            leftFlags = 0;
            leftSeq.clear();
            for (unsigned char base = 0; base < 4; ++base) {
                leftKmer.setForward(leftKmer.kmerLength - 1, base);
                if (counter.findValueGap(leftKmer.forward) >= minCoverage) {
                    leftFlags |= 1 << base;
                }
            }
            if (leftFlags != 0)
                leftSeq.push_back(platanus::Flag2Base(leftFlags));
           // search parent kmer node
           // example kmer "ATTGCAGC", find "TTGCAGC*"
            rightKmer.forward = startKmer.forward;
            rightKmer.forward <<= 2;
            rightKmer.maskForward(mask);
            rightSeq.clear();
            // search whether parent kmer node contains valueTable.
            rightFlags = 0;
            for (unsigned char base = 0; base < 4; ++base) {
                rightKmer.setForward(0, base);
                if (counter.findValueGap(rightKmer.forward) >= minCoverage) {
                    rightFlags |= 1 << base;
                }
            }
            if (rightFlags != 0)
                rightSeq.push_back(platanus::Flag2Base(rightFlags));

            if ((leftSeq.size() > 0 && leftSeq[0] == 4) || (rightSeq.size() > 0 && rightSeq[0] == 4)) {
                kmerStack.pop();

                for (unsigned char base = 0; base < 4; ++base) {
                    if ((rightFlags & (1 << base)) != 0) {
                        rightKmer.setForward(0, base);
                        kmerStack.push(rightKmer.forward);
                    }
                    if ((leftFlags & (1 << base)) != 0) {
                        leftKmer.setForward(leftKmer.kmerLength - 1, base);
                        kmerStack.push(leftKmer.forward);
                    }
                }

                leftJunction.coverage = sum;
                leftJunction.out = (leftFlags << 4) | rightFlags;
                setJunctionTable(leftJunction, startKmer.forward);
                continue;
            }


            if (rightSeq.size() > 0) {
                // set kmer key
                kmer.forward = startKmer.forward;
                kmer.forward <<= 2;
                kmer.maskForward(mask);
                kmer.setForward(0, rightSeq[0]);
                key = kmer.forward;
                while (1) {
                    // check whether only one node is found in reverse relation.
                    // i.e. we check one son "ATGGCAG" in before phase
                    // ATGGCAG -> AATGGCA, check AATGGCA <- only ATGGCAG, not exist ATGGCAA, ATGGCAC, ATGGCAT
                    rightKmer = kmer;
                    kmer.forward >>= 2;
                    unsigned char tmpFlags = 0;
                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(kmer.kmerLength - 1, base);
                        if (counter.findValueGap(kmer.forward) >= minCoverage) {
                            tmpFlags |= 1 << base;
                        }
                    }
                    if (platanus::FlagCount(tmpFlags) > 1 || counter.findValueGap(key) == UINT16_MAX) {
                        rightFlags = 0 | (1 << rightSeq.back());
                        rightSeq.pop_back();
                        break;
                    }

                    // search parent kmer node "iterative"
                    // example kmer "ATTGCAGC", find "*ATTGCAG"
                    kmer.forward <<= 4;
                    kmer.maskForward(mask);
                    kmer.setForward(1, rightSeq.back());
                    rightFlags = 0;

                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(0, base);
                        tmpKey = kmer.forward;
                        kmerValue = counter.findValueGap(tmpKey);
                        if (kmerValue >= minCoverage) {
                            rightFlags |= 1 << base;
                            nextKey = tmpKey;
                        }
                    }

                    if (rightFlags == 0) {
                        sum += counter.findValueAndOverwriteGap(key, UINT16_MAX);
                        break;
                    } else if (platanus::FlagCount(rightFlags) > 1) {
                        rightFlags = 0 | (1 << rightSeq.back());
                        rightSeq.pop_back();
                        break;
                    }

                    rightSeq.push_back(platanus::Flag2Base(rightFlags));
                    sum += counter.findValueAndOverwriteGap(key, UINT16_MAX);
                    key = nextKey;
                }
            }


            if (leftSeq.size() > 0) {
                kmer = startKmer;
                kmer.forward >>= 2;
                kmer.setForward(kmer.kmerLength - 1, leftSeq[0]);
                key = kmer.forward;
                while (1) {
                    leftKmer = kmer;
                    kmer.forward <<= 2;
                    kmer.maskForward(mask);
                    unsigned char tmpFlags = 0;
                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(0, base);
                        if (counter.findValueGap(kmer.forward) >= minCoverage) {
                            tmpFlags |= 1 << base;
                        }
                    }
                    if (platanus::FlagCount(tmpFlags) > 1 || counter.findValueGap(key) == UINT16_MAX) {
                        leftFlags = 0 | (1 << leftSeq.back());
                        leftSeq.pop_back();
                        break;
                    }
                    kmer.forward >>= 4;
                    kmer.setForward(kmer.kmerLength - 2, leftSeq.back());
                    leftFlags = 0;

                    for (unsigned char base = 0; base < 4; ++base) {
                        kmer.setForward(kmer.kmerLength - 1, base);
                        tmpKey = kmer.forward;
                        kmerValue = counter.findValueGap(tmpKey);
                        if (kmerValue >= minCoverage) {
                            leftFlags |= 1 << base;
                            nextKey = tmpKey;
                        }
                    }
                    if (leftFlags == 0) {
                        sum += counter.findValueAndOverwriteGap(key, UINT16_MAX);
                        break;
                    } else if (platanus::FlagCount(leftFlags) > 1) {
                        leftFlags = 0 | (1 << leftSeq.back());
                        leftSeq.pop_back();
                        break;
                    }

                    leftSeq.push_back(platanus::Flag2Base(leftFlags));
                    sum += counter.findValueAndOverwriteGap(key, UINT16_MAX);
                    key = nextKey;
                }
            }

            kmerStack.pop();

            if (leftFlags != 0) {
                leftKmer.setForward(leftKmer.kmerLength - 1, platanus::Flag2Base(leftFlags));
                kmerStack.push(leftKmer.forward);
            }

            if (rightFlags != 0) {
                rightKmer.setForward(0, platanus::Flag2Base(rightFlags));
                kmerStack.push(rightKmer.forward);
            }


            leftStraight.length = leftSeq.size() + rightSeq.size() + 1;
            binstr_t seq(leftStraight.length + kmerLength);
            for (unsigned j = 0, n = leftSeq.size(); j < n; ++j) {
                seq.set(leftSeq.size() - 1 - j, leftSeq[leftSeq.size() - 1 - j]);
            }
            seq <<= 2 * kmerLength;
            for (unsigned j = 0; j < kmerLength; ++j) {
                seq.set(kmerLength - 1 - j, (startKmer.forward >> (2 * (kmerLength - 1 - j))) & 0x3);
            }
            seq <<= 2 * rightSeq.size();
            for (unsigned j = 0, n = rightSeq.size(); j < n; ++j) {
                seq.set(rightSeq.size() - 1 - j, rightSeq[j]);
            }
            leftStraight.coverage = (static_cast<double>(sum) / static_cast<double>(leftStraight.length) + 0.5);
            leftStraight.out = (leftFlags << 4) | rightFlags;
            seq.convertToVector(leftStraight.seq);
            this->insertStraight(leftStraight);

        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// divide straight node into kmer
// this kmer length is higher one for gap close
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
template <typename NEXT>
void BruijnGraph<KMER>::saveLargeKmerForGapClose(Counter<NEXT> &counter)
{
    unsigned nextLength = counter.getKmerLength();
    typename NEXT::keyType key(nextLength);
    for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
        const Straight &straight = it->second;
        if (straight.coverage == UINT16_MAX || straight.length + this->kmerLength - 1 < nextLength) continue;

        for (unsigned i = 0; i < nextLength - 1; ++i) {
            key.set(i + 1, straight.get(i));
        }
        for (unsigned i = 0; i < straight.length + this->kmerLength - nextLength; ++i) {
            key >>= 2;
            key.set(nextLength - 1, straight.get(nextLength - 1 + i));
            counter.setOccurrenceValueGap(key, straight.coverage);
        }
    }
}



template <typename KMER>
double BruijnGraph<KMER>::getAverageCoverageExcludingBubble(void)
{
    typedef unsigned long long u64_t;

    std::vector<Junction*> junction(4);
    std::vector<Straight*> straight(4);
    Junction *rootJunction;
	KMER kmer(kmerLength);
    typename KMER::keyType key(kmerLength);
    std::unordered_map<typename KMER::keyType, bool, typename KMER::hasher> isBubble;

    for (auto it = junctionTable.begin(), end = junctionTable.end(); it != end; ++it) {
        rootJunction = &(it->second);
        if (rootJunction->coverage == UINT16_MAX || platanus::FlagCount(rootJunction->out & 0x0f) < 2) continue;

        searchBubbleStructure(rootJunction, junction, straight, it->first);
        for (unsigned char base1 = 0; base1 < 3; ++base1) {
            if (straight[base1] == NULL) continue;
            for (unsigned char base2 = base1 + 1; base2 < 4; ++base2) {
                if (straight[base2] == NULL || straight[base1]->coverage == UINT16_MAX || straight[base2]->coverage == UINT16_MAX || junction[base1] == NULL || junction[base1] != junction[base2]) continue;

				for (unsigned i = 0; i < kmerLength; ++i)
					kmer.setForward(i, straight[base1]->get(i));
				key = kmer.forward;
				isBubble[key] = true;

				for (unsigned i = 0; i < kmerLength; ++i)
					kmer.setForward(i, straight[base2]->get(i));
				key = kmer.forward;
				isBubble[key] = true;
			}
        }
    }

	u64_t sum = 0;
	u64_t num = 0;

	for (auto it = leftStraightTable.begin(), end = leftStraightTable.end(); it != end; ++it) {
		if (it->second.coverage == UINT16_MAX) continue;
		for (unsigned i = 0; i < kmerLength; ++i)
			kmer.setForward(i, it->second.get(i));
		key = kmer.forward;
		if (isBubble.find(key) == isBubble.end()) {
			num += it->second.length;
			sum += it->second.length * it->second.coverage;
		}
	}

    for (auto it = junctionTable.begin(), end = junctionTable.end(); it != end; ++it) {
		if (it->second.coverage != UINT16_MAX) {
			++num;
			sum += it->second.coverage;
		}
	}

    return static_cast<double>(sum) / static_cast<double>(num);
}








#endif
