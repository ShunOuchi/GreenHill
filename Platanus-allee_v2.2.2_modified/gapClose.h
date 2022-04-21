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

#ifndef GAP_CLOSE_H
#define GAP_CLOSE_H

#include "baseCommand.h"
#include "common.h"
#include "seqlib.h"
#include "mapper.h"
#include <vector>
#include <string>
#include <set>



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// GAP class
// this class has GAP information
struct GAP
{
    enum STATE {CLOSED, PART_CLOSED, UNCLOSED};

    // constructor and destructor
    GAP(): scaffoldID(-1), start(0), end(0), length(0), closingLength(0), remainedGap(0), totalSeqBase(0), seq(), headSeq(), tailSeq(), closingSeq(), state(UNCLOSED), tooManyReads(false) {}
    GAP(const GAP &) = default;
    GAP &operator=(const GAP &) = default;
    ~GAP() = default;


    int scaffoldID;
    int start;
    int end;
    int length;
    int closingLength;
    int remainedGap;
    unsigned long long totalSeqBase;
    std::vector<std::vector<char> > seq;
    std::vector<char> headSeq;
    std::vector<char> tailSeq;
    std::vector<char> closingSeq;
    STATE state;
    bool tooManyReads;

// < operator assign > operator's mean on purpose
    bool operator()(GAP *gapLeft, GAP *gapRight)
    {   return gapLeft->seq.size() > gapRight->seq.size(); }

    struct seqGreater
    {
        bool operator()(const std::vector<char> &seq1, const std::vector<char> &seq2)
        {
            if (seq1.size() != seq2.size()) return seq1.size() < seq2.size();

            for (unsigned long long i = 0; i < seq1.size(); ++i) {
                if (seq1[i] != seq2[i]) return seq1[i] < seq2[i];
            }
            return 0;
        }
    };

    void sortSeq(void)
    {
		if (seq.size() > 3)
			std::sort(seq.begin() + 2, seq.end(), seqGreater());
    }
};
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
// gap close execute class
//////////////////////////////////////////////////////////////////////////////////////
class GapClose : public BaseCommand
{
private:
    static const unsigned short SD_RATIO_MAPPED_GAP;
	static const int MIN_OVERLAP_FOR_CIRCLE;

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Gloser class
// this class exec gap close actually
// maybe this class doesn't need because this class can be integrated GapClose Class
    class Closer
    {
    private:

		struct GapClosedSeq
		{
			std::string seq;
			std::string baseName;
			int coverage;
			int leftExtendedLen;
			int rightExtendedLen;
			int heteroCounterPartIndex;
			bool circularFlag;
			bool heteroFlag;
			bool redundantFlag;

			GapClosedSeq(): coverage(0), leftExtendedLen(0), rightExtendedLen(0), heteroCounterPartIndex(0), circularFlag(false), heteroFlag(false), redundantFlag(false) {}
			~GapClosedSeq() = default;
		};

        static const unsigned short MIN_NUM_READS_COVERING_SMALL_GAP;
        static const unsigned BRUIJN_MIN_KMER;
        static const unsigned BRUIJN_MAX_KMER;
        static const unsigned MIN_COVERAGE;
        static const unsigned short SD_RATIO_TOO_MANY_READS_ON_GAP;

        const long numScaffold;
        const std::vector<platanus::SEQ> scaffold;
        const std::vector<unsigned short> scaffoldCoverage;
        const std::vector<std::string> scaffoldName;
        std::vector<GAP> gap;
        std::unordered_map<unsigned long long, long> gapTable;
		std::vector<std::string> scaffoldFileName;
		std::vector<long> numScaffoldInFile;

        const unsigned long long olcThreshold;
        const unsigned minOverlapOLC;
        const unsigned minOverlapDe;
        const unsigned maxEditDistance;
        const double maxMissRate;
        const double minConsensusRate;

		std::vector<GapClosedSeq> gapClosedSeq;

        void insertGap(const int id, const int offset, const long value);
        long findGapID(const int id, const int offset);
        void judgePairReadMappedNearGap(platanus::Position &mappedPosition, platanus::SEQ &pairSeq, const long averageInsSize, const long tolerance, FILE *tmpFP);
        void clearGap(unsigned gapID);
        double tooManyReadsGathered(const GAP * const gap, const long sd) const;
		long selfOverlap(std::string &seq, const long minOverlap);

    public:
        // constructor and destructor
        Closer(const platanus::Contig &contig, const unsigned long long threshold, const unsigned minOLC, const unsigned minDe, const unsigned ed, const double miss, const double con): numScaffold(contig.numSeq), scaffold(contig.seq), scaffoldCoverage(contig.coverage), scaffoldName(contig.name), gap(), gapTable(), olcThreshold(threshold), minOverlapOLC(minOLC), minOverlapDe(minDe), maxEditDistance(ed), maxMissRate(miss), minConsensusRate(con) {}
        ~Closer() = default;
        Closer() = delete;
        Closer(const Closer &) = delete;
        Closer &operator=(const Closer &) = delete;


        void makeGapTable(void);
        FILE *saveGapCoveringReads(const std::vector<SeqLib> &library, const long numThread);
        void loadLocalReads(FILE *gapSeqFP);
        void gapCloseUsingPairReads(const long keyLength, const long sd, const long numThread, const bool noPartialFlag, const bool final=false);
        void saveUnusedReads(const unsigned long long maxBase, FILE *unusedFP);
        void loadUnusedReads(FILE *unusedFP);
        void closeSmallGaps(const double consensusRate, FILE *gapSeqFP);
        bool decideConsensusFromReads(const std::vector<std::pair<std::vector<char>, int> > &seq, const double threshold, std::pair<std::vector<char>, int> &consensus, const unsigned long long id);
        bool decideConsensusFromReads(const std::multiset<std::pair<std::vector<char>, int> > &seq, const double threshold, std::pair<std::vector<char>, int> &consensus, const unsigned long long id);
		void generateGapClosedSeq(const bool addN, const long contigMaxK);
		void printGapClosedSeq(const std::string &outPrefix, const unsigned long long readLength, const unsigned long long contigMaxK, const bool extentionFlag);
		void setHeteroInfoOfGapClosedSeq();
		void findCircularGapClosedSeq(const long numThread);
        void printSeq(const std::string &outName, const unsigned long long readLength, const bool addN);
        void removeEdgeNInformation(void);
        unsigned long long calcMaxGapSeqBase(const long numThread) const;
        FILE *localAssemble(const unsigned long long numThread);
		void setFastaFileInfo(const std::vector<std::string> fastaFileName);
		void markRedundantSeq(const long maxPrefixLength, const unsigned long long numThread);
        void printContig(const std::string &outputFilename, FILE *contigFP, const double averageReadLength, const unsigned long long readLength, const unsigned long long contigMaxK) const
        {
            platanus::printContig(outputFilename, contigFP, averageReadLength / (averageReadLength - BRUIJN_MAX_KMER + 1.0), readLength, contigMaxK, "scaffold");
        }
    };
//////////////////////////////////////////////////////////////////////////////////////


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].size() == 0) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        if (optionMultiArgs["-f"].size() == 0 && optionPairFile.size() == 0) {
            std::cerr << "Error: not specified read file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    static const unsigned short HEAD_TAIL_SEQ_LEN;
    static const unsigned short GAP_CANDIDATE_SEQ_1ST;


    // constructor and destructor
    GapClose();
    ~GapClose() = default;
    GapClose(const GapClose &) = delete;
    GapClose &operator=(const GapClose &) = delete;

    virtual void usage(void) const;
    virtual void exec(void);

    void readLibrary(Mapper &mapper, platanus::Contig &contig, std::vector<std::vector<SeqLib> > &libraryMT, std::vector<SeqLib> &singleLibrary, std::vector<int> &numFilePerLibraryID, const int numLibrary, const int numThread);
};

#endif
