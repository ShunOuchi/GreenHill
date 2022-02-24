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

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#include "baseCommand.h"
#include "counter.h"
#include "graph.h"


//////////////////////////////////////////////////////////////////////////////////////
// assemble class
//////////////////////////////////////////////////////////////////////////////////////
class Assemble : public BaseCommand
{
private:
    // constant parameter
    static const unsigned SMOOTHING_WINDOW;
    static const double MAX_COVERAGE_CUT_DIFF_RATE;
	static const double REPEAT_MODE_CUTOFF_FACTOR;
	static const double REPEAT_MODE_BUBBLE_IDENTITY_THRESHOLD;
	static const int MAX_COVERAGE_CUTOFF_FACTOR;


	bool repeatModeFlag;
    unsigned long long lengthStep;
    unsigned long long minCoverage;
    double averageLength;
    double averageCoverage;
    double minLogPJoin;
    double maxKmerLengthRatio;
    unsigned long long numThread;
    std::vector<unsigned long long> coverageCutoff;
    unsigned long long lengthCutoff;
    unsigned long long numKmer;
    unsigned long long memory;
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
    unsigned long long doubleHashSize;



    void initializeParameters(void);
    template <typename KMER> void initialKmerAssemble(std::vector<unsigned> &k, FILE **readFP, const unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const double coverageCutoffFactor);
    template <typename KMER> void saveAndRedoAssemble(const unsigned k, const unsigned previousk, FILE **readFP, FILE **allReadFP, const unsigned long long numThread, BruijnGraph<KMER> &lowKmerGraph, Counter<KMER> &lowKmerCounter, const unsigned long long position);
    template <typename KMER> void saveGraph(const unsigned k, FILE **readFP, const unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const bool readFPCloseFlag);
    template <typename KMER> void redoAssemble(const unsigned k, const unsigned previousk, FILE **readFP, FILE **allReadFP, unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const unsigned long long position);
    template <typename KMER> void outputAndAfterTreatment(const unsigned k, FILE ** readFP, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const unsigned long long numThread);
	template <typename KMER> void mergeContig(platanus::Contig &contig, Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength, FILE **mergedContigFP);
    unsigned long long extendKmer(const double minLogPJoin, const double averageCoverage, const double averageLength, const unsigned long long minCoverage, std::vector<unsigned> &kmer, const unsigned long long lengthStep);
    void showKmerExtension(const double minLogPJoin, const double averageCoverage, const double averageLength, const unsigned long long minCoverage, unsigned long long kmer, const unsigned long long lengthStep, unsigned long long coverageCutoff);
    double calcLogProbabilityJoin(const unsigned long long coverageCutoff, const double averageCoverage, const double averageLength, const unsigned long long largeKmer, const unsigned long long smallKmer) const;
    unsigned long long decreaseCoverageCutoff(const unsigned long long coverageCutoff, const double averageCoverage, const double averageLength, const double minLogPJoin, const unsigned long long largeKmer, const unsigned long long smallKmer) const;
    unsigned long long calcMaxKmerLength(const double minimumLogPJoin, const double averageCoverage, const double averageLength, const unsigned long long minimumCoverage, unsigned long long kmerLength, const unsigned long long lengthStep) const;
    template <typename KEY, typename KMER> void  updateDoubleHashSize(const BruijnGraph<KMER> &graph, const unsigned k, const unsigned previousk);
    void readInputFile(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;
    void readFasta(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;
    void readFastaUncompressed(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;
    void readFastaCompressed(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;
    void readFastq(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;
    void readFastqUncompressed(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;
    void readFastqCompressed(const std::string &inputFilename, FILE **outputMT, const unsigned long long numThread) const;


    void initializeGraph(const double bubble, const double branch)
    {
        kmer31Graph.setBubbleAndBranch(bubble, branch);
        kmer63Graph.setBubbleAndBranch(bubble, branch);
        kmer95Graph.setBubbleAndBranch(bubble, branch);
        kmer127Graph.setBubbleAndBranch(bubble, branch);
        kmer159Graph.setBubbleAndBranch(bubble, branch);
        kmerNGraph.setBubbleAndBranch(bubble, branch);
    }

    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-f"].size() == 0) {
            std::cerr << "Error: not specified input file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    Assemble();
    Assemble(const Assemble &) = delete;
    Assemble &operator=(const Assemble &) = delete;
    ~Assemble() {}
    virtual void usage() const;
    virtual void exec();


};



#endif
