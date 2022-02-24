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

#include "assemble.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <climits>
#include <cfloat>

using std::string;
using std::vector;
using std::unordered_map;
using std::cerr;
using std::endl;
using std::ifstream;
using platanus::SEQ;

//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const unsigned Assemble::SMOOTHING_WINDOW = 1;
const double Assemble::MAX_COVERAGE_CUT_DIFF_RATE = 0.25;
const double Assemble::REPEAT_MODE_CUTOFF_FACTOR = 1.75;
const double Assemble::REPEAT_MODE_BUBBLE_IDENTITY_THRESHOLD = 0.95;
const int Assemble::MAX_COVERAGE_CUTOFF_FACTOR = 0;


//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
Assemble::Assemble()
: lengthStep(0), minCoverage(0), averageLength(0), averageCoverage(0), minLogPJoin(0), numThread(0), coverageCutoff(), lengthCutoff(0), numKmer(0), memory(0), kmer31Graph(), kmer63Graph(), kmer95Graph(), kmer127Graph(), kmer159Graph(), kmerNGraph(), kmer31Counter(), kmer63Counter(), kmer95Counter(), kmer127Counter(), kmer159Counter(), kmerNCounter(), doubleHashSize()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-k"] = "32";
    optionSingleArgs["-K"] = "0.5";
    optionSingleArgs["-s"] = "20";
    optionSingleArgs["-n"] = "0";
    optionSingleArgs["-c"] = "2";
    optionSingleArgs["-a"] = "10.0";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-d"] = "0.5";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-m"] = "16";
    optionMultiArgs["-f"] = vector<string>();
    optionSingleArgs["-tmp"] = ".";
    optionBool["-repeat"] = false;
}



//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::usage(void) const
{
    std::cerr << "\nUsage platanus_allee assemble [Options]\n"
              << "Options:\n"
              << "    -o STR               : prefix of output files (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -f FILE1 [FILE2 ...] : reads file (fasta or fastq, number <= "<< platanus::ConstParam::MAX_FILE_NUM << ")\n"
              << "    -k INT               : initial k-mer size (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -K FLOAT             : maximum-k-mer factor (maximum-k = FLOAT*read-length, default  " << optionSingleArgs.at("-K") << ")\n"
              << "    -s INT               : step size of k-mer extension (>= 1, default " << optionSingleArgs.at("-s") << ")\n"
              << "    -n INT               : initial k-mer coverage cutoff (default " << optionSingleArgs.at("-n") << ", 0 means auto)\n"
              << "    -c INT               : minimun k-mer coverage (default " << optionSingleArgs.at("-c") << ")\n"
              << "    -a FLOAT             : k-mer extension safety level (default " << optionSingleArgs.at("-a") << ")\n"
//              << "    -u FLOAT             : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -d FLOAT             : maximum difference for branch cutting (coverage ratio, default " << optionSingleArgs.at("-d") << ")\n"
              << "    -e FLOAT             : k-mer coverage depth (k = initial k-mer size specified by -k) of homozygous region (default auto)\n"
              << "    -t INT               : number of threads (<= " << platanus::ConstParam::MAX_THREAD << ", default " << optionSingleArgs.at("-t") << ")\n"
              << "    -m INT               : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "    -tmp DIR             : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
//              << "    -repeat              : mode to assemble repetitive sequences (e.g. 16s rRNA) (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -f option.\n"
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_contig.fa (final)\n"
              << "    PREFIX_kmerFrq.tsv\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initilaize parametors
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::initializeParameters(void)
{
	this->repeatModeFlag = optionBool["-repeat"];
    this->minLogPJoin = log(1.0 - pow(10.0, -(atof(optionSingleArgs["-a"].c_str()))));
    this->numThread = atoi(optionSingleArgs["-t"].c_str());
    this->lengthStep = atoi(optionSingleArgs["-s"].c_str());
    this->maxKmerLengthRatio = atof(optionSingleArgs["-K"].c_str());
    this->minCoverage = atoi(optionSingleArgs["-c"].c_str());
    this->memory = atoi(optionSingleArgs["-m"].c_str()) * static_cast<unsigned long long>(1000000000);

	double bubble = atof(optionSingleArgs["-u"].c_str());
	if (this->repeatModeFlag)
		bubble = REPEAT_MODE_BUBBLE_IDENTITY_THRESHOLD;

    const double branch = atof(optionSingleArgs["-d"].c_str());
    initializeGraph(bubble, branch);
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}

//////////////////////////////////////////////////////////////////////////////////////
// exec assemble
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::exec(void)
{
    initializeParameters();
    vector<unsigned> k;
    FILE *readFP[numThread];
    FILE *allReadFP[numThread];
    string outputFilename;

    omp_set_num_threads(numThread);

	// read file
	for (unsigned long long i = 0; i < numThread; ++i)
		readFP[i] = platanus::makeTemporaryFile();
	for (unsigned i = 0; i < optionMultiArgs["-f"].size(); ++i)
		readInputFile(optionMultiArgs["-f"][i], readFP, numThread);
	for (unsigned long long i = 0; i < numThread; ++i)
		allReadFP[i] = readFP[i];

	k.clear();
	k.push_back(atoi(optionSingleArgs["-k"].c_str()));

	for (unsigned long long i = 0; i < numThread; ++i)
		readFP[i] = allReadFP[i];

	// execute initial assemble
	if (k[0] <= 32)
		initialKmerAssemble(k, readFP, numThread, kmer31Graph, kmer31Counter, 0);
	else if (k[0] <= 64)
		initialKmerAssemble(k, readFP, numThread, kmer63Graph, kmer63Counter, 0);
	else if (k[0] <= 96)
		initialKmerAssemble(k, readFP, numThread, kmer95Graph, kmer95Counter, 0);
	else if (k[0] <= 128)
		initialKmerAssemble(k, readFP, numThread, kmer127Graph, kmer127Counter, 0);
	else if (k[0] <= 160)
		initialKmerAssemble(k, readFP, numThread, kmer159Graph, kmer159Counter, 0);
	else
		initialKmerAssemble(k, readFP, numThread, kmerNGraph, kmerNCounter, 0);

	// kmer extension and assemble iterative
	for (unsigned i = 1; i < numKmer; ++i) {
		if (k[i - 1] <= 32)
			saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer31Graph, kmer31Counter, i);
		else if (k[i - 1] <= 64)
			saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer63Graph, kmer63Counter, i);
		else if (k[i - 1] <= 96)
			saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer95Graph, kmer95Counter, i);
		else if (k[i - 1] <= 128)
			saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer127Graph, kmer127Counter, i);
		else if (k[i - 1] <= 160)
			saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmer159Graph, kmer159Counter, i);
		else
			saveAndRedoAssemble(k[i], k[i - 1], readFP, allReadFP, numThread, kmerNGraph, kmerNCounter, i);
	}

	if (k[numKmer - 1] <= 32) {
		outputAndAfterTreatment(k[numKmer - 1], allReadFP, kmer31Graph, kmer31Counter, numThread);
	} else if (k[numKmer - 1] <= 64) {
		outputAndAfterTreatment(k[numKmer - 1], allReadFP, kmer63Graph, kmer63Counter, numThread);
	} else if (k[numKmer - 1] <= 96) {
		outputAndAfterTreatment(k[numKmer - 1], allReadFP, kmer95Graph, kmer95Counter, numThread);
	} else if (k[numKmer - 1] <= 128) {
		outputAndAfterTreatment(k[numKmer - 1], allReadFP, kmer127Graph, kmer127Counter, numThread);
	} else if (k[numKmer - 1] <= 160) {
		outputAndAfterTreatment(k[numKmer - 1], allReadFP, kmer159Graph, kmer159Counter, numThread);
	} else {
		outputAndAfterTreatment(k[numKmer - 1], allReadFP, kmerNGraph, kmerNCounter, numThread);
	}

	for (unsigned long long i = 0; i < numThread; ++i)
		fclose(readFP[i]);

    cerr << "assemble completed!" << endl;
}


template <typename KMER>
void Assemble::mergeContig(platanus::Contig &contig, Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength, FILE **mergedContigFP)
{
    graph.setBubbleAndBranch(atof(optionSingleArgs["-u"].c_str()), atof(optionSingleArgs["-d"].c_str()));
    counter.setKmerLength(kmerLength);
    graph.setKmerLength(kmerLength);
	double averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
    unsigned long long doubleHashSize = counter.makeKmerReadDistributionFromContig(contig, kmerLength, 1, memory);

    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(0);
    counter.loadKmer(0, doubleHashSize);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(numThread);
    graph.crushBubbleIterative(averageCoverage);

    graph.saveContigSimple(*mergedContigFP, 1.0);
}




//////////////////////////////////////////////////////////////////////////////////////
// initial assemble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::initialKmerAssemble(vector<unsigned> &k, FILE **readFP, const unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const double coverageCutoffFactor)
{
    cerr << "K = " << k[0] << ", saving kmers from reads..." << endl;
    FILE *sortedKeyFP;
    std::ostringstream oss;

    counter.setKmerLength(k[0]);
    graph.setKmerLength(k[0]);
	coverageCutoff.clear();

    // make kmer distribution
    this->doubleHashSize = counter.makeKmerReadDistributionMT(k[0], readFP, memory, numThread);

    // calculate various value and decide kmer extension
	if (!repeatModeFlag)
		coverageCutoff.push_back(std::max((atoi(optionSingleArgs["-n"].c_str()) != 0 ? atoi(optionSingleArgs["-n"].c_str()) : counter.getLeftLocalMinimalValue(SMOOTHING_WINDOW) / 2), 2ull));
	else
		coverageCutoff.push_back(std::max((atoi(optionSingleArgs["-n"].c_str()) != 0 ? atoi(optionSingleArgs["-n"].c_str()) : counter.getLeftLocalMinimalValue(SMOOTHING_WINDOW)), 2ull));

	averageCoverage = counter.calcOccurrenceDistributionAverage(coverageCutoff[0], counter.getMaxOccurrence());
	if (averageCoverage * coverageCutoffFactor > coverageCutoff.back())
		coverageCutoff.back() = averageCoverage * coverageCutoffFactor;
	averageCoverage = counter.calcOccurrenceDistributionAverage(coverageCutoff[0], counter.getMaxOccurrence());

	if (optionSingleArgs["-e"] != "")
		averageCoverage = atof(optionSingleArgs["-e"].c_str());

    averageLength = counter.calcLengthDistributionAverage(0, platanus::ConstParam::MAX_READ_LEN);
	averageCoverage = averageCoverage * averageLength / (averageLength - k[0] + 1.0);

    cerr << "AVE_READ_LEN=" << averageLength << endl;
    numKmer = this->extendKmer(minLogPJoin, averageCoverage, averageLength, minCoverage, k, lengthStep);

    // output kmer frequency distribution
    oss << optionSingleArgs["-o"] << '_' << k[0] << "merFrq.tsv";
    string kmerFrqFilename = oss.str();
    counter.outputOccurrenceDistribution(kmerFrqFilename);

    // load kmer from temporary file
    sortedKeyFP = counter.sortedKeyFromKmerFile(coverageCutoff[0]);
    this->doubleHashSize = counter.loadKmer(coverageCutoff[0], this->doubleHashSize);

    // make bruijn graph
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(numThread);

/*
    graph.cutBranchIterative(numThread);
    graph.deleteErroneousStraightNodeIterative(2 * k[0], 2 *coverageCutoff[0], numThread);
    graph.cutBranchIterative(numThread);
*/
/*
    unsigned kmerCoverage = averageCoverage * (averageLength - k[0] + 1.0) / averageLength + 0.5;
	vector<unsigned long long> lengthThreshold{std::min(2 * k[0], static_cast<unsigned>(averageLength)), std::max(2 * k[0], static_cast<unsigned>(averageLength))};
	vector<unsigned long long> coverageThreshold{std::min(2 * coverageCutoff[0], static_cast<unsigned long long>(kmerCoverage / 2.0 + 0.5)), std::max(2 * coverageCutoff[0], static_cast<unsigned long long>(kmerCoverage / 2.0 + 0.5))};
	for (unsigned i = 0; i < lengthThreshold.size(); ++i) {
		graph.cutBranchIterative(numThread);
		graph.deleteErroneousStraightNodeIterative(lengthThreshold[i], coverageThreshold[i], numThread);
		graph.cutBranchIterative(numThread);
	}
*/
	if (optionSingleArgs["-e"] == "")
		averageCoverage = graph.getAverageCoverageExcludingBubble();
	else
		averageCoverage = atof(optionSingleArgs["-e"].c_str());

	if (repeatModeFlag) {
		graph.deleteErroneousStraightNodeIterative(ULONG_MAX, REPEAT_MODE_CUTOFF_FACTOR * averageCoverage + 0.5, numThread);
		graph.crushBubbleIterative(DBL_MAX);
	}

    averageCoverage = averageCoverage * averageLength / (averageLength - k[0] + 1.0);
}


//////////////////////////////////////////////////////////////////////////////////////
// save contig and kmer extension assemble
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::saveAndRedoAssemble(const unsigned k, const unsigned previousk, FILE **readFP, FILE **allReadFP, const unsigned long long numThread, BruijnGraph<KMER> &lowKmerGraph, Counter<KMER> &lowKmerCounter, const unsigned long long position)
{

//    lowKmerGraph.divideStraightNode(allReadFP, coverageCutoff[position - 1], this->doubleHashSize, numThread);

	if (position <= 1)
		saveGraph(k, readFP, numThread, lowKmerGraph, lowKmerCounter, false);
	else
		saveGraph(k, readFP, numThread, lowKmerGraph, lowKmerCounter, true);

    if (k <= 32) {
        updateDoubleHashSize<Kmer31::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<Kmer31::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<Kmer31>(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer31Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer31Graph, kmer31Counter, position);
    } else if (k <= 64) {
        updateDoubleHashSize<KmerN<Binstr63>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr63>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr63> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer63Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer63Graph, kmer63Counter, position);
    } else if (k <= 96) {
        updateDoubleHashSize<KmerN<Binstr95>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr95>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr95> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer95Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer95Graph, kmer95Counter, position);
    } else if (k <= 128) {
        updateDoubleHashSize<KmerN<Binstr127>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr127>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr127> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer127Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer127Graph, kmer127Counter, position);
    } else if (k <= 160) {
        updateDoubleHashSize<KmerN<Binstr159>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<Binstr159>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<Binstr159> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmer159Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmer159Graph, kmer159Counter, position);
    } else {
        updateDoubleHashSize<KmerN<binstr_t>::keyType>(lowKmerGraph, k, previousk);
        DoubleHash<KmerN<binstr_t>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        lowKmerGraph.template saveContig<KmerN<binstr_t> >(k, (averageLength - k + 1.0) / (averageLength - previousk + 1.0), tmpOccurrenceTable);
        lowKmerGraph.deleteAllTable();
        kmerNCounter.swapOccurrenceTable(k, tmpOccurrenceTable);
        redoAssemble(k, previousk, readFP, allReadFP, numThread, kmerNGraph, kmerNCounter, position);
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// save contig made brujin graph and prepare next kmer step
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::saveGraph(const unsigned k, FILE **readFP, const unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const bool readFPCloseFlag)
{
    unsigned long long i;
    counter.setOccurrenceTableSize(this->doubleHashSize);
    graph.saveEdgeKmer(counter, k);
    cerr << "extracting reads (containing kmer used in contig assemble)..." << endl;
    # pragma omp parallel for schedule(static, 1)
    for (i = 0; i < numThread; ++i)
        counter.pickupReadMatchedEdgeKmer(&readFP[i], readFPCloseFlag);
    counter.deleteAllTable();
}



//////////////////////////////////////////////////////////////////////////////////////
// load new kmer and make de brujin graph
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::redoAssemble(const unsigned k, const unsigned previousk, FILE **readFP, FILE **allReadFP, unsigned long long numThread, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const unsigned long long position)
{
    counter.makeKmerReadDistributionConsideringPreviousGraph(k, readFP, memory, numThread);
    cerr << "COVERAGE_CUTOFF = " << coverageCutoff[position] << endl;

    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(coverageCutoff[position]);
    this->doubleHashSize = counter.loadKmer(coverageCutoff[position], this->doubleHashSize);
    graph.setKmerLength(k);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();
    graph.cutBranchIterative(numThread);
/*
    graph.cutBranchIterative(numThread);
    graph.deleteErroneousStraightNodeIterative(2 * k, 2 * coverageCutoff[position], numThread);
    graph.cutBranchIterative(numThread);
*/
    unsigned kmerCoverage = averageCoverage * (averageLength - k + 1.0) / averageLength + 0.5;
/*
	vector<unsigned long long> lengthThreshold{std::min(2 * k, static_cast<unsigned>(averageLength)), std::max(2 * k, static_cast<unsigned>(averageLength))};
	vector<unsigned long long> coverageThreshold{std::min(2 * coverageCutoff[position], static_cast<unsigned long long>(kmerCoverage / 2.0 + 0.5)), std::max(2 * coverageCutoff[position], static_cast<unsigned long long>(kmerCoverage / 2.0 + 0.5))};
	for (unsigned i = 0; i < lengthThreshold.size(); ++i) {
		graph.cutBranchIterative(numThread);
		graph.deleteErroneousStraightNodeIterative(lengthThreshold[i], coverageThreshold[i], numThread);
		graph.cutBranchIterative(numThread);
	}
*/
	if (repeatModeFlag) {
		graph.deleteErroneousStraightNodeIterative(ULONG_MAX, REPEAT_MODE_CUTOFF_FACTOR *  kmerCoverage + 0.5, numThread);
		graph.crushBubbleIterative(DBL_MAX);
	}
}

//////////////////////////////////////////////////////////////////////////////////////
// after treatment and output function
// after treatment contains below
// 1. delete erroneous straight node (using coverage and length)
// 2. crush bubble (using coverage)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Assemble::outputAndAfterTreatment(const unsigned k, FILE **readFP, BruijnGraph<KMER> &graph, Counter<KMER> &counter, const unsigned long long numThread)
{
    string outputFilename;

    // graph.cutBranchIterative(numThread);

    unsigned long long leftMinimalCoverage = graph.getLeftMinimalCoverage();
    unsigned long long lengthCutoff = k * 2;
    cerr << "LENGTH_CUTOFF = " << lengthCutoff << endl;
    cerr << "COVERAGE_CUTOFF = " << leftMinimalCoverage << endl;
    graph.deleteErroneousStraightNodeIterative(lengthCutoff, leftMinimalCoverage, numThread);

    averageCoverage = averageCoverage * (averageLength - k + 1.0) / averageLength;
    cerr << "AVE_KMER_COV_REMOVING_BUBBLE=" << averageCoverage << endl;
    outputFilename = optionSingleArgs["-o"];
    outputFilename += "_contigBubble.fa";
    graph.crushBubbleIterative(averageCoverage);
	if (ceil(atoi(optionSingleArgs["-u"].c_str())) > 0)
		graph.printBubble(outputFilename, averageLength / (averageLength - k + 1.0));

//    graph.divideStraightNode(readFP, coverageCutoff[numKmer - 1], this->doubleHashSize, numThread);

    saveGraph(k, readFP, numThread, graph, counter, true);
    if (k <= 32) {
        DoubleHash<Kmer31::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        graph.template saveContig<Kmer31>(k, 1.0, tmpOccurrenceTable);
        graph.deleteAllTable();
        kmer31Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
    } else if (k <= 64) {
        long base = sizeof(std::pair<typename KmerN<Binstr63>::keyType, unsigned short>);
        unsigned long long tmpMemory = memory / base;
        base = log(tmpMemory) / log(2);
        tmpMemory = std::pow(2, base);
        if (tmpMemory > memory) {
            while (tmpMemory > memory) {
                tmpMemory >>= 1;
            }
        }
        this->doubleHashSize = tmpMemory;
        DoubleHash<KmerN<Binstr63>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        graph.template saveContig<KmerN<Binstr63> >(k, 1.0, tmpOccurrenceTable);
        graph.deleteAllTable();
        kmer63Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
    } else if (k <= 96) {
        long base = sizeof(std::pair<typename KmerN<Binstr95>::keyType, unsigned short>);
        unsigned long long tmpMemory = memory / base;
        base = log(tmpMemory) / log(2);
        tmpMemory = std::pow(2, base);
        if (tmpMemory > memory) {
            while (tmpMemory > memory) {
                tmpMemory >>= 1;
            }
        }
        this->doubleHashSize = tmpMemory;
        DoubleHash<KmerN<Binstr95>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        graph.template saveContig<KmerN<Binstr95> >(k, 1.0, tmpOccurrenceTable);
        graph.deleteAllTable();
        kmer95Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
    } else if (k <= 128) {
        long base = sizeof(std::pair<typename KmerN<Binstr127>::keyType, unsigned short>);
        unsigned long long tmpMemory = memory / base;
        base = log(tmpMemory) / log(2);
        tmpMemory = std::pow(2, base);
        if (tmpMemory > memory) {
            while (tmpMemory > memory) {
                tmpMemory >>= 1;
            }
        }
        this->doubleHashSize = tmpMemory;
        DoubleHash<KmerN<Binstr127>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        graph.template saveContig<KmerN<Binstr127> >(k, 1.0, tmpOccurrenceTable);
        graph.deleteAllTable();
        kmer127Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
    } else if (k <= 160) {
        long base = sizeof(std::pair<typename KmerN<Binstr159>::keyType, unsigned short>);
        unsigned long long tmpMemory = memory / base;
        base = log(tmpMemory) / log(2);
        tmpMemory = std::pow(2, base);
        if (tmpMemory > memory) {
            while (tmpMemory > memory) {
                tmpMemory >>= 1;
            }
        }
        this->doubleHashSize = tmpMemory;
        DoubleHash<KmerN<Binstr159>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        graph.template saveContig<KmerN<Binstr159> >(k, 1.0, tmpOccurrenceTable);
        graph.deleteAllTable();
        kmer159Counter.swapOccurrenceTable(k, tmpOccurrenceTable);
    } else {
        long base = sizeof(std::pair<typename KmerN<binstr_t>::keyType, unsigned short>);
        unsigned long long tmpMemory = memory / base;
        base = log(tmpMemory) / log(2);
        tmpMemory = std::pow(2, base);
        if (tmpMemory > memory) {
            while (tmpMemory > memory) {
                tmpMemory >>= 1;
            }
        }
        this->doubleHashSize = tmpMemory;
        DoubleHash<KmerN<binstr_t>::keyType, unsigned short> tmpOccurrenceTable(this->doubleHashSize);
        graph.template saveContig<KmerN<binstr_t> >(k, 1.0, tmpOccurrenceTable);
        graph.deleteAllTable();
        kmerNCounter.swapOccurrenceTable(k, tmpOccurrenceTable);
    }

    counter.makeKmerReadDistributionConsideringPreviousGraph(k, readFP, memory, numThread);
    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(coverageCutoff[numKmer - 1]);
    this->doubleHashSize = counter.loadKmer(coverageCutoff[numKmer - 1], this->doubleHashSize);
    graph.setKmerLength(k);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);

    outputFilename = optionSingleArgs["-o"];
    outputFilename += "_contig.fa";

    FILE *contigFP = platanus::makeTemporaryFile();
    graph.saveContigSimple(contigFP, averageLength / (averageLength - k + 1.0));
    platanus::printContig(outputFilename, contigFP, 1.0, averageLength, k, "seq");
	fclose(contigFP);

    contigFP = platanus::makeTemporaryFile();
    graph.saveJunction(contigFP, averageLength / (averageLength - k + 1.0));
    platanus::printContig(outputFilename, contigFP, 1.0, averageLength, k, "junction", true);
	fclose(contigFP);
}


//////////////////////////////////////////////////////////////////////////////////////
// decide how extend kmer
//////////////////////////////////////////////////////////////////////////////////////
/*
unsigned long long Assemble::extendKmer(const double minLogPJoin, const double averageCoverage, const double averageLength, const unsigned long long minCoverage, vector<unsigned> &kmer, const unsigned long long lengthStep)
{
    unsigned long long i;

    unsigned long long maxK = static_cast<long>(averageLength * maxKmerLengthRatio + 0.5);

    std::cerr << "\nKMER_EXTENSION:" << std::endl;
//    cerr << "K=" << kmer[0] << ", KMER_COVERAGE=" << averageCoverage * (averageLength - kmer[0] + 1.0) / averageLength;
//    cerr << " (>= " << coverageCutoff[0] << "), COVERAGE_CUTOFF=" << coverageCutoff[0] << endl;
    cerr << "K=" << kmer[0] << ", KMER_COVERAGE=" << averageCoverage * (averageLength - kmer[0] + 1.0) / averageLength << endl;

    for (i = 1; kmer[i - 1] <= maxK; ++i) {
        kmer.push_back(0);
        coverageCutoff.push_back(0);
		kmer[i] = std::min(maxK, kmer[i - 1] + lengthStep);
        coverageCutoff[i] = std::max(minCoverage, (this->decreaseCoverageCutoff(coverageCutoff[i - 1], averageCoverage, averageLength, minLogPJoin, kmer[i], kmer[i - 1])) / 2);
        if (kmer[i] == kmer[i - 1]) {
			kmer.pop_back();
            break;
		}

        cerr << "K=" << kmer[i] << ", KMER_COVERAGE=" << averageCoverage * (averageLength - kmer[i] + 1.0) / averageLength;
        cerr << ", COVERAGE_CUTOFF=" << coverageCutoff[i] << ", PROB_SPLIT=10e" << log10(1.0 - exp(calcLogProbabilityJoin(coverageCutoff[i], averageCoverage, averageLength, kmer[i], kmer[i - 1]))) << endl;
    }

    return i;
}
*/
unsigned long long Assemble::extendKmer(const double minLogPJoin, const double averageCoverage, const double averageLength, const unsigned long long minCoverage, vector<unsigned> &kmer, const unsigned long long lengthStep)
{
    unsigned long long i, j;

    unsigned long long minMaxK = static_cast<long>(averageLength * maxKmerLengthRatio + 0.5);

    std::cerr << "\nKMER_EXTENSION:" << std::endl;
    cerr << "K=" << kmer[0] << ", KMER_COVERAGE=" << averageCoverage * (averageLength - kmer[0] + 1.0) / averageLength;
    cerr << " (>= " << coverageCutoff[0] << "), COVERAGE_CUTOFF=" << coverageCutoff[0] << endl;

    for (i = 1; kmer[i - 1] <= averageLength; ++i) {
        kmer.push_back(0);
        coverageCutoff.push_back(0);
        for (j = 1; j <= lengthStep + 1; ++j) {
            kmer[i] = kmer[i - 1] + j;
            coverageCutoff[i] = this->decreaseCoverageCutoff(coverageCutoff[i - 1], averageCoverage, averageLength, minLogPJoin, kmer[i], kmer[i - 1]);
            if (coverageCutoff[i] < minCoverage)
                coverageCutoff[i] = minCoverage;

            if (kmer[i - 1] + j > minMaxK && calcLogProbabilityJoin(coverageCutoff[i], averageCoverage, averageLength, kmer[i], kmer[i - 1]) < minLogPJoin) {
                break;
            }
        }

        --(kmer[i]);
        coverageCutoff[i] = this->decreaseCoverageCutoff(coverageCutoff[i - 1], averageCoverage, averageLength, minLogPJoin, kmer[i], kmer[i - 1]);
		if (coverageCutoff[i] < minCoverage)
			coverageCutoff[i] = minCoverage;

        if (kmer[i] == kmer[i - 1])
            break;

        cerr << "K=" << kmer[i] << ", KMER_COVERAGE=" << averageCoverage * (averageLength - kmer[i] + 1.0) / averageLength;
        cerr << ", COVERAGE_CUTOFF=" << coverageCutoff[i] << ", PROB_SPLIT=10e" << log10(1.0 - exp(calcLogProbabilityJoin(coverageCutoff[i], averageCoverage, averageLength, kmer[i], kmer[i - 1]))) << endl;
    }

    return i;
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate something value (I don't know detail...)
//////////////////////////////////////////////////////////////////////////////////////
double Assemble::calcLogProbabilityJoin(const unsigned long long coverageCutoff, const double averageCoverage, const double averageLength, const unsigned long long largeKmer, const unsigned long long smallKmer) const
{
    double p;
    double s = 0;
    double tmpAverageCoverage = averageCoverage * (averageLength - largeKmer + 1.0) / averageLength;

    for (unsigned long long i = 0; i < coverageCutoff; ++i) {
        p = 0;
        for (unsigned long long j = 1; j <= i; ++j)
            p += log(tmpAverageCoverage) - log(static_cast<double>(j));
        s += exp(p);
    }
    s = exp(-tmpAverageCoverage + log(s));

    return ((largeKmer - smallKmer) + 1.0) * (-s);
}




//////////////////////////////////////////////////////////////////////////////////////
// decrease coverage cutoff value
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long Assemble::decreaseCoverageCutoff(const unsigned long long coverageCutoff, const double averageCoverage, const double averageLength, const double minLogPJoin, const unsigned long long largeKmer, const unsigned long long smallKmer) const
{

    unsigned long long i;

    if (coverageCutoff <= 1)
        return 1;

    for (i = coverageCutoff; i > 1; --i) {
        if (calcLogProbabilityJoin(i, averageCoverage, averageLength, largeKmer, smallKmer) > minLogPJoin)
            break;
    }
    return i;
}



//////////////////////////////////////////////////////////////////////////////////////
// calculate max kmer extension
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long Assemble::calcMaxKmerLength(const double minimumLogPJoin, const double averageCoverage, const double averageLength, const unsigned long long minimumCoverage, unsigned long long kmerLength, const unsigned long long lengthStep) const
{
    unsigned long long i;

    for (kmerLength += lengthStep; kmerLength < (unsigned long long)(averageLength + 0.5); kmerLength += lengthStep) {
        if (calcLogProbabilityJoin(minimumCoverage, averageCoverage, averageLength, kmerLength, kmerLength - lengthStep) < minLogPJoin)
            break;
    }
    kmerLength -= lengthStep;

    for (i = 1; i < lengthStep; ++i) {
        if (calcLogProbabilityJoin(minimumCoverage, averageCoverage, averageLength, kmerLength + i, kmerLength) < minLogPJoin)
            break;
    }
    --i;

    return kmerLength + i;
}

template <typename KEY, typename KMER>
void Assemble::updateDoubleHashSize(const BruijnGraph<KMER> &graph, const unsigned k, const unsigned previousk)
{
    unsigned long long estimate = graph.estimateNumKmerOnStraight();
    unsigned long long size = log(estimate / platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR) / log(2);
    size = std::pow(2, size + 1);

    long base = sizeof(std::pair<KEY, unsigned short>);
    unsigned long long tmpMemory = this->memory / base;
    base = log(tmpMemory) / log(2);
    this->doubleHashSize = std::pow(2, base);
    if (this->doubleHashSize > this->memory) {
        this->doubleHashSize >>= 1;
    }

    if (size > this->doubleHashSize) {
        platanus::MemoryAlert();
    }
    this->doubleHashSize = std::max(size, this->doubleHashSize);

}


// below this line there are only read function

//////////////////////////////////////////////////////////////////////////////////////
// read functions
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::readInputFile(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
	platanus::FILETYPE fileType = checkFileFormat(inputFilename);
	if (fileType == platanus::FILETYPE::FASTA)
		readFasta(inputFilename, outputMT, numThread);
	else if (fileType > platanus::FILETYPE::FASTA)
		readFastq(inputFilename, outputMT, numThread);
	else
		throw platanus::ReadError("Read file exception!!\nRead file is not FASTA/FASTQ format.");
}


//////////////////////////////////////////////////////////////////////////////////////
// read FASTA
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::readFasta(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
	platanus::FILECOMPRESSION format = platanus::checkFileCompression(inputFilename);

	if (format == platanus::FILECOMPRESSION::UNCOMPRESSED)
		readFastaUncompressed(inputFilename, outputMT, numThread);
	else
		readFastaCompressed(inputFilename, outputMT, numThread);
}


void Assemble::readFastaUncompressed(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    ifstream ifs(inputFilename);
    SEQ seq;
    string oneLine;
    string read = "";
    unsigned i = 0;

    for (unsigned long long i = 0; i < numThread; ++i)
        fseek(outputMT[i], 0, SEEK_END);
    while (ifs && getline(ifs, oneLine)) {
        if (oneLine[0] == '>') {
            break;
        }
    }

    while (ifs && getline(ifs, oneLine)) {
        if (oneLine[0] != '>') {
            read += oneLine;
        } else {
            if (read.length() != 0) {
                seq.convertFromString(read);
                seq.writeTemporaryFile(outputMT[i]);
                read = "";
                i = (i + 1) % numThread;
            }
        }
    }
    seq.convertFromString(read);
    seq.writeTemporaryFile(outputMT[i]);

    ifs.close();
}


void Assemble::readFastaCompressed(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    SEQ seq;
    string oneLine;
    string read = "";
    unsigned i = 0;

	FILE* fp = platanus::openFileAllowingCompression(inputFilename, "r");

    for (unsigned long long i = 0; i < numThread; ++i)
        fseek(outputMT[i], 0, SEEK_END);

    while (platanus::getlineFILE(oneLine, fp) != NULL) {
        if (oneLine[0] == '>') {
            break;
        }
    }

    while (platanus::getlineFILE(oneLine, fp) != NULL) {
        if (oneLine[0] != '>') {
            read += oneLine;
        } else {
            if (read.length() != 0) {
                seq.convertFromString(read);
                seq.writeTemporaryFile(outputMT[i]);
                read = "";
                i = (i + 1) % numThread;
            }
        }
    }
    seq.convertFromString(read);
    seq.writeTemporaryFile(outputMT[i]);

	platanus::closeFileAllowingCompression(fp, inputFilename);
}


//////////////////////////////////////////////////////////////////////////////////////
// read FASTQ
//////////////////////////////////////////////////////////////////////////////////////
void Assemble::readFastq(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
	platanus::FILECOMPRESSION format = platanus::checkFileCompression(inputFilename);

	if (format == platanus::FILECOMPRESSION::UNCOMPRESSED)
		readFastqUncompressed(inputFilename, outputMT, numThread);
	else
		readFastqCompressed(inputFilename, outputMT, numThread);
}


void Assemble::readFastqUncompressed(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    ifstream ifs(inputFilename);
    SEQ seq;
    string oneLine;
    string read = "";
    unsigned i = 0;
    bool flag = true;

    for (unsigned long long i = 0; i < numThread; ++i)
        fseek(outputMT[i], 0, SEEK_END);
    while (ifs && getline(ifs, oneLine)) {
        if (oneLine[0] == '@') {
            break;
        }
    }

    while (ifs && getline(ifs, oneLine)) {
        if (oneLine.length() != 0) {
            if (oneLine[0] != '@') {
                if (flag && oneLine[0] != '+')
                    read += oneLine;
                else
                    flag = false;
            } else {
                if (read.length() != 0) {
                    seq.convertFromString(read);
                    seq.writeTemporaryFile(outputMT[i]);
                    read = "";
                    i = (i + 1) % numThread;
                }
                flag = true;
            }
        }
    }
    seq.convertFromString(read);
    seq.writeTemporaryFile(outputMT[i]);

    ifs.close();

}


void Assemble::readFastqCompressed(const string &inputFilename, FILE **outputMT, const unsigned long long numThread) const
{
    SEQ seq;
    string oneLine;
    string read = "";
    unsigned i = 0;
    bool flag = true;

	FILE* fp = platanus::openFileAllowingCompression(inputFilename, "r");

    for (unsigned long long i = 0; i < numThread; ++i)
        fseek(outputMT[i], 0, SEEK_END);

    while (platanus::getlineFILE(oneLine, fp) != NULL) {
        if (oneLine[0] == '@') {
            break;
        }
    }

    while (platanus::getlineFILE(oneLine, fp) != NULL) {
        if (oneLine.length() != 0) {
            if (oneLine[0] != '@') {
                if (flag && oneLine[0] != '+')
                    read += oneLine;
                else
                    flag = false;
            } else {
                if (read.length() != 0) {
                    seq.convertFromString(read);
                    seq.writeTemporaryFile(outputMT[i]);
                    read = "";
                    i = (i + 1) % numThread;
                }
                flag = true;
            }
        }
    }
    seq.convertFromString(read);
    seq.writeTemporaryFile(outputMT[i]);

	platanus::closeFileAllowingCompression(fp, inputFilename);
}
