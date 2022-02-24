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

#include "merge.h"




//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
ContigMerger::ContigMerger()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-c"] = "1";
    optionSingleArgs["-k"] = "1.0";
    optionSingleArgs["-l"] = "2.0";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-d"] = "0.5";
    optionSingleArgs["-m"] = "16";
    optionMultiArgs["-f"] = std::vector<std::string>();
    optionSingleArgs["-tmp"] = ".";
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void ContigMerger::usage(void) const
{
    std::cerr << "\nUsage: platanus_allee merge [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -f FILE1 [FILE2 ...]               : contig or scaffold_file (fasta format)\n"
              << "    -c INT                             : minimum coverage (default " << optionSingleArgs.at("-c") << ")\n"
              << "    -k FLOAT                           : k-mer size factor (k = FLOAT * read_length) (default read from fasta-header)\n"
              << "    -l FLOAT                           : minimum length factor (minimum_length = FLOAT * read_length) (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -d FLOAT                           : maximum difference for branch cutting (coverage ratio, default " << optionSingleArgs.at("-d") << ")\n"
              << "    -m INT                             : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Output:\n"
              << "    PREFIX_merged.fa\n"
              << std::endl;
}




//////////////////////////////////////////////////////////////////////////////////////
// exec divide
//////////////////////////////////////////////////////////////////////////////////////
void ContigMerger::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
    this->readLength = this->contig.getReadLengthFromFastaHeader(optionMultiArgs["-f"][0]);
    unsigned long long kmerLength = atof(optionSingleArgs["-k"].c_str()) * readLength + 0.5;

    if (kmerLength <= 32) {
        Counter<Kmer31> counter;
        BruijnGraph<Kmer31> graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 64) {
        Counter<KmerN<Binstr63> > counter;
        BruijnGraph<KmerN<Binstr63> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 96) {
        Counter<KmerN<Binstr95> > counter;
        BruijnGraph<KmerN<Binstr95> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 128) {
        Counter<KmerN<Binstr127> > counter;
        BruijnGraph<KmerN<Binstr127> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else if (kmerLength <= 160) {
        Counter<KmerN<Binstr159> > counter;
        BruijnGraph<KmerN<Binstr159> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    } else {
        Counter<KmerN<binstr_t> > counter;
        BruijnGraph<KmerN<binstr_t> > graph;
        this->exec2_ForIntegrateKmer(counter, graph, kmerLength);
    }
}



template <typename KMER>
void ContigMerger::exec2_ForIntegrateKmer(Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength)
{
    unsigned long long memory = atoi(optionSingleArgs["-m"].c_str()) * static_cast<unsigned long long>(1000000000);
    unsigned long long contigReadLength = platanus::Contig::getReadLengthFromFastaHeader(optionMultiArgs["-f"][0]);
    unsigned long long minContigLength = atof(optionSingleArgs["-l"].c_str()) * contigReadLength + 0.5;
    graph.setBubbleAndBranch(atof(optionSingleArgs["-u"].c_str()), atof(optionSingleArgs["-d"].c_str()));
    counter.setKmerLength(kmerLength);
    graph.setKmerLength(kmerLength);
    for (auto itr = optionMultiArgs["-f"].begin(), end = optionMultiArgs["-f"].end(); itr != end; ++itr) {
        contig.readFastaCoverageCutN(*itr, minContigLength);
        //contig.readFastaCoverage(*itr);
    }
	double averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
    unsigned long long doubleHashSize = counter.makeKmerReadDistributionFromContig(contig, minContigLength, atoi(optionSingleArgs["-c"].c_str()), memory);
    contig.clear();

    FILE *sortedKeyFP = counter.sortedKeyFromKmerFile(0);
    counter.loadKmer(0, doubleHashSize);
    graph.makeInitialBruijnGraph(counter, sortedKeyFP);
    fclose(sortedKeyFP);
    counter.deleteAllTable();

    graph.cutBranchIterative(1);
    graph.crushBubbleIterative(averageCoverage);
    std::string outputFilename = optionSingleArgs["-o"];
    outputFilename += "_merged.fa";
    graph.printContig(outputFilename, 1.0, this->readLength);

    FILE *contigFP = platanus::makeTemporaryFile();
    graph.saveJunction(contigFP, 1.0);
    outputFilename = optionSingleArgs["-o"];
    outputFilename += "_mergedJunctionKmer.fa";
    platanus::printContig(outputFilename, contigFP, 1.0, this->readLength, kmerLength, "junction");
	fclose(contigFP);
}
