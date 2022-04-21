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

#include "solveDBG.h"
#include "seqlib.h"
#include "kmer.h"
#include <cmath>
#include <algorithm>
#include <sstream>
#include <climits>
#include <cfloat>
#include <iomanip>

#include <mcheck.h>

using std::vector;
using std::string;
using std::unordered_map;
using std::cerr;
using std::endl;


//////////////////////////////////////////////////////////////////////////////////////
// const parameter define
//////////////////////////////////////////////////////////////////////////////////////
const double SolveDBG::MIN_LONG_READ_LENGTH_CUTOFF_FACTOR = 1;
const double SolveDBG::MAX_LONG_READ_LENGTH_CUTOFF_FACTOR = 4;
const double SolveDBG::LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR = 0.1;
const long SolveDBG::LONG_READ_MIN_ALIGNMENT_LENGTH = 1000;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_COVERAGE = 0.8;
const double SolveDBG::LONG_READ_MIN_ALIGNMENT_IDENTITY = 0.8;

//////////////////////////////////////////////////////////////////////////////////////
// default constructor
//////////////////////////////////////////////////////////////////////////////////////
SolveDBG::SolveDBG()
: Scaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-k"] = "1";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-unlink_list"] = "";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();
    optionMultiArgs["-hic"] = vector<string>(); //added by ouchi
    optionMultiArgs["-HIC"] = vector<string>(); //added by ouchi

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    optionMultiArgs["-S"] = vector<string>(1);
    optionMultiArgs["-S"][0] = "20";

    optionBool["-no_scaffold"] = false;
    optionBool["-unphase"] = false;
    optionBool["-reduce_redundancy"] = false;
	optionBool["-divide_only"] = false;
	optionBool["-simplify_junction"] = false;

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";

    optionSingleArgs["-mapper"] = "";
    optionBool["-minimap2_sensitive"] = false;
    optionBool["-minimap2"] = true;
    optionBool["-minialign"] = false;
    optionBool["-kmer_align"] = false;

    optionBool["-fastg"] = false;
    optionBool["-aggressive"] = false;
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::usage(void) const
{
    std::cerr << "\nUsage: platanus_allee solveDBG [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (reads in 1 file, fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : tagged_pair_files (10x Genomics) (reads in 1 file, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : tagged_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)\n"
              << "    -hic PAIR1 [PAIR2 ...]             : HiC_pair_files (reads in 1 file, fasta or fastq)\n" //added by ouchi
              << "    -HIC FWD1 REV1 [FWD2 REV2 ...]     : HiC_pair_files (reads in 2 files, fasta or fastq)\n" //added by ouchi
              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
              << "    -a{INT} INT                        : lib_id average_insert_size\n"
              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -e FLOAT                           : coverage depth of homozygous region (default auto)\n"
              << "    -L INT                             : maximum fragment length of tag (10x Genomics) (default " << optionSingleArgs.at("-L") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1] << " " << optionMultiArgs.at("-s")[2]  << ")\n"
//              << "    -S INT1 [INT2 ...]                 : mapping seed length for long reads (default " << optionMultiArgs.at("-S")[0] << ", only effective with -kmer_align option)\n"
              << "    -k INT                             : minimum number of links to phase variants (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -unphase                           : not phase heterozygous regions and construct consensus scaffolds (default false)\n"
              << "    -simplify_junction                 : only simplify junction-contigs based on DBG-overlaps (default false)\n"
              << "    -mapper FILE                       : path of mapper executable file (default, minimap2; only effective with -p option)\n"
              << "    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)\n"
              << "    -reduce_redundancy                 : reduce redundant sequences that exactly matche others (default, off)\n"
              << "    -divide_only                       : only divide input sequences (default, off)\n"
              << "    -unlink_list                       : reduce redundant sequences that exactly matche others (default, off)\n"
//              << "    -minimap2                          : use minimap2 insterd of minimap (default) for alignment of long reads\n"
//              << "    -minialign                         : use minialign insterd of minimap (default) for alignment of long reads\n"
//              << "    -kmer_align                        : use built-in alignment based on k-mer match insterd of minimap (default) for long reads\n"
//              << "    -fastg                             : output only fastg files of graphs for Bandage (default off)\n"
              << "    -aggressive                        : aggressively extend scaffolds (default off)\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x, -X, -hic and -HIC options.\n" //edited by ouchi
              << "\n\n"


              << "Outputs:\n"
			  << "    PREFIX_*.fa\n"
              << "\n"
              << "Outputs (-fastg):\n"
			  << "    PREFIX_lib#_graph.fastg\n"
             << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// initialize parameters
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::initializeParameters(void)
{
	seedLength = platanus::ConstParam::SCAFFOLD_HASH_OVERLAP;
    multiSeedLengthForShortRead.clear();
	if (!(optionPairFile.empty())) {
		for (auto itr = optionMultiArgs["-s"].begin(); itr != optionMultiArgs["-s"].end(); ++itr) {
			unsigned length = std::stoi(*itr);
			multiSeedLengthForShortRead.push_back(length);
			if (seedLength > length)
				seedLength = length;
		}
	}
	if (optionBool["-kmer_align"] && !(optionMultiArgs["-p"].empty())) {
		for (auto itr = optionMultiArgs["-S"].begin(); itr != optionMultiArgs["-S"].end(); ++itr) {
			unsigned length = std::stoi(*itr);
			multiSeedLengthForLongRead.push_back(length);
			if (seedLength > length)
				seedLength = length;
		}
	}

    keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);
    bubbleThreshold = atof(optionSingleArgs["-u"].c_str());
    minLink = atoi(optionSingleArgs["-l"].c_str());
    minLinkToPhase = atoi(optionSingleArgs["-k"].c_str());
    minOverlapForScaffolding = atoi(optionSingleArgs["-v"].c_str());
    numThread = atoi(optionSingleArgs["-t"].c_str());
    pairedDBG.setSeedLength(seedLength);
    pairedDBG.setMinTolerenceFactor(MIN_TOL_FACTOR);
    pairedDBG.setMaxFragmentLengthOfTag(atoi(optionSingleArgs["-L"].c_str()));

    sort(optionPairFile.begin(), optionPairFile.end());
    numFilePerLibraryID.resize(optionPairFile.size());
    libraryIDList.resize(optionPairFile.size());
    unsigned numLibrary = 0;
    for (unsigned i = 0; i < optionPairFile.size(); ++i) {
        ++(numFilePerLibraryID[numLibrary]);
        libraryIDList[numLibrary] = optionPairFile[i].libraryID;
        if (i + 1 < optionPairFile.size()) { //added by ouchi
            if (optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
                ++numLibrary;
            }
        } else { //added by ouchi
            ++numLibrary; //added by ouchi
        } //added by ouchi
    }

    libraryMT.resize(numLibrary);
    omp_set_num_threads(numThread);

	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
}


//////////////////////////////////////////////////////////////////////////////////////
// exec scaffold
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::exec(void)
{
    initializeParameters();
    mapLibraryAndInitGraph(numThread);

	if (optionBool["-fastg"]) {
		generateGraphFastg();
		cerr << "scaffold completed!" << endl;
		return;
	}

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.makeGraph(numThread);

	if (optionBool["-simplify_junction"]) {
		pairedDBG.joinUnambiguousNodePairIterative(numThread);

		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "seq");
		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_simplified.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_simplifiedComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	if (optionBool["-unphase"]) {
		if (optionSingleArgs["-e"] == "")
			pairedDBG.calculateHeteroAndAverageCoverageUnphase();
		else
			pairedDBG.setHeteroCoverage(atof(optionSingleArgs["-e"].c_str()));

		pairedDBG.clearEdges();

		extendConsensus(4, !(optionBool["-aggressive"]));

		if (!(libraryMT.empty()))
			pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
		else
			pairedDBG.setTolerence(this->contigMaxK);
		pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");

		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_consensusScaffold.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_consensusScaffoldComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}


	if (optionBool["-divide_only"]) {
		if (optionSingleArgs["-e"] == "")
			pairedDBG.calculateHeteroAndAverageCoverageUnphase();
		else
			pairedDBG.setHeteroCoverage(atof(optionSingleArgs["-e"].c_str()));

		pairedDBG.clearEdges();

		extendConsensusToEstimateInsertSize();
		pairedDBG.resetGraph();

		for (unsigned libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex)
			pairedDBG.divideErroneousNodeBaseLevel(libraryIndex, libraryIndex + 1, this->contigReadLength, false, true, true, numThread);

		pairedDBG.loadDividedContigResultSeq(this->contigMaxK, this->contigReadLength);
		if (optionBool["-reduce_redundancy"])
			pairedDBG.markRedundantResultSeq(numThread);

		outputGraph("_divided.fa");
		pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_dividedComponent.bed");
		cerr << "solve_DBG completed!" << endl;
		return;
	}

//	pairedDBG.calcPercentilValue(95, numThread); //added by ouchi

	pairedDBG.extractDBGBubbleInformation();
	pairedDBG.clearEdges();

	extendConsensusToEstimateInsertSize();
	pairedDBG.resetGraph();


	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setCutoffLength(0);
	pairedDBG.makeGraph(numThread);


	pairedDBG.extractDBGBubbleInformation();
//	pairedDBG.setOppositeForkContigIDOverlapped(numThread);
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.setOppositeBubbleContigIDByEndMatch();
//	pairedDBG.setOppositeBubbleContigIDByOneEndMatch();

//	pairedDBG.setForkJunctionContigIDOverlapped();
	pairedDBG.setBubbleJunctionContigIDOverlapped();

	pairedDBG.clearEdges();


	//////added by ouchi for test

		//pairedDBG.calcIdealContact(numThread);
		//pairedDBG.divideNodeBasedOnHiC(numThread);
//		exit(0);

		//pairedDBG.initConsensusNode(numThread);
		//std::cerr << "mergeHeteroNode" << std::endl;
		//pairedDBG.mergeHeteroNode(numThread);


/*		pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
		pairedDBG.setMinLink(1);
		pairedDBG.makeConsensusGraph(numThread);
		std::ostringstream outStream5;
		outStream5 << optionSingleArgs["-o"] << "outIteration_consensus.fastg";
		pairedDBG.outputConsensusFastg(outStream5.str());
		exit(0);*/

/*
		pairedDBG.setMinOverlap(minOverlapForScaffolding);

		pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
		pairedDBG.unsetMode(PairedDBG::BUBBLE_AWARE_MODE);
		for (long libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
			pairedDBG.setTargetLibraryIndex(libraryIndex);
			for (long iteration = 0; iteration < 2; ++iteration) {
				if (iteration == 0)
					pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = "<< libraryMT[libraryIndex][0].getSDInsSize() << endl;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

					pairedDBG.joinUnambiguousConsensusPairIterative(numThread);
					pairedDBG.setMinLink(1);
					pairedDBG.makeConsensusGraph(numThread, false);
					if (iteration == 0)
						pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
					else
						pairedDBG.setMinLink(minLink);
					pairedDBG.deleteErroneousConsensusEdgebyHiC(numThread);
					pairedDBG.detectConflictConsensus();
					pairedDBG.consensusScaffolding();
				}
			}
		}

		if (libraryMT.size() > 0) {
			for (long iteration = 0; iteration < 2; ++iteration) {
				pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
				long linkThreashold = (iteration + 1) * minLink;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[libraryMT.size() - 1][0].getSDInsSize(), 0.1 * libraryMT[libraryMT.size() - 1][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.joinUnambiguousConsensusPairIterative(numThread);
					pairedDBG.setMinLink(1);
					pairedDBG.makeConsensusGraph(numThread);
					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.deleteErroneousConsensusEdgebyHiC(numThread);
					pairedDBG.detectConflictConsensus();
					pairedDBG.consensusScaffolding();
				}
			}
		}


		pairedDBG.phasing(numThread);
		exit(0);


		pairedDBG.setMinLink(1);
		pairedDBG.setCutoffLength(0);
		pairedDBG.makeConsensusGraph(numThread);

//		pairedDBG.makeConsensusPathEdge(numThread);

		//pairedDBG.deleteErroneousConsensusEdgebyHiC(numThread);
		std::ostringstream outStream4;
		outStream4 << optionSingleArgs["-o"] << "outIteration_consensus_scaffolded_corrected.fastg";
		pairedDBG.outputConsensusFastg(outStream4.str());

        pairedDBG.setMinLink(minLink);
        pairedDBG.setCutoffLength(20000);
        pairedDBG.setTolerence(10000);
        pairedDBG.HiC_Scaffolding(numThread);
		exit(0); */
	/////



	for (long outerIteration = 0; outerIteration < 4; ++outerIteration) {
		for (long iteration = 0; iteration < 2; ++iteration) {
			if (optionBool["-no_scaffold"])
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::BUBBLE_AWARE_MODE);
			else if (iteration == 0)
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			pairedDBG.setCutoffLength(0);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setTargetLibraryIndex(i);
				unsigned tolerenceFactor = MAX_TOL_FACTOR;
				pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
				cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
			}

			if (longReadLibraryMT.size() > 0) {
				if (optionBool["-no_scaffold"])
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::BUBBLE_AWARE_MODE);
				else if (iteration == 0)
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				else
					pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize() << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			if (optionBool["-no_scaffold"])
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::BUBBLE_AWARE_MODE);
			else if (iteration == 0 || optionBool["-no_scaffold"])
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE);
			else
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::NON_DBG_OVERLAP_MODE);
			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
			}

			if (optionBool["-no_scaffold"]) {
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
				pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			}
		}

		if (optionBool["-no_scaffold"]) {
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
			pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
			continue;
		}


		pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
		pairedDBG.setCutoffLength(0);
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);
		pairedDBG.clearEdges();


		for (long iteration = 0; iteration < 2; ++iteration) {
			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();
					if (iteration > 0)
						pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
//						pairedDBG.joinUnambiguousUniqueNodePairGappedIterative(numThread); //added by ouchi for test

					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

				}
			}

			if (longReadLibraryMT.size() > 0) {
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
				cerr << "[LONG_READ_LIBRARY]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageLength() << endl;
				pairedDBG.setTolerence(2 * this->contigMaxK);
				if (iteration == 0)
					pairedDBG.setMinLink(minLinkToPhase);
				else
					pairedDBG.setMinLink(minLink);

				for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
					pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR * longReadLibraryMT[0].getAverageInsSize());
					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::SCORE, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
				}
			}

			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE);

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MIN_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					if (iteration == 0)
						pairedDBG.setMinLink(minLinkToPhase);
					else {
						pairedDBG.setMinLink(minLink);
						pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					}

					pairedDBG.trimSparseEnd();
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

				}
			}
		}

		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleContigInNonHeteroNode();
		pairedDBG.divideBubbleJunctionNode(false);


		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = 0; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				long linkThreshold;
				if (iteration % 2 == 0)
					linkThreshold = minLink;
				else
					linkThreshold = std::max(minLink, pairedDBG.estimateLink());

				cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize() << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
					cerr << "LENGTH_CUTOFF = " << pairedDBG.getCutoffLength() << endl;

					pairedDBG.trimSparseEnd();

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
//					pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					//pairedDBG.detectRepeatSeverely(pairedDBG.getHeteroCoverage()); //edited by ouchi
					pairedDBG.makeScaffold();

					pairedDBG.joinUnambiguousNodePairGappedIterative(numThread);
					pairedDBG.setMinLink(minLink);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
					//pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::HIC, numThread); //added by ouchi

					pairedDBG.setMinLink(linkThreshold);
					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
//					pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					//pairedDBG.detectRepeatSeverely(pairedDBG.getHeteroCoverage()); //edited by ouchi
					pairedDBG.makeScaffold();

//					if (outerIteration < 2)
//						pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
				}
			}
			if (outerIteration < 2)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeUsingBubbleContigPair(numThread);
		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(false, numThread);


		pairedDBG.setMinOverlap(minOverlapForScaffolding);
		for (long iteration = 0; iteration < 2; ++iteration) {
			long linkThreashold = (iteration + 1) * minLink;
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
				pairedDBG.setTargetLibraryIndex(i);
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTolerenceFactor(tolerenceFactor);
					pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);

					pairedDBG.trimSparseEnd();

					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
//					pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					//pairedDBG.detectRepeatSeverely(pairedDBG.getHeteroCoverage()); //edited by ouchi
					pairedDBG.makeScaffold();


					pairedDBG.setMinLink(minLink);
					pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
					pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
					//pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::HIC, numThread); //added by couhi


					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.makeGraphAllLibraries(numThread);
					pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
//					pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
					pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
					pairedDBG.deleteConflictingBubbleEdge(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
					//pairedDBG.detectRepeatSeverely(pairedDBG.getHeteroCoverage()); //edited by ouchi
					pairedDBG.makeScaffold();


//					if (outerIteration < 2)
//						pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
				}
			}
			if (outerIteration < 2)
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
			else
				pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		}
		pairedDBG.setMinOverlap(this->contigMaxK - 1);


		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMinOverlap(minOverlapForScaffolding);

			pairedDBG.divideNodeUsingBubbleContigPair(numThread);

			for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
				if (outerIteration == 0)
					pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
				else
					pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				pairedDBG.setTolerence(2 * this->contigMaxK);
				pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());

				pairedDBG.trimSparseEnd();

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();

				pairedDBG.setMinLink(minLink);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
				pairedDBG.solveSimpleGappedCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
//				pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
				pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();

				pairedDBG.divideNodeUsingBubbleContigPair(numThread);
				pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
				pairedDBG.setMinLink(minLinkToPhase);
				pairedDBG.setCutoffLength(0);
				pairedDBG.joinUnambiguousNodePairIterative(numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
				pairedDBG.solveSimpleCrossStructureIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);

//				if (outerIteration < 2)
//					pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
			}

			pairedDBG.setMinOverlap(this->contigMaxK - 1);
		}

		/*//added by ouchi
		if (HiCLibraryMT.size() > 0) {
			for (long iteration = 0; iteration < 2; ++iteration) {
				long linkThreashold = (iteration + 1) * minLink;
				pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::BUBBLE_AWARE_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				for (unsigned i = libraryMT.size() - 1; i < libraryMT.size(); ++i) {
					pairedDBG.setTargetLibraryIndex(i);
					for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
						pairedDBG.setTolerenceFactor(tolerenceFactor);
						pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);

						pairedDBG.trimSparseEnd();

						pairedDBG.setMinLink(linkThreashold);
						pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
						pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
						pairedDBG.makeGraphAllLibraries(numThread);
						pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
						pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
						pairedDBG.deleteRepeatEdge();
						pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
						pairedDBG.makeScaffold();

						pairedDBG.setMinLink(minLink);
						pairedDBG.joinUnambiguousNodePairGappedIterativeAllLibraries(numThread);
						pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
						pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::LINK, numThread);
						pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::TAG, numThread);
						pairedDBG.solveSimpleGappedCrossStructureAllLibrariesIterative(PairedDBG::CROSS_RESOLUTION_MODE::HIC, numThread); //added by ouchi

						pairedDBG.setMinLink(linkThreashold);
						pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[i][0].getSDInsSize(), 0.1 * libraryMT[i][0].getAverageInsSize()));
						pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
						pairedDBG.makeGraphAllLibraries(numThread);
						pairedDBG.deleteLongEdge(libraryMT[i][0].getAverageInsSize());
						pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
						pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
						pairedDBG.deleteDifferentBubbleEdgeIterative(numThread);
						pairedDBG.deleteConflictingBubbleEdge(numThread);
						pairedDBG.deleteRepeatEdge();
						pairedDBG.detectRepeat(pairedDBG.getHeteroCoverage());
						pairedDBG.makeScaffold();

//						if (outerIteration < 2)
//							pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
					}
				}
				if (outerIteration < 2)
					pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::GAP, numThread, this->contigMaxK);
				else
					pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
			}
		}
		// */


		pairedDBG.trimSparseEnd();
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
		pairedDBG.divideBubbleJunctionNode(true);


/*
        //added by ouchi for test
		pairedDBG.initConsensusNode(numThread);
		std::cerr << "mergeHeteroNode" << std::endl;
		pairedDBG.mergeHeteroNode(numThread);
		pairedDBG.makeConsensusGraph(numThread);

		std::ostringstream outStream2;
		outStream2 << optionSingleArgs["-o"] << "outIteration" << outerIteration << "_consensus.fastg";
		pairedDBG.outputConsensusFastg(outStream2.str());

		pairedDBG.setMinOverlap(minOverlapForScaffolding);

		pairedDBG.unsetMode(PairedDBG::BUBBLE_AWARE_MODE);
		for (long libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
			pairedDBG.setTargetLibraryIndex(libraryIndex);
			for (long iteration = 0; iteration < 2; ++iteration) {
				if (iteration == 0)
					pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = "<< libraryMT[libraryIndex][0].getSDInsSize() << endl;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

					pairedDBG.joinUnambiguousConsensusPairIterative(numThread);
					pairedDBG.makeConsensusGraph(numThread, false);
					pairedDBG.deleteErroneousConsensusEdgebyHiC(numThread);
					pairedDBG.detectConflictConsensus();
					pairedDBG.consensusScaffolding();
				}
			}
		}

		if (libraryMT.size() > 0) {
			for (long iteration = 0; iteration < 2; ++iteration) {
				pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
				long linkThreashold = (iteration + 1) * minLink;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
					pairedDBG.setMinLink(linkThreashold);
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * std::min(1.0 * libraryMT[libraryMT.size() - 1][0].getSDInsSize(), 0.1 * libraryMT[libraryMT.size() - 1][0].getAverageInsSize()));
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					pairedDBG.joinUnambiguousConsensusPairIterative(numThread);
					pairedDBG.makeConsensusGraph(numThread);
					pairedDBG.deleteErroneousConsensusEdgebyHiC(numThread);
					pairedDBG.detectConflictConsensus();
					pairedDBG.consensusScaffolding();
				}
			}
		}


		pairedDBG.setMinLink(1);
		//pairedDBG.setCutoffLength(0);
		pairedDBG.makeConsensusGraph(numThread);

		pairedDBG.setMinLink(minLink);
		std::ostringstream outStream3;
		outStream3 << optionSingleArgs["-o"] << "outIteration" << outerIteration << "_consensus_scaffolded.fastg";
		pairedDBG.outputConsensusFastg(outStream3.str());

//		pairedDBG.makeConsensusPathEdge(numThread);

		//pairedDBG.deleteErroneousConsensusEdgebyHiC(numThread);
		//std::ostringstream outStream4;
		//outStream4 << optionSingleArgs["-o"] << "outIteration" << outerIteration << "_consensus_scaffolded_corrected.fastg";
		//pairedDBG.outputConsensusFastg(outStream4.str());

        pairedDBG.setMinLink(minLink);
        pairedDBG.setCutoffLength(20000);
        //pairedDBG.HiC_Scaffolding(numThread);
*/
		//



		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::SWITCH, numThread);
		pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);


		if (outerIteration < 2) {
			pairedDBG.joinUnambiguousNodePairIterative(numThread);
			pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
		}
	}


	if (!optionBool["-aggressive"])
		pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, true, true, false, numThread);

	pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);


	pairedDBG.divideBubbleContigInNonHeteroNode();

	pairedDBG.setMode(0);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	pairedDBG.copyAllNodes(phasedGraph);

	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();


	extendConsensus(1, false);

	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");
	outputGraph("_preliminaryConsensusScaffold.fa");
	pairedDBG.clearResultSeq();

	pairedDBG.setMode(0);
	pairedDBG.remakeGraphRecoveringSecondaryBubble(phasedGraph);
	pairedDBG.makeGraph(numThread);
	pairedDBG.divideNodeBasedOnBubblesIterative(true, numThread);
	pairedDBG.makeGraph(numThread);
	pairedDBG.adjustOppositeBubbleNodeIDDirection();
	if (!(libraryMT.empty()))
		pairedDBG.setTolerence(MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
	else
		pairedDBG.setTolerence(this->contigMaxK);

    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");

	if (optionBool["-reduce_redundancy"])
		pairedDBG.markRedundantResultSeq(numThread);

	pairedDBG.outputResultSeqWithBubble(optionSingleArgs["-o"], "_primaryBubble.fa", "_secondaryBubble.fa", "_primaryFork.fa", "_secondaryFork.fa", "_nestedBubble.fa", "_nonBubbleOther.fa", "_bubbleRelation.tsv", this->contigMaxK, this->contigReadLength);
    pairedDBG.loadResultSeq(this->contigMaxK, this->contigReadLength, "scaffold");
	pairedDBG.outputResultSeqComponent(optionSingleArgs["-o"], "_phasedScaffoldComponent.bed");
	cerr << "solve_DBG completed!" << endl;
	return;
}


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::mapLibraryAndInitGraph(const int numThread)
{
    platanus::Contig contig;

    std::unique_ptr<HeteroMapper> mapper(new HeteroMapper(seedLength, keyLength));

    readLibrary(mapper, contig, numThread);
    cerr << "CONTIG_AVERAGE_COVERAGE = " << averageCoverage << endl;

	if (optionBool["-fastg"]) {
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_nameIndex.csv";
		contig.printNameIndexCSV(outStream.str());
	}

	mapper->setMultiSeedLength(multiSeedLengthForShortRead);
    unsigned nowFileNumber = 0;
    for (unsigned i = 0; i < libraryMT.size(); ++i) {
		libraryMT[i][0].setAverageInsSize(0);
        int nowLibraryID = optionPairFile[nowFileNumber].libraryID;
        cerr << "[LIBRARY " << libraryIDList[i] << "]" << endl;
        // set average length and minimum insert size
        // estimate insert size
        libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));
        long minInsertion = optionMinIns.find(nowLibraryID) == optionMinIns.end() ? 0 : optionMinIns[nowLibraryID];

        mapper->contigMap.mapPairAndSaveReadLink(libraryMT[i], minInsertion, this->contigMaxK + 1, numThread);

		for (int j = 0; j < numThread; ++j) {
			fclose(libraryMT[i][j].pairFP);
			libraryMT[i][j].pairFP = NULL;
		}

        libraryMT[i][0].setInsCutoffRate(optionInsCutoffRate.find(nowLibraryID) == optionInsCutoffRate.end() ? DEFAULT_INS_CUTOFF_RATE : optionInsCutoffRate[nowLibraryID]);
        if (optionAveIns.find(nowLibraryID) != optionAveIns.end() || optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
            if (optionAveIns.find(nowLibraryID) != optionAveIns.end()) {
                libraryMT[i][0].setAverageInsSize(optionAveIns[nowLibraryID]);
                std::cerr << "Average insert size specified: AVE = " << libraryMT[i][0].getAverageInsSize() << std::endl;
            }
            if (optionSDIns.find(nowLibraryID) != optionSDIns.end()) {
                libraryMT[i][0].setSDInsSize(optionSDIns[nowLibraryID]);
            } else {
                libraryMT[i][0].setSDInsSize(static_cast<long>(static_cast<double>(libraryMT[i][0].getAverageInsSize()) / 10.0 + 0.5));
            }
        }
        nowFileNumber += numFilePerLibraryID[i];
    }

	if (longReadLibraryMT.size() > 0) {
		cerr << "[LONG_READ LIBRARY]" << endl;

		string alignerOutFilename(optionSingleArgs["-o"]);
		alignerOutFilename += "_longReadAlignment.tsv";

		string aligner = optionSingleArgs["-mapper"];

		if (optionBool["-kmer_align"]) {
			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			mapper->setMultiSeedLength(multiSeedLengthForLongRead);
			mapper->contigMap.mapLongReadAndSaveReadLink(longReadLibraryMT, this->contigMaxK, numThread);
		}
		else if (optionBool["-minialign"]) {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "minialign";

			string alignerTargetFilename(optionSingleArgs["-o"]);
			alignerTargetFilename += "_longReadTarget.fa";
			std::ostringstream targetOss;
			targetOss << "cat ";
			for (auto itr = optionMultiArgs["-c"].begin(); itr != optionMultiArgs["-c"].end(); ++itr)
				targetOss << " " << *itr;
			for (auto itr = optionMultiArgs["-b"].begin(); itr != optionMultiArgs["-b"].end(); ++itr)
				targetOss << " " << *itr;
			targetOss << " >" << alignerTargetFilename;
			system(targetOss.str().c_str());

			string alignerQueryFilename(optionSingleArgs["-o"]);
			alignerQueryFilename += "_longReadQuery.txt";
			std::ostringstream queryOss;
			queryOss << "cat ";
			for (auto itr = optionMultiArgs["-p"].begin(); itr != optionMultiArgs["-p"].end(); ++itr)
				queryOss << " " << *itr;
			queryOss << " >" << alignerQueryFilename;
			system(queryOss.str().c_str());

			execMinialign(alignerTargetFilename, alignerQueryFilename, alignerOutFilename , numThread, aligner);

			std::ostringstream rmOss;
			rmOss << "rm " << alignerTargetFilename << " " << alignerQueryFilename;
			system(rmOss.str().c_str());

			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			mapper->contigMap.readLongReadPAFfileAndSaveLink(alignerOutFilename, longReadLibraryMT, LONG_READ_MIN_ALIGNMENT_LENGTH, LONG_READ_MIN_ALIGNMENT_COVERAGE, LONG_READ_MIN_ALIGNMENT_IDENTITY, contigMaxK, numThread);
		}
		else if (optionBool["-minimap2"]) {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "minimap2";
			execMinimap2(optionMultiArgs["-c"], optionMultiArgs["-b"], optionMultiArgs["-p"], alignerOutFilename, numThread, aligner);
			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			mapper->contigMap.readLongReadPAFfileAndSaveLink(alignerOutFilename, longReadLibraryMT, LONG_READ_MIN_ALIGNMENT_LENGTH, LONG_READ_MIN_ALIGNMENT_COVERAGE, LONG_READ_MIN_ALIGNMENT_IDENTITY, contigMaxK, numThread);
		}
		else {
			if (optionSingleArgs["-mapper"].empty())
				aligner = "minimap";
			execMinimap(optionMultiArgs["-c"], optionMultiArgs["-b"], optionMultiArgs["-p"], alignerOutFilename, contigMaxK/2 , numThread, aligner);
			longReadLibraryMT[0].setAverageLength((long)((double)(longReadLibraryMT[0].getTotalLength()) / (2 * longReadLibraryMT[0].getNumPair()) + 0.5));
			mapper->contigMap.readLongReadPAFfileAndSaveLink(alignerOutFilename, longReadLibraryMT, LONG_READ_MIN_ALIGNMENT_LENGTH, LONG_READ_MIN_ALIGNMENT_COVERAGE, 0.0, contigMaxK, numThread);
		}

		for (int j = 0; j < numThread; ++j) {
			fclose(longReadLibraryMT[j].pairFP);
			longReadLibraryMT[j].pairFP = NULL;
		}

		vector<long> insSizeDistribution;
		longReadLibraryMT[0].readInsertSizeFile(insSizeDistribution);
		longReadLibraryMT[0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_longReadLibrary" << "_readDistribution.tsv";
		longReadLibraryMT[0].printInsertSizeFreq(outStream.str());
		cerr << "[LONG_READ_LIBRARY " << 1 << "]\nAVE_READ_LENGTH = " << longReadLibraryMT[0].getAverageInsSize()
			 << ", SD_READ_LENGTH = " << longReadLibraryMT[0].getSDInsSize() << endl;
	}

	if (tagLibraryMT.size() > 0) {
		mapper->setMultiSeedLength(multiSeedLengthForShortRead);

		cerr << "[TAG LIBRARY]" << endl;
		mapper->contigMap.mapTagPairMT(tagLibraryMT, numThread);

		for (int j = 0; j < numThread; ++j) {
			fclose(tagLibraryMT[j].pairFP);
			tagLibraryMT[j].pairFP = NULL;
		}
	}

	//added by ouchi
	if (HiCLibraryMT.size() > 0) {
		mapper->setMultiSeedLength(multiSeedLengthForShortRead);

		cerr << "[HIC LIBRARY]" << endl;
		mapper->contigMap.mapHiCPairMT(HiCLibraryMT, numThread);

		for (int j = 0; j< numThread; ++j) {
			fclose(HiCLibraryMT[j].pairFP);
			HiCLibraryMT[j].pairFP = NULL;
		}
	}
	//

	if (libraryMT.size() > 0) {
		pairedDBG.setAllLibraryMT(&libraryMT);
		pairedDBG.setTargetLibraryIndex(0);
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);
		pairedDBG.setContigNameIndex(contig.nameIndex);

		if (optionSingleArgs["-unlink_list"].size() > 0)
			pairedDBG.setContigUnlinkFlags(optionSingleArgs["-unlink_list"]);

		vector<long> insSizeDistribution;
		libraryMT[0][0].readInsertSizeFile(insSizeDistribution);
		pairedDBG.insertSizeDistribution(libraryMT[0], insSizeDistribution, numThread);

		vector<long> seqLengths;
		pairedDBG.scaffoldLengthList(seqLengths);

		if (libraryMT[0][0].getAverageInsSize() == 0)
			libraryMT[0][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << 1 << "_insFreq.tsv";
		libraryMT[0][0].printInsertSizeFreq(outStream.str());
		cerr << "[LIBRARY " << 1 << "]\nAVE_INS = " << libraryMT[0][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[0][0].getSDInsSize() << endl;
	}
	else {
		pairedDBG.initScaffolding(contig.coverage, mapper->contigMap, averageCoverage, bubbleThreshold);
		pairedDBG.setContigName(contig.name);
	}

    pairedDBG.setContigMaxK(this->contigMaxK);
    pairedDBG.setMinOverlap(this->contigMaxK - 1);
    pairedDBG.saveOverlap(mapper->contigMap, this->contigMaxK - 1, this->contigMaxK, numThread);
    pairedDBG.classifyNode();

	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setLongReadLibraryMT(&longReadLibraryMT);
	}

	if (tagLibraryMT.size() > 0) {
		pairedDBG.setTagLibraryMT(&tagLibraryMT);
		pairedDBG.countMappedTagForEachContig(numThread);
	}

	//added by ouchi
	if (HiCLibraryMT.size() > 0) {
		pairedDBG.setHiCLibraryMT(&HiCLibraryMT);
	}
	//

    cerr << "destructing mapper objects..." << std::endl;
}


void SolveDBG::readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread)
{
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(static, 1)
    for (int i = -4; i < static_cast<int>(libraryMT.size()); ++i) { //edidted by ouchi (from i=-3)
        try {
            // load contig file
            if (i == -4) {  //edited by ouchi (from i == -3)
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-c"][i]);

				long j = contig.numSeq;
                for (unsigned i = 0; i < optionMultiArgs["-b"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-b"][i]);

				pairedDBG.setNumInputBubbleContig(contig.numSeq - j);
				contig.setNameIndex();


                this->contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
                this->contigReadLength = contig.getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
                if (this->contigReadLength == 0)
                    this->contigReadLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;

				if (optionSingleArgs["-e"] == "")
					averageCoverage = contig.calculateAverageCoverageExcludingOutlier(contig.getSeqLengthMedian());
				else
					averageCoverage = atof(optionSingleArgs["-e"].c_str());

                mapper->setContigMap(contig);
                mapper->makeKmerTableContigMap();
            }
			else if (i == -3) { //edited by ouchi (from i==-2)
                if (optionMultiArgs["-p"].size() == 0)
					continue;

				vector<string> filenames(optionMultiArgs["-p"]);

				longReadLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					longReadLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-p"].size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-p"][j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleMT(longReadLibraryMT, optionMultiArgs["-p"][j], numThread, false, isFastq, true);
				}
            }
			else if (i == -2) { //edited by ouchi (from i== -1)
                if (optionMultiArgs["-x"].size() == 0 && optionMultiArgs["-X"].size() == 0)
					continue;

				vector<string> filenames(optionMultiArgs["-x"]);
				filenames.insert(filenames.end(), optionMultiArgs["-X"].begin(), optionMultiArgs["-X"].end());

				std::unordered_map<string, int> tagStringConverter;
				setTagStringConverter(filenames, tagStringConverter);

				tagLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					tagLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-x"].size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-x"][j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleTaggedMT(tagLibraryMT, optionMultiArgs["-x"][j], numThread, false, isFastq, true, tagStringConverter);
				}

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-X"].size()); j += 2) {
					platanus::FILETYPE fileFormat1 = checkFileFormat(optionMultiArgs["-X"][j]);
					platanus::FILETYPE fileFormat2 = checkFileFormat(optionMultiArgs["-X"][j + 1]);
					if (fileFormat1 != fileFormat2) {
						throw platanus::FormatError("Different file type in paired-file (-X).");
					}
					bool isFastq;
					switch (fileFormat1) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaPairTaggedMT(tagLibraryMT, optionMultiArgs["-X"][j] , optionMultiArgs["-X"][j + 1], numThread, false, isFastq, tagStringConverter);
				}
            }
			//added by ouchi
			else if (i == -1) {
                if (optionMultiArgs["-hic"].size() == 0 && optionMultiArgs["-HIC"].size() == 0)
					continue;

				HiCLibraryMT.resize(numThread);
				for (long j = 0; j < numThread; ++j)
					HiCLibraryMT[j].pairFP = platanus::makeTemporaryFile();

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-hic"].size()); ++j) {
					platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-hic"][j]);
					bool isFastq;
					switch (fileFormat) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaSingleMT(HiCLibraryMT, optionMultiArgs["-hic"][j], numThread, false, isFastq);
				}

				for (unsigned j = 0; j < static_cast<long>(optionMultiArgs["-HIC"].size()); j += 2) {
					platanus::FILETYPE fileFormat1 = checkFileFormat(optionMultiArgs["-HIC"][j]);
					platanus::FILETYPE fileFormat2 = checkFileFormat(optionMultiArgs["-HIC"][j + 1]);
					if (fileFormat1 != fileFormat2) {
						throw platanus::FormatError("Different file type in paired-file (-HIC).");
					}
					bool isFastq;
					switch (fileFormat1) {
					case platanus::UNKNOWN:
						throw platanus::FormatError();
						break;
					case platanus::FASTA:
						isFastq = false;
						break;
					default:
						isFastq = true;
						break;
					}
					ReadFastaPairMT(HiCLibraryMT, optionMultiArgs["-HIC"][j] , optionMultiArgs["-HIC"][j + 1], numThread, false, isFastq);
				}
			}
			////
			else {
                unsigned nowFileNumber = 0;
                libraryMT[i].resize(numThread);
                for (int j = 0; j < i; ++j)
                    nowFileNumber += numFilePerLibraryID[j];

                for (int j = 0; j < numThread; ++j) {
                    libraryMT[i][j].pairFP = platanus::makeTemporaryFile();
                }
                for (int j = 0; j < numFilePerLibraryID[i]; ++j) {
                    bool isFastq, isMate;
                    isMate = optionPairFile[nowFileNumber+j].libraryType == "-op"
                           || optionPairFile[nowFileNumber+j].libraryType == "-OP" ? true : false;
                    if (optionPairFile[nowFileNumber+j].libraryType == "-ip" || optionPairFile[nowFileNumber+j].libraryType == "-op") {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        switch (fileFormat) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaSingleMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, numThread, isMate, isFastq);
                    }
					else {
                        platanus::FILETYPE fileFormat1 = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        platanus::FILETYPE fileFormat2 = checkFileFormat(optionPairFile[nowFileNumber+j].fileSecond);
                        if (fileFormat1 != fileFormat2) {
                            throw platanus::FormatError("Different file type in paired-file.");
                        }
                        switch (fileFormat1) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaPairMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, optionPairFile[nowFileNumber+j].fileSecond, numThread, isMate, isFastq);
                    }
                }
            }
        } catch (platanus::ErrorBase &e) {
            e.showErrorMessage();
            exit(e.getID());
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// destroy graph
//////////////////////////////////////////////////////////////////////////////////////
void SolveDBG::outputAndAfterTreatment(void)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    string componentFilename(outFilename);
    outFilename += "_solvedContig.fa";
    componentFilename += "_solvedContigComponent.tsv";
    pairedDBG.cutAndPrintSeq(this->contigMaxK, this->contigReadLength, outFilename, componentFilename);
}

void SolveDBG::outputGraph(const char *suffix)
{
    cerr << "writing scaffold files..." << endl;
    string outFilename(optionSingleArgs["-o"]);
    outFilename += suffix;
    pairedDBG.printResultSeq(outFilename);
}

void SolveDBG::updateAndWriteInsertSize(const long libraryIndex)
{
	pairedDBG.updateInsertLengthFP(libraryMT[libraryIndex], numThread);
	vector<long> insSizeDistribution;
	libraryMT[libraryIndex][0].readInsertSizeFile(insSizeDistribution);
	pairedDBG.insertSizeDistribution(libraryMT[libraryIndex], insSizeDistribution, numThread);
	if (libraryIndex > 0)
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, libraryMT[libraryIndex - 1][0].getAverageInsSize(), platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);
	else
		libraryMT[libraryIndex][0].estimateInsSize(insSizeDistribution, 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

	std::ostringstream outStream;
	outStream << optionSingleArgs["-o"] << "_lib" << (libraryIndex + 1) << "_insFreq.tsv";
	printInsertSizeFreq(outStream.str(), insSizeDistribution);
}

void SolveDBG::extendConsensus(const long numOuterIteration, const bool divisionFlag)
{
	pairedDBG.clearContigPreviousParentNodeID();

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);

	pairedDBG.setMinOverlap(minOverlapForScaffolding);

	for (long outerIteration = 0; outerIteration < numOuterIteration; ++outerIteration) {
		long numIteration = 2;
		for (long iteration = 0; iteration < numIteration; ++iteration) {
			pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

			for (long libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
				pairedDBG.setTargetLibraryIndex(libraryIndex);

				if (libraryMT[libraryIndex][0].getAverageInsSize() <= 0)
					updateAndWriteInsertSize(libraryIndex);

				if (iteration == 0)
					pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
				else
					pairedDBG.setMinLink(minLink);

				cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
				for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
					pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
					pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
					cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

					pairedDBG.makeGraph(numThread);
					pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
//					pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
					if (iteration > 0)
						pairedDBG.deleteErroneousEdgeIterative(numThread);
					pairedDBG.deleteRepeatEdge();
					pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
					pairedDBG.makeScaffold();

					if (numIteration >= 4 && outerIteration < 2)
						pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);
				}
			}

			if (divisionFlag)
				pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
		}


		if (libraryMT.size() > 0) {
			for (long iteration = 0; iteration < numIteration; ++iteration) {
				pairedDBG.setTargetLibraryIndex(libraryMT.size() - 1);
				long linkThreashold = (1 + iteration) * minLink;
				pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

				pairedDBG.setTolerenceFactor(MAX_TOL_FACTOR);
				pairedDBG.setCutoffLengthFactor(MAX_LONG_READ_LENGTH_CUTOFF_FACTOR * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR);
				pairedDBG.setMinLink(linkThreashold);
				pairedDBG.setCutoffLength(2.0 * MAX_TOL_FACTOR * libraryMT.back()[0].getSDInsSize());
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				pairedDBG.makeGraphAllLibraries(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();

				if (numIteration >= 4 && outerIteration < 2)
					pairedDBG.divideGappedNode((outerIteration + 1) * this->contigMaxK);

				if (divisionFlag)
					pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
			}
		}


		if (longReadLibraryMT.size() > 0) {
			pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);
			pairedDBG.setTolerence(0.1 * longReadLibraryMT[0].getAverageInsSize());

			for (unsigned cutoffFactor = MIN_LONG_READ_LENGTH_CUTOFF_FACTOR; cutoffFactor <= MAX_LONG_READ_LENGTH_CUTOFF_FACTOR; ++cutoffFactor) {
				pairedDBG.setCutoffLength(cutoffFactor * LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR *  longReadLibraryMT[0].getAverageInsSize());
				pairedDBG.setMinLink(minLink);
				pairedDBG.makeGraph(numThread);
				pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
//				pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
				pairedDBG.deleteErroneousEdgeNumLinkRateIterative(numThread);
				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();
			}

			if (divisionFlag)
				pairedDBG.divideErroneousNodeBaseLevel(0, libraryMT.size(), this->contigReadLength, false, true, false, numThread);
		}
	}

	pairedDBG.setMinOverlap(this->contigMaxK - 1);
}

void SolveDBG::execMinialign(const string targetFilename, const string &readFilename, const string &outFilename, const long numThread, const string minialignExecutable)
{
	std::ostringstream oss;

	oss << minialignExecutable <<  " -x pacbio -m 0 -O paf" << " -t " << numThread << " " << targetFilename << " " << readFilename << " >" << outFilename;

	std::cerr << "Executing minialign ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minialign fineshed." << endl;
}

void SolveDBG::execMinimap(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long minAlignmentLength, const long numThread, const string minimapExecutable)
{
	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimapExecutable <<  " -t " << numThread <<  " -L " << minAlignmentLength << " - ";

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
		oss << " " << *itr;
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename;

	std::cerr << "Executing minimap ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	cerr << "minimap fineshed." << endl;
}

void SolveDBG::execMinimap2(const vector<string> &contigFilenames, const vector<string> &bubbleFilenames, const vector<string> &readFilenames, const string &outFilename, const long numThread, const string minimap2Executable)
{
	std::ostringstream oss;

	oss << "cat ";
	for (auto itr = contigFilenames.begin(); itr != contigFilenames.end(); ++itr)
		oss << " " << *itr;
	for (auto itr = bubbleFilenames.begin(); itr != bubbleFilenames.end(); ++itr)
		oss << " " << *itr;

	oss << " | " << minimap2Executable << " -c " << " -t " << numThread << " - ";
	if (optionBool["-minimap2_sensitive"])
		oss << " -p 0";


	bool bzip2Flag = false;
	char bzip2TempFileName[platanus::globalTmpFileDir.size() + 7 + 1];
	platanus::FILECOMPRESSION format;

	for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
		format = platanus::checkFileCompression(*itr);
		if (format == platanus::FILECOMPRESSION::BZIP2) {
			bzip2Flag = true;
			break;
		}
	}

	if (bzip2Flag) {
		strcpy(bzip2TempFileName, platanus::globalTmpFileDir.c_str());
		strcat(bzip2TempFileName, "/XXXXXX");

        int fd = mkstemp(bzip2TempFileName);
        if (fd == -1) {
            throw platanus::TMPError();
		}
        FILE *fp = fdopen(fd, "w+");
		fclose(fp);

		std::ostringstream bzip2Oss;
		bzip2Oss << "bzip2 -cd ";

		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr) {
			format = platanus::checkFileCompression(*itr);
			if (format == platanus::FILECOMPRESSION::BZIP2)
				bzip2Oss << " " << *itr;
			else
				oss << " " << *itr;
		}

		bzip2Oss << " >"  << bzip2TempFileName;
		if (system(bzip2Oss.str().c_str()) != 0) {
			throw platanus::AlignerError();
		}

		oss << " " << bzip2TempFileName;
	}
	else {
		for (auto itr = readFilenames.begin(); itr != readFilenames.end(); ++itr)
			oss << " " << *itr;
	}
	oss << " | cut -f1-11 ";
	oss << " >" << outFilename;

	std::cerr << "Executing minimap2 ..." << std::endl;
	std::cerr << oss.str() << std::endl << std::endl << std::endl;

	if (system(oss.str().c_str()) != 0) {
		throw platanus::AlignerError();
	}

	if (bzip2Flag) {
        unlink(bzip2TempFileName);
	}

	cerr << "minimap2 finished." << endl;
}




void SolveDBG::extendConsensusToEstimateInsertSize(void)
{
	pairedDBG.clearContigPreviousParentNodeID();

	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.setMinLink(minLink);
	pairedDBG.makeGraph(numThread);
	pairedDBG.joinUnambiguousNodePairIterative(numThread);

	pairedDBG.makeGraph(numThread);
	pairedDBG.setOppositeBubbleContigIDOverlapped(numThread);
	pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);
	pairedDBG.clearEdges();
	pairedDBG.makeScaffold();
	pairedDBG.joinUnambiguousNodePairIterative(numThread);


	pairedDBG.setMinOverlap(minOverlapForScaffolding);

	pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE | PairedDBG::LENGTH_CUTOFF_MODE | PairedDBG::PREVIOUS_DIVISION_AWARE_MODE);

	for (long iteration = 0; iteration < 2; ++iteration) {
		for (unsigned libraryIndex = 0; libraryIndex < libraryMT.size(); ++libraryIndex) {
			pairedDBG.setTargetLibraryIndex(libraryIndex);
			updateAndWriteInsertSize(libraryIndex);

			if (iteration == 0)
				pairedDBG.setMinLink(std::max(minLink, pairedDBG.estimateLink()));
			else
				pairedDBG.setMinLink(minLink);

			cerr << "[LIBRARY " << libraryIndex + 1 << "]\nAVE_INS = " << libraryMT[libraryIndex][0].getAverageInsSize() << ", SD_INS = " << libraryMT[libraryIndex][0].getSDInsSize() << endl;
			for (unsigned tolerenceFactor = MIN_TOL_FACTOR; tolerenceFactor <= MAX_TOL_FACTOR; ++tolerenceFactor) {
				pairedDBG.setCutoffLength(2.0 * tolerenceFactor * libraryMT[libraryIndex][0].getSDInsSize());
				pairedDBG.setTolerence(pairedDBG.getCutoffLength() / 2);
				cerr << "TOLERENCE_LEVEL_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;

				pairedDBG.makeGraph(numThread);
				pairedDBG.setOppositeBubbleContigIDGapped(numThread);
				pairedDBG.deleteSecondaryBubbleNodeAndEdge(numThread);

				pairedDBG.deleteErroneousEdgeNumTagRateIterative(numThread);
//				pairedDBG.deleteErroneousEdgebyHiCIterative(numThread); //added by ouchi
				if (iteration > 0)
					pairedDBG.deleteErroneousEdgeIterative(numThread);

				pairedDBG.deleteRepeatEdge();
				pairedDBG.detectRepeat(2 * pairedDBG.getHeteroCoverage());
				pairedDBG.makeScaffold();
			}
		}

		pairedDBG.divideErroneousNode(minLink, PairedDBG::DIVISION_MODE::MIS, numThread);
	}

	pairedDBG.setMinOverlap(this->contigMaxK - 1);
}


void SolveDBG::generateGraphFastg()
{
	pairedDBG.setMode(PairedDBG::OVERLAP_MODE);
	pairedDBG.makeGraph(numThread);

	std::ostringstream outStream;
	outStream << optionSingleArgs["-o"] << "_overlap.fastg";
	pairedDBG.outputFastg(outStream.str());
	pairedDBG.clearEdges();


	for (unsigned i = 0; i < libraryMT.size(); ++i) {
		pairedDBG.setTargetLibraryIndex(i);

		pairedDBG.setMode(PairedDBG::PAIRED_END_LINK_MODE);
		if (i > 0 && libraryMT[i][0].getAverageInsSize() == 0) {
			vector<long> insSizeDistribution;
			libraryMT[i][0].readInsertSizeFile(insSizeDistribution);
			pairedDBG.insertSizeDistribution(libraryMT[i], insSizeDistribution, numThread);
			libraryMT[i][0].estimateInsSize(insSizeDistribution, libraryMT[i - 1][0].getAverageInsSize(), platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

			std::ostringstream outStream;
			outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_insFreq.tsv";
			printInsertSizeFreq(outStream.str(), insSizeDistribution);
		}
		cerr << "[LIBRARY " << i + 1 << "]\nAVE_INS = " << libraryMT[i][0].getAverageInsSize()
			 << ", SD_INS = " << libraryMT[i][0].getSDInsSize() << endl;

		unsigned tolerenceFactor = MAX_TOL_FACTOR;
		pairedDBG.setTolerence(tolerenceFactor * libraryMT[i][0].getSDInsSize());
		cerr << "TOLERENCE_OF_CONTIG_OVERLAP = " << pairedDBG.getTolerence() << endl;
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);

		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_paired_end.fastg";
		pairedDBG.outputFastg(outStream.str());
		}
		pairedDBG.clearEdges();


		pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::PAIRED_END_LINK_MODE);
		pairedDBG.makeGraph(numThread);
		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_lib" << (i + 1) << "_paired_end_overlap.fastg";
		pairedDBG.outputFastg(outStream.str());
		}
		pairedDBG.clearEdges();
	}


	if (longReadLibraryMT.size() > 0) {
		pairedDBG.setMode(PairedDBG::LONG_READ_LINK_MODE);
		pairedDBG.setTolerence(2 * this->contigMaxK);
		pairedDBG.setMinLink(minLink);
		pairedDBG.makeGraph(numThread);
		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_long_read.fastg";
		pairedDBG.outputFastg(outStream.str());
		pairedDBG.clearEdges();
		}

		pairedDBG.setMode(PairedDBG::OVERLAP_MODE | PairedDBG::LONG_READ_LINK_MODE);
		pairedDBG.makeGraph(numThread);
		{
		std::ostringstream outStream;
		outStream << optionSingleArgs["-o"] << "_long_read_overlap.fastg";
		pairedDBG.outputFastg(outStream.str());
		}
		pairedDBG.clearEdges();
	}
}
