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

#include "phase.h"
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Phase::Phase()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-k"] = "1";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    optionSingleArgs["-i"] = "2";
    optionSingleArgs["-tmp"] = ".";

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();
    optionMultiArgs["-hic"] = vector<string>(); //added by ouchi
    optionMultiArgs["-HIC"] = vector<string>(); //added by ouchi

    optionSingleArgs["-mapper"] = "";
    optionBool["-minimap2_sensitive"] = false;
    optionBool["-minimap2"] = false;
    optionBool["-minialign"] = false;
    optionBool["-kmer_align"] = false;

    optionBool["-reduce_redundancy"] = true;

    optionBool["-aggressive"] = false;
    optionBool["-no_partial"] = false;
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Phase::usage(void) const
{

    std::cerr << "\nUsage: platanus_allee phase [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)\n"
              << "    -hic PAIR1 [PAIR2 ...]             : HiC-reads files (paired-ends) (interleaved, fasta or fastq)\n" //added by ouchi
              << "    -HIC FWD1 REV1 [FWD2 REV2 ...]     : HiC-reads files (paired-ends) (separate forward and reverse files, fasta or fastq)\n" //added by ouchi
              << "    -i INT                             : number of iterations (default " << optionSingleArgs.at("-i") << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -k INT                             : minimum number of links to phase variants (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1] << " " << optionMultiArgs.at("-s")[2]  << ")\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -mapper FILE                       : path of mapper executable file (default minimap2, only effective with -p option)\n"
              << "    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)\n"
              << "    -aggressive                        : aggressively extend scaffolds (default off)\n"
//              << "    -no_partial                        : not close gaps partially, i.e. only close ones completely (default, off)\n"
//              << "    -reduce_redundancy                 : reduce redundant sequences that exactly matche others (default, off)\n"
//              << "    -minialign                         : use minialign insterd of minimap2 (default) for alignment of long reads\n"
//              << "    -kmer_align                        : use built-in alignment based on k-mer match insterd of minimap (default) for long reads\n"
              << "\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x, -X, -hic and -HIC options.\n" //edited by ouchi
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_allPhasedScaffold.fa (including sequences below)\n"
              << "    PREFIX_primaryBubble.fa\n"
              << "    PREFIX_secondaryBubble.fa\n"
              << "    PREFIX_nonBubbleHetero.fa\n"
              << "    PREFIX_nonBubbleOther.fa\n"
              << "    PREFIX_consensusInput.fa (for \"consensus\" command (-c))\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


void Phase::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	string intermediateDirectoryName = optionSingleArgs["-o"];
	intermediateDirectoryName += "_intermediateResults";
	this->setDirectoryName(intermediateDirectoryName);
	this->createDirectory();

    unsigned long numIterate = atoi(optionSingleArgs["-i"].c_str());
    if (!(optionMultiArgs["-p"].empty())) {
		++numIterate;
	}

    for (unsigned long times = 1; times <= numIterate; ++times) {
		std::ostringstream oss;
		oss << intermediateDirectoryName << "/round" << times;
		this->setDirectoryName(oss.str());

        this->createDirectory();
		this->execSolveDBG(times, numIterate);
		this->execGapClose(times, numIterate);
    }

	this->moveAndConcatenateFinalRoundResult(intermediateDirectoryName, numIterate);

	std::cerr << "phase completed!" << std::endl;
}


void Phase::setDirectoryName(const string newDirectoryName)
{
	this->previousDirectoryName = this->directoryName;
	this->directoryName = newDirectoryName;
}


void Phase::createDirectory(void) const
{
    struct stat st;

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
            throw platanus::CreateDirError(directoryName);
        }
    }
}


void Phase::execSolveDBG(const unsigned long times, const unsigned long numIterate)
{
	vector<string> boolArgOptionList = {"-minimap2_sensitive", "-aggressive", "-reduce_redundancy"};
	vector<string> singleArgOptionList = {"-t", "-tmp", "-l", "-k", "-mapper"};
	vector<string> multiArgOptionList = {"-x", "-X", "-hic", "-HIC", "-s"}; //edited by ouchi
	if (times > 1) {
		multiArgOptionList.push_back("-p");
	}

	string inputPrefix = this->previousDirectoryName;
	inputPrefix += "/";
	inputPrefix += optionSingleArgs["-o"];

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];


    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG";

	if (times == 1) {
		multiArgOptionList.push_back("-c");
		multiArgOptionList.push_back("-b");
	}
	else {
        oss << " -c " << inputPrefix << "_gapClosed_nonBubbleOther.fa"
			<< " -b " << inputPrefix << "_gapClosed_primaryBubble.fa" << " " << inputPrefix << "_gapClosed_secondaryBubble.fa" << " " << inputPrefix << "_gapClosed_primaryFork.fa" << " " << inputPrefix << "_gapClosed_secondaryFork.fa";
	}

	for (auto optItr = boolArgOptionList.begin(); optItr != boolArgOptionList.end(); ++optItr) {
		if (optionBool[*optItr])
			oss << " " << *optItr;
	}

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

	for (auto optItr = singleArgOptionList.begin(); optItr != singleArgOptionList.end(); ++optItr) {
		if (!(optionSingleArgs[*optItr].empty())) {
			oss << " " << *optItr << " " << optionSingleArgs[*optItr];
		}
	}

	for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
		if (!(optionMultiArgs[*optItr].empty())) {
			oss << " " << *optItr;
			for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

	oss << " -o " << outputPrefix;
    oss << " 2>" << outputPrefix << ".solveDBGLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }
}


void Phase::execGapClose(const unsigned long times, const unsigned long numIterate)
{
	string inputPrefix = this->directoryName;
	inputPrefix += "/";
	inputPrefix += optionSingleArgs["-o"];

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];

	vector<string> singleArgOptionList = {"-t", "-tmp"};
	vector<string> multiArgOptionList = {"-s"};
	vector<string> boolArgOptionList = {"-no_partial", "-reduce_redundancy"};

    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close";

	oss << " -c"
		<< " " << inputPrefix << "_primaryBubble.fa"
		<< " " << inputPrefix << "_secondaryBubble.fa"
		<< " " << inputPrefix << "_primaryFork.fa"
		<< " " << inputPrefix << "_secondaryFork.fa"
		<< " " << inputPrefix << "_nestedBubble.fa"
		<< " " << inputPrefix << "_nonBubbleOther.fa";

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

	for (auto optItr = singleArgOptionList.begin(); optItr != singleArgOptionList.end(); ++optItr) {
		if (!(optionSingleArgs[*optItr].empty())) {
			oss << " " << *optItr << " " << optionSingleArgs[*optItr];
		}
	}

	for (auto optItr = multiArgOptionList.begin(); optItr != multiArgOptionList.end(); ++optItr) {
		if (!(optionMultiArgs[*optItr].empty())) {
			oss << " " << *optItr;
			for (auto argItr = optionMultiArgs[*optItr].begin(); argItr != optionMultiArgs[*optItr].end(); ++argItr) {
				oss << " " << *argItr;
			}
		}
	}

	for (auto optItr = boolArgOptionList.begin(); optItr != boolArgOptionList.end(); ++optItr) {
		if (optionBool[*optItr])
			oss << " " << *optItr;
	}

	oss << " -o " << outputPrefix;
    oss << " 2>" << outputPrefix << ".gapCloseLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }
}

void Phase::moveAndConcatenateFinalRoundResult(const string intermediateDirectoryName, const unsigned long times)
{
    std::ostringstream prefixOss;
	prefixOss << intermediateDirectoryName << "/round" << times << "/" << optionSingleArgs["-o"];
	string finalRoundPrefix = prefixOss.str();

    std::ostringstream cmdOss;
    cmdOss << "mv " << finalRoundPrefix << "_gapClosed_primaryBubble.fa" << " " << optionSingleArgs["-o"] << "_primaryBubble.fa" << "; "
		<< "mv " << finalRoundPrefix << "_gapClosed_secondaryBubble.fa" << " " << optionSingleArgs["-o"] << "_secondaryBubble.fa" << "; "
		<< "mv " << finalRoundPrefix << "_gapClosed_nonBubbleOther.fa" << " " << optionSingleArgs["-o"] << "_nonBubbleOther.fa" << "; "
		<< "cat " << optionSingleArgs["-o"] << "_primaryBubble.fa" << " "
		<< optionSingleArgs["-o"] << "_secondaryBubble.fa" << " "
		<< finalRoundPrefix << "_gapClosed_primaryFork.fa" << " "
		<< finalRoundPrefix << "_gapClosed_secondaryFork.fa" << " "
		<< finalRoundPrefix << "_gapClosed_nestedBubble.fa" << " "
		<< optionSingleArgs["-o"] << "_nonBubbleOther.fa" << " "
		<< ">" << optionSingleArgs["-o"] <<  "_allPhasedScaffold.fa" << "; "

		<< "cat " << finalRoundPrefix << "_gapClosed_primaryFork.fa" << " "
		<< finalRoundPrefix << "_gapClosed_secondaryFork.fa" << " "
		<< finalRoundPrefix << "_gapClosed_nestedBubble.fa" << " "
		<< ">" << optionSingleArgs["-o"] <<  "_nonBubbleHetero.fa" << "; "

		<< "cat " << optionSingleArgs["-o"] << "_primaryBubble.fa" << " "
		<< finalRoundPrefix << "_gapClosed_primaryFork.fa" << " "
		<< optionSingleArgs["-o"] << "_nonBubbleOther.fa" << " "
		<< ">" << optionSingleArgs["-o"] <<  "_consensusInput.fa";

    if (system(cmdOss.str().c_str()) != 0) {
        throw platanus::CreateLinkError();
    }
}
