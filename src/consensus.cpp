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

#include "consensus.h"
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Consensus::Consensus()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-e"] = "";
    optionSingleArgs["-v"] = "32";
    optionSingleArgs["-k"] = "1";
    optionSingleArgs["-l"] = "3";
    optionSingleArgs["-u"] = "0";
    optionSingleArgs["-L"] = "200000";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-b"] = vector<string>();
    optionMultiArgs["-p"] = vector<string>();
    optionMultiArgs["-x"] = vector<string>();
    optionMultiArgs["-X"] = vector<string>();

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    optionMultiArgs["-S"] = vector<string>(1);
    optionMultiArgs["-S"][0] = "20";


    optionBool["-no_scaffold"] = false;
    optionBool["-unphase"] = false;

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";

    optionSingleArgs["-mapper"] = "";
    optionBool["-minimap2_sensitive"] = false;
    optionSingleArgs["-unlink_list"] = "";

    optionBool["-minimap2"] = false;
    optionBool["-minialign"] = false;
    optionBool["-kmer_align"] = false;

    optionBool["-aggressive"] = false;
    optionBool["-no_partial"] = false;
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Consensus::usage(void) const
{

    std::cerr << "\nUsage: platanus_allee consensus [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : input_scaffolds (fasta format)\n"
//              << "    -b FILE1 [FILE2 ...]               : bubble_seq_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -p PAIR1 [PAIR2 ...]               : long-read file (PacBio, Nanopore) (fasta or fastq)\n"
              << "    -x PAIR1 [PAIR2 ...]               : linked-reads files (paired-ends, 10x Genomics) (interleaved, fasta or fastq)\n"
              << "    -X FWD1 REV1 [FWD2 REV2 ...]       : linked-reads files (paired-ends, 10x Genomics) (separate forward and reverse files, fasta or fastq)\n"
//              << "    -n{INT} INT                        : lib_id minimum_insert_size\n"
//              << "    -a{INT} INT                        : lib_id average_insert_size\n"
//              << "    -d{INT} INT                        : lib_id SD_insert_size\n"
              << "    -e FLOAT                           : coverage depth of homozygous region (default auto)\n"
              << "    -L INT                             : maximum fragment length of tag (10x Genomics) (default " << optionSingleArgs.at("-L") << ")\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1] << " " << optionMultiArgs.at("-s")[2]  << ")\n"
              << "    -l INT                             : minimum number of links to scaffold (default " << optionSingleArgs.at("-l") << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "    -mapper FILE                       : path of mapper executable file (default, minimap, only effective with -p option)\n"
              << "    -minimap2_sensitive                : sensitive mode for minimap2 (default, off; only effective with -p option)\n"
              << "    -unlink_list                       : reduce redundant sequences that exactly matche others (default, off)\n"
//              << "    -minimap2                          : use minimap2 insterd of minimap (default) for alignment of long reads\n"
//              << "    -minialign                         : use minialign insterd of minimap (default) for alignment of long reads\n"
//              << "    -kmer_align                        : use built-in alignment based on k-mer match insterd of minimap (default) for long reads\n"
              << "    -aggressive                        : aggressively extend scaffolds (default off)\n"
              << "    -no_partial                        : not close gaps partially, i.e. only close ones completely (default, off)\n"
              << "\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op, -OP, -p, -x and -X options.\n"
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_consensusScaffold.fa (final)\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


void Consensus::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	std::cerr << "executing platanus_allee solve_DBG internally ..." << std::endl;

	string intermediateDirectoryName = optionSingleArgs["-o"];
	intermediateDirectoryName += "_consensusIntermediateResults";
	this->setDirectoryName(intermediateDirectoryName);
	this->createDirectory();

	this->execSolveDBG();
	this->execGapClose();

    std::ostringstream cmdOss;
    cmdOss << "mv "  << intermediateDirectoryName << "/" << optionSingleArgs["-o"]  << "_gapClosed_consensusScaffold.fa " << optionSingleArgs["-o"] << "_consensusScaffold.fa";
    system(cmdOss.str().c_str());

	std::cerr << "consensus completed!" << std::endl;
}


void Consensus::setDirectoryName(const string newDirectoryName)
{
	this->directoryName = newDirectoryName;
}


void Consensus::createDirectory(void) const
{
    struct stat st;

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
            throw platanus::CreateDirError(directoryName);
        }
    }
}


void Consensus::execSolveDBG(void)
{
	vector<string> boolArgOptionList = {"-minimap2_sensitive", "-aggressive"};
	vector<string> singleArgOptionList = {"-t", "-tmp", "-e", "-l", "-L", "-mapper", "-unlink_list"};
	vector<string> multiArgOptionList = {"-c", "-b", "-p", "-x", "-X", "-s"};

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];

    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -unphase";

    for (auto itr = optionPairFile.begin(); itr != optionPairFile.end(); ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (!(itr->fileSecond.empty())) {
            oss << " " << itr->fileSecond;
        }
    }

	for (auto optItr = boolArgOptionList.begin(); optItr != boolArgOptionList.end(); ++optItr) {
		if (optionBool[*optItr])
			oss << " " << *optItr;
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
    oss << " 2>" << outputPrefix << ".consensusSolveDBGLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }
}


void Consensus::execGapClose(void)
{
	vector<string> boolArgOptionList = {"-no_partial"};
	vector<string> singleArgOptionList = {"-t", "-tmp"};
	vector<string> multiArgOptionList = {"-s"};

	string inputPrefix = this->directoryName;
	inputPrefix += "/";
	inputPrefix += optionSingleArgs["-o"];

	string outputPrefix = this->directoryName;
	outputPrefix += "/";
	outputPrefix += optionSingleArgs["-o"];

    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close -reduce_redundancy"
		<< " -c " << inputPrefix << "_consensusScaffold.fa";

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
    oss << " 2>" << outputPrefix  << ".consensusGapCloseLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }
}
