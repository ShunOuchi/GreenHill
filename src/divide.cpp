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

#include "divide.h"
#include <sys/stat.h>
#include <sys/types.h>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Divide::Divide()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    pairedEndSingleFileType.push_back("-ip");
    pairedEndSingleFileType.push_back("-op");
    pairedEndPairFileType.push_back("-IP");
    pairedEndPairFileType.push_back("-OP");
    optionSingleArgs["-tmp"] = ".";
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Divide::usage(void) const
{

    std::cerr << "\nUsage: platanus_allee divide [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (interleaved file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (interleaved, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (separate forward and reverse files, fasta or fastq)\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1] << " " << optionMultiArgs.at("-s")[2]  << ")\n"
              << "    -t INT                             : number of threads (<= " << optionSingleArgs.at("-t") << ", default 1)\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n"
              << "\n\n"

              << "Input format:\n"
              << "    Uncompressed and compressed (gzip or bzip2) files are accepted for -c, -ip, -IP, -op and -OP options.\n"
              << "\n\n"

              << "Outputs:\n"
              << "    PREFIX_divided.fa\n"
              << "    PREFIX_dividedComponent.bed\n"
              << "\n"
              << "Note, PREFIX is specified by -o\n"
              << std::endl;
}


void Divide::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

	std::cerr << "executing platanus_allee solve_DBG internally ..." << std::endl;
	this->execSolveDBG();
	std::cerr << "phase completed!" << std::endl;
}


void Divide::execSolveDBG(void)
{
	vector<string> singleArgOptionList = {"-o", "-t", "-tmp"};
	vector<string> multiArgOptionList = {"-c", "-s"};

    std::ostringstream oss;
    oss << this->platanusBinary << " solve_DBG -divide_only";

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

    if (system(oss.str().c_str()) != 0) {
        throw platanus::SolveDBGError();
    }
}
