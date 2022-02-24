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

#include "iterate.h"
#include <sys/stat.h>
#include <sys/types.h>


const std::string IterateScaffold::CONTIG_FOOTER = "_contig.fa";
const std::string IterateScaffold::SCAF_FOOTER = "_scaffold.fa";
const std::string IterateScaffold::POLISH_FOOTER = "_polished.fa";
const std::string IterateScaffold::GAP_FOOTER = "_gapClosed.fa";
const std::string IterateScaffold::EX_FOOTER = "_extraContig.fa";
const std::string IterateScaffold::MERGE_FOOTER = "_merged.fa";
const std::string IterateScaffold::ITERATION_FOOTER = "_iterativeAssembly.fa";


//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
IterateScaffold::IterateScaffold()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-i"] = "6";
    optionSingleArgs["-m"] = "1";
    optionSingleArgs["-l"] = "";
    optionSingleArgs["-c"] = "";
    optionSingleArgs["-b"] = "";
    optionSingleArgs["-t"] = "1";
    optionSingleArgs["-u"] = "0";
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionSingleArgs["-tmp"] = ".";
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void IterateScaffold::usage(void) const
{

    std::cerr << "\nUsage: platanus_allee iterate [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file and directory (do not use \"/\", default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : contig (or scaffold) file (fasta format)\n"
              << "    -b FILE1 [FILE2 ...]               : bubble-seqquence file (or scaffold) file (fasta format)\n"
              << "    -i INT                             : number of iterations (default " << optionSingleArgs.at("-i") << ")\n"
              << "    -l INT                             : -l value of \"scaffold\" step\n"
              << "    -u FLOAT                           : maximum difference for bubble crush (identity, default " << optionSingleArgs.at("-u") << ")\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -m INT                             : memory limit for making kmer distribution (GB, >=1, default " << optionSingleArgs.at("-m") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"
              << "Outputs:\n"
              << "    PREFIX_iterativeAssembly.fa\n"
              << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// exec gap close
//////////////////////////////////////////////////////////////////////////////////////
void IterateScaffold::exec(void)
{
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());

    const unsigned long long iterateTimes = atoi(optionSingleArgs["-i"].c_str());

    for (unsigned long long times = 1; times <= iterateTimes; ++times) {
        this->setDirectoryName(times);
        this->createDirectory();
        this->createContig(times);
        this->execScaffoldMode();
        this->execPolishMode();
        if (times == iterateTimes) {
            this->execGapCloseMode("-noextend");
        } else {
            this->execGapCloseMode("");
        }
    }
	this->execFinalPolishMode();
    std::ostringstream oss;
    oss << "mv " <<  optionSingleArgs["-o"] << POLISH_FOOTER
        << " " << optionSingleArgs["-o"] << ITERATION_FOOTER;
    if (system(oss.str().c_str()) != 0) {
        throw platanus::CreateLinkError();
    }
}


void IterateScaffold::createDirectory(void) const
{
    struct stat st;

    if (stat(directoryName.c_str(), &st) != 0) {
        int returnValue = mkdir(directoryName.c_str(), 0755);
        if (returnValue != 0) {
            throw platanus::CreateDirError(directoryName);
        }
    }
}


void IterateScaffold::createContig(const unsigned long long times)
{
    if (times == 1) {
        std::ostringstream oss;
        oss << "cp " << this->optionSingleArgs["-c"]
            << " " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER;
        if (system(oss.str().c_str()) != 0) {
            throw platanus::CreateLinkError();
        }
    } else {
        this->execMergeMode(times);
        std::ostringstream oss;
        oss << "mv -f " << this->directoryName << "/" << optionSingleArgs["-o"] << MERGE_FOOTER
            << " " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER;
        if (system(oss.str().c_str()) != 0) {
            throw platanus::CreateLinkError();
        }
    }
}


void IterateScaffold::execMergeMode(const unsigned long long times)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " merge"
        << " -m " << optionSingleArgs["-m"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -f " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
        << " " << this->previousDirectoryName << "/" << optionSingleArgs["-o"] << EX_FOOTER
        << " -k " << 1.0 + 0.5 * ((times - 1) / 3)
        << " -l " << 1.0 + 0.5 * ((times - 1) / 3)
        << " -u " <<  optionSingleArgs["-u"]
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
        << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".mergeLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::MergeError();
    }
}


void IterateScaffold::execScaffoldMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " scaffold"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << CONTIG_FOOTER
        << " -b " << optionSingleArgs["-b"]
        << " -u " <<  optionSingleArgs["-u"]
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"];
    if (optionSingleArgs["-l"] != "") {
        oss << " -l " << optionSingleArgs["-l"];
    }
    if (optionSingleArgs["-r"] != "") {
        oss << " -r " << optionSingleArgs["-r"];
    }

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".scafLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::ScaffoldError();
    }
}


void IterateScaffold::execPolishMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " polish"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << SCAF_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"];

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".polLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::PolishError();
    }
}


void IterateScaffold::execFinalPolishMode(void)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " polish"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << GAP_FOOTER
        << " -o " << optionSingleArgs["-o"];

    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }

    oss << " 2>" << optionSingleArgs["-o"] << ".finalPolishLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::PolishError();
    }
}


void IterateScaffold::execGapCloseMode(const std::string &mode)
{
    std::ostringstream oss;
    oss << this->platanusBinary << " gap_close"
        << " -t " << optionSingleArgs["-t"]
        << " -tmp " << optionSingleArgs["-tmp"]
        << " -c " << this->directoryName << "/" << optionSingleArgs["-o"] << POLISH_FOOTER
        << " -o " << this->directoryName << "/" << optionSingleArgs["-o"]
        << " " << mode;
    for (auto itr = optionPairFile.begin(), end = optionPairFile.end(); itr != end; ++itr) {
        oss << " " << itr->libraryType << itr->libraryID << " " << itr->fileFirst;
        if (itr->fileSecond != "") {
            oss << " " << itr->fileSecond;
        }
    }
    oss << " 2>" << this->directoryName << "/" << optionSingleArgs["-o"] << ".gapLog";

    if (system(oss.str().c_str()) != 0) {
        throw platanus::GapError();
    }
}
