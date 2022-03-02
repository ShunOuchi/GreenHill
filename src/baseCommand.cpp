/*
Copyright (C) 2022 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GreenHill.

GreenHill is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GreenHill is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GreenHill; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "baseCommand.h"
#include <string>
#include <fstream>


//////////////////////////////////////////////////////////////////////////////////////
// check input file format either FASTA or FASTQ
//////////////////////////////////////////////////////////////////////////////////////
platanus::FILETYPE BaseCommand::checkFileFormat(const std::string &filename) const
{
	FILE* fp = platanus::openFileAllowingCompression(filename, "r");

    std::string line[4];
	for (unsigned i = 0; i < 4; ++i)
		platanus::getlineFILE(line[i], fp);

	platanus::closeFileAllowingCompression(fp, filename);

	if (line[0][0]=='>' && strspn(line[1].c_str(), "ACGTN")==line[1].size()) {
		return platanus::FILETYPE::FASTA;
	}
	if (line[0][0]=='@'
		&& strspn(line[1].c_str(), "ACGTN")==line[1].size()
		&& line[2][0]=='+') {

		return platanus::FILETYPE::FASTQ;
	}

    return platanus::FILETYPE::UNKNOWN;
}


//////////////////////////////////////////////////////////////////////////////////////
// parse option argument
//////////////////////////////////////////////////////////////////////////////////////
bool BaseCommand::parseArgs(int argc, char **argv)
{
    int optInd = 1; //edited by ouchi (2->1)
    while (optInd < argc) {
        if (strcmp(argv[optInd], "-h") == 0) {
            return false;
        }
        // parse boolean option (has no other argument)
        if (optionBool.find(argv[optInd]) != optionBool.end()) {
            optionBool[argv[optInd]] = true;
            ++optInd;
        // parse single argument option
        } else if (optionMultiArgs.find(argv[optInd]) != optionMultiArgs.end()) {
            int optIndNext = optInd + 1;
			optionMultiArgs[argv[optInd]].clear();
            while (optIndNext < argc && argv[optIndNext][0] != '-') {
                std::string tmpArgs = argv[optIndNext];
                optionMultiArgs[argv[optInd]].push_back(tmpArgs);
                ++optIndNext;
            }
            optInd = optIndNext;
        // parse pair file (paired-end or mate-pair read file)
        } else if (optionSingleArgs.find(argv[optInd]) != optionSingleArgs.end()) {
            optionSingleArgs[argv[optInd]] = argv[optInd + 1];
            optInd += 2;
        // parse multiple argument option
        } else {
            auto pairFileIterator = pairedEndSingleFileType.begin();
            for (; pairFileIterator != pairedEndSingleFileType.end(); ++pairFileIterator) {
                if (strstr(argv[optInd], (*pairFileIterator).c_str()) == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum == 0) {
                        return false;
                    }
                    int optIndNext = optInd + 1;
                    while (optIndNext < argc && argv[optIndNext][0] != '-') {
                        std::string tmpArgs = argv[optIndNext];
                        PairFile tmpPair(pairNum, *pairFileIterator, tmpArgs);
                        optionPairFile.push_back(tmpPair);
                        ++optIndNext;
                    }
                    optInd = optIndNext;
                    break;
                }
            }
            if (pairFileIterator == pairedEndSingleFileType.end()) {
                pairFileIterator = pairedEndPairFileType.begin();
                for (; pairFileIterator != pairedEndPairFileType.end(); ++pairFileIterator) {
                    if (strstr(argv[optInd], (*pairFileIterator).c_str()) == argv[optInd]) {
                        int pairNum = divideArgvInt(argv[optInd]);
                        if (pairNum == 0) {
                            return false;
                        }
                        int optIndNext = optInd + 1;
                        while (optIndNext < argc && argv[optIndNext][0] != '-') {
                            std::string tmpArgs1 = argv[optIndNext];
                            ++optIndNext;
                            if (optIndNext > argc) return false;
                            if (argv[optIndNext][0] == '-') return false;
                            std::string tmpArgs2 = argv[optIndNext];
                            PairFile tmpPair(pairNum, *pairFileIterator, tmpArgs1, tmpArgs2);
                            optionPairFile.push_back(tmpPair);
                            ++optIndNext;
                        }
                        optInd = optIndNext;
                        break;
                    }
                }
                if (pairFileIterator == pairedEndPairFileType.end()) {
                    int add = checkOtherOption(argv[optInd]);
                    if (add == 0) {
                        return false;
                    } else {
                        optInd += add;
                    }
                }
            }
        }
    }
    return checkFileEnough();
}






//////////////////////////////////////////////////////////////////////////////////////
// recognise option number (ex: -IP1 <- get 1)
//////////////////////////////////////////////////////////////////////////////////////
int BaseCommand::divideArgvInt(char *args) const
{
    char *moveArgs = args;
    for (;moveArgs != '\0'; ++moveArgs) {
        if (isdigit(*moveArgs)) {
            return atoi(moveArgs);
        }
    }
    return 0;
}


