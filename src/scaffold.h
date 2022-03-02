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

#ifndef SCAFFOLD_H
#define SCAFFOLD_H

#include "baseCommand.h"
#include "mapper.h"
#include "scaffoldGraph.h"
#include <cmath>


class Scaffold : public BaseCommand
{
protected:

    static const unsigned MIN_SCAFFOLD_LEN;
    static const unsigned MIN_TOL_FACTOR;
    static const unsigned MAX_TOL_FACTOR;
    static const double DEFAULT_INS_CUTOFF_RATE;
	static const double SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR;
	static const double SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR;
	static const double LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR;
	static const double LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR;


    std::unordered_map<int, int> optionMinIns;
    std::unordered_map<int, int> optionAveIns;
    std::unordered_map<int, int> optionSDIns;
    std::unordered_map<int, double> optionInsCutoffRate;
    unsigned numThread;
    unsigned seedLength;
    int keyLength;
    long minLink;
    double averageCoverage;
    double bubbleThreshold;
    unsigned long long contigMaxK;
    unsigned long long contigReadLength;
    std::vector<std::vector<SeqLib> > libraryMT;
    std::vector<int> numFilePerLibraryID;
    std::vector<int> libraryIDList;
	std::vector<SeqLib> tagLibraryMT;
	std::vector<SeqLib> longReadLibraryMT;
	std::vector<SeqLib> HiCLibraryMT;
    ScaffoldGraph scaffoldGraph;

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            int optInd = 2;
            while (optInd < argc) {
                if (strstr(argv[optInd], "-n") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum == 0) {
                        return false;
                    }
                    int minInsertSize = atoi(argv[optInd + 1]);
                    optionMinIns[pairNum] = minInsertSize;
                    optInd += 2;
                } else if (strstr(argv[optInd], "-a") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum == 0) {
                        return false;
                    }
                    int minInsertSize = atoi(argv[optInd + 1]);
                    optionAveIns[pairNum] = minInsertSize;
                    optInd += 2;
                } else if (strstr(argv[optInd], "-d") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum == 0) {
                        return false;
                    }
                    int minInsertSize = atoi(argv[optInd + 1]);
                    optionSDIns[pairNum] = minInsertSize;
                    optInd += 2;
                } else if (strstr(argv[optInd], "-z") == argv[optInd]) {
                    int pairNum = divideArgvInt(argv[optInd]);
                    if (pairNum == 0) {
                        return false;
                    }
                    double insCutoffRate = atof(argv[optInd + 1]);
                    optionInsCutoffRate[pairNum] = insCutoffRate;
                    optInd += 2;
                } else {
                    ++optInd;
                }
            }
        } else return false;
        return true;
    }

    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].size() == 0) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        if (optionPairFile.size() == 0) {
            std::cerr << "Error: not specified pair_read file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        if (strstr(argv, "-n") == argv)
            return 2;
        if (strstr(argv, "-a") == argv)
            return 2;
        if (strstr(argv, "-d") == argv)
            return 2;
        if (strstr(argv, "-z") == argv)
            return 2;
        else
            return 0;
    }


public:
    Scaffold();
    Scaffold(const Scaffold &) = delete;
    Scaffold &operator=(const Scaffold &) = delete;
    ~Scaffold() = default;

    virtual void usage(void) const;
    virtual void exec(void);
    void initializeParameters(void);
    virtual void mapLibraryAndInitGraph(const int numThread);
    virtual void readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, platanus::Contig &bubble, const int numThread);
    virtual void outputAndAfterTreatment(void);
    void printInsertSizeFreq(const std::string &outputFilename, const std::vector<long>& distribution);
	virtual void generateGraphFastg();
};




#endif
