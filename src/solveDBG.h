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

#ifndef SOLVE_DBG_H
#define SOLVE_DBG_H

#include "baseCommand.h"
#include "mapper.h"
#include "scaffold.h"
#include "scaffoldGraph.h"
#include "pairedDBG.h"
#include <cmath>


//class SolveDBG : public BaseCommand
class SolveDBG : public Scaffold
{
private:

	static const double MIN_LONG_READ_LENGTH_CUTOFF_FACTOR;
	static const double MAX_LONG_READ_LENGTH_CUTOFF_FACTOR;
    static const double LONG_READ_LENGTH_CUTOFF_UNIT_FACTOR;
	static const long LONG_READ_MIN_ALIGNMENT_LENGTH;
	static const double LONG_READ_MIN_ALIGNMENT_COVERAGE;
	static const double LONG_READ_MIN_ALIGNMENT_IDENTITY;

    PairedDBG pairedDBG;
    PairedDBG unphasedGraph;
    PairedDBG phasedGraph;
	std::vector<int> multiSeedLengthForShortRead;
	std::vector<int> multiSeedLengthForLongRead;
	long minLinkToPhase;
	long minOverlapForScaffolding;

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            int optInd = 1; //edited by ouchi (2->1)
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
        if (optionMultiArgs["-c"].size() == 0 && optionMultiArgs["-cph"].size() == 0 && optionSingleArgs["-cgfa"] == "") { //edited by ouchi
            std::cerr << "Error: not specified contig file!!" << std::endl;
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
    SolveDBG();
    SolveDBG(const SolveDBG &) = delete;
    SolveDBG &operator=(const SolveDBG &) = delete;
    ~SolveDBG() = default;

    virtual void usage(void) const;
    virtual void exec(void);
    virtual void initializeParameters(void);
    virtual void readLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread);
    void mapLibrary(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread); //added by ouchi
    virtual void mapLibraryAndInitGraph(const int numThread);
    virtual void outputAndAfterTreatment(void);
	virtual void generateGraphFastg();
	void outputGraph(const char *suffix);
	void updateAndWriteInsertSize(const long libraryIndex);
	void extendConsensusToEstimateInsertSize(void);
	void extendConsensus(const long numOuterIteration, const bool divisionFlag);
	void execMinialign(const std::string targetFilename, const std::string &readFilename, const std::string &outFilename, const long numThread, const std::string minialignExecutable);
	void execMinimap(const std::vector<std::string> &contigFilenames, const std::vector<std::string> &bubbleFilenames, const std::vector<std::string> &readFilenames, const std::string &outFilename, const long minAlignmentLength, const long numThread, const std::string minimapExecutable);
	void execMinimap2(const std::vector<std::string> &contigFilenames, const std::vector<std::string> &bubbleFilenames, const std::vector<std::string> &readFilenames, const std::string &outFilename, const long numThread, const std::string minimap2Executable);
	void mincingBubble(std::unique_ptr<HeteroMapper> &mapper, platanus::Contig &contig, const int numThread); //added by ouchi
};




#endif
