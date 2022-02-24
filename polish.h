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

#ifndef POLISH_H
#define POLISH_H

#include "baseCommand.h"
#include "common.h"
#include "seqlib.h"
#include "mapper.h"
#include <vector>
#include <string>
#include <unordered_map>
#include <array>


//////////////////////////////////////////////////////////////////////////////////////
// polish execute class
//////////////////////////////////////////////////////////////////////////////////////
class Polish : public BaseCommand
{
private:

    static const unsigned long MIN_COVERAGE;
	
	struct pileupRecord
	{
//		std::array<unsigned, 4> numBase;
		unsigned numRead;
		unsigned numOddRead;

		pileupRecord(): numRead(0), numOddRead(0) {}
        ~pileupRecord() = default;
	};

    platanus::Contig contig;
    unsigned numThread;
    unsigned seedLength;
    int keyLength;
    double averageCoverage;
    unsigned long long contigReadLength;
    unsigned long long contigMaxK;
    std::vector<std::vector<SeqLib> > libraryMT;
    std::vector<int> numFilePerLibraryID;
    std::vector<int> libraryIDList;
    std::vector<std::vector<pileupRecord> > pileup;
	std::vector<std::string> contigFileName;
	std::vector<long> numContigInFile;

    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].size() == 0) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        if (optionPairFile.size() == 0) {
            std::cerr << "Error: not specified read file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }

	void pileupMappedReads(const std::vector<SeqLib> &library, const double minIdentity, const long numThread);
	void pileupSingleRead(const platanus::SEQ &read, const platanus::Position &position, std::vector<std::vector<pileupRecord> > &result);
	void pileupSingleOddRead(const platanus::SEQ &read, const platanus::Position &position, std::vector<std::vector<pileupRecord> > &result);
	void maskErrorBases(const double minNumOddReadRatio);
	void maskLowCoverageEdge(const unsigned long minCoverage);
	void maskShortContig(platanus::SEQ &seq, const long minContigLength);
	void trimEdgeN(platanus::SEQ &seq);
	void printSeq(const std::string &outPrefix, const unsigned long long readLength, const unsigned long long contigMaxK);
	void setFastaFileInfo(const std::vector<std::string> fastaFileName);


public:
    Polish();
    Polish(const Polish &) = delete;
    Polish &operator=(const Polish &) = delete;
    ~Polish() = default;

    virtual void usage(void) const;
    virtual void exec(void);
    void initializeParameters(void);
    void readLibrary(Mapper &mapper, platanus::Contig &contig, std::vector<std::vector<SeqLib> > &libraryMT, std::vector<int> &numFilePerLibraryID, const int numLibrary, const int numThread);
};


#endif
