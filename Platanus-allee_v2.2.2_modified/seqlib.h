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

#ifndef SEQLIB_H
#define SEQLIB_H

#include "common.h"
#include <unordered_map>
#include <cctype>

//////////////////////////////////////////////////////////////////////////////////////
// library file class
//////////////////////////////////////////////////////////////////////////////////////
class SeqLib
{
private:
    static const double INS_DISTR_TRUNC;
    static const double INS_DISTR_TRUNC_SD_RATE;
    static const unsigned INS_DISTR_TRUNC_NUM_ITERATION;
    static const double INS_CUTOFF_RATE_TO_PEAK;
    static const long INS_PEAK_WINDOW;

    // member variable
    long numPair;
    long totalLength;
    long averageLength;
    long averageInsSize;
    long sdInsSize;
    double averageCoverage;
    double insCutoffRate;

    template <typename T>
    void truncateDistribution(std::vector<T>& distribution, const double edge) const;
    template <typename T>
    void truncateDistributionByNumber(std::vector<T>& distribution, const double edge) const;
    template <typename T>
    void truncateDistributionBySD(std::vector<T>& distribution, const double edge);
    template <typename T>
    double calcDistributionAverage(const std::vector<T> &distribution, const long lowerBound, const long upperBound) const;
    template <typename T>
    double calcDistributionSD(const std::vector<T> &distribution, const long lowerBound, const long upperBound) const;
    template <typename T>
    long findDistributionPeak(const std::vector<T> &distribution, const long windowSize, long minPeakThreshold) const;
    void normalizeDistribution(const std::vector<long>& preDist, std::vector<double>& postDist, const std::vector<long>& seqLengths) const;

public:
    FILE *pairFP;
    FILE *insertLengthFP;
    FILE *mappedSameContigFP;
    FILE *mappedFP;
    FILE *mappedReadFP;
    enum FORREV {FORWARD = 0, REVERSE, FWTOTAL};

    SeqLib():numPair(0), totalLength(0), averageLength(0), averageInsSize(0), sdInsSize(0), averageCoverage(0), pairFP(NULL), insertLengthFP(NULL), mappedSameContigFP(NULL), mappedFP(NULL), mappedReadFP(NULL) {}
    void estimateInsSize(const std::vector<long>& distribution, const long minPeakThreshold, const double lowerBoundFactor, const double upperBoundFactor);
    void estimateInsSizeNormalized(const std::vector<long>& insSizeDistribution, std::vector<long>& seqLengths);
    inline void makeTempPairFP(void)
    { pairFP = platanus::makeTemporaryFile(); }


    inline void setAverageInsSize(const long num) { averageInsSize = num; }
    inline void setAverageLength(const long num) { averageLength = num; }
    inline void setSDInsSize(const long num) { sdInsSize = num; }
    inline void setAverageCoverage(const double num) { averageCoverage = num; }
    inline void setInsCutoffRate(const double num) { insCutoffRate = num; }

    inline long getAverageInsSize(void) const { return averageInsSize; }
    inline double getAverageCoverage(void) const { return averageCoverage; }

    inline long getSDInsSize(void) const {return sdInsSize; }
    inline long getAverageLength(void) const {return averageLength; }
    inline long getTotalLength(void) const {return totalLength; }
    inline long getNumPair(void) const {return numPair; }
    void printInsertSizeFreq(const std::string &outputFilename);
    void deleteErrorneousMappedPair(const std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink);

    inline void addTotalLength(const long length) { totalLength += length; }
    inline void addNumPair(const long num) { numPair += num; }

    void readInsertSizeFile(std::vector<long>& distribution);
};

class SeqLibInput
{
private:
    enum PAIRFORMAT {UNDEF = 0, PE1, PE2, MP1, MP2};
    std::string filename[2];
    platanus::FILETYPE seqFileType;
    PAIRFORMAT pairFileFormat;
    long id;
public:
    SeqLibInput():filename(), seqFileType(platanus::FILETYPE::UNKNOWN), pairFileFormat(UNDEF), id(0) {};

};



inline bool operator<(const std::vector<SeqLib> &a, const std::vector<SeqLib> &b)
{
    return a[0].getAverageInsSize() - b[0].getAverageInsSize() < 0;
}

void ReadFastaSingleMT(std::vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate=false, const bool isFastq=false, const bool notPair=false);
void ReadFastaSingleUncompressedMT(std::vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate=false, const bool isFastq=false, const bool notPair=false);
void ReadFastaSingleCompressedMT(std::vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate=false, const bool isFastq=false, const bool notPair=false);
void ReadFastaPairMT(std::vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate=false, const bool isFastq=false);
void ReadFastaPairUncompressedMT(std::vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate=false, const bool isFastq=false);
void ReadFastaPairCompressedMT(std::vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate=false, const bool isFastq=false);
void ReadFastaSingleTaggedMT(std::vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaSingleTaggedUncompressedMT(std::vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaSingleTaggedCompressedMT(std::vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaPairTaggedMT(std::vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaPairTaggedUncompressedMT(std::vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaPairTaggedCompressedMT(std::vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaHeaderUncompressed(std::ifstream &ifs);
void ReadFastqHeaderUncompressed(std::ifstream &ifs);
void ReadFastaHeaderCompressed(FILE *fp);
void ReadFastqHeaderCompressed(FILE *fp);
void ReadFastaSeqUncompressed(platanus::SEQ &seq, std::ifstream &ifs);
void ReadFastqSeqUncompressed(platanus::SEQ &seq, std::ifstream &ifs);
void ReadFastaSeqTaggedUncompressed(platanus::SEQ &seq, int &tagID, std::ifstream &ifs, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastqSeqTaggedUncompressed(platanus::SEQ &seq, int &tagID, std::ifstream &ifs, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastaSeqCompressed(platanus::SEQ &seq, FILE *fp);
void ReadFastqSeqCompressed(platanus::SEQ &seq, FILE *fp);
void ReadFastaSeqTaggedCompressed(platanus::SEQ &seq, int &tagID, FILE *fp, std::unordered_map<std::string, int> &tagStringConverter);
void ReadFastqSeqTaggedCompressed(platanus::SEQ &seq, int &tagID, FILE *fp, std::unordered_map<std::string, int> &tagStringConverter);
void setTagStringConverter(const std::vector<std::string> &filenames, std::unordered_map<std::string, int> &tagStringConverter);
bool tagPositionInline(std::string line, size_t &start, size_t &end);


#endif
