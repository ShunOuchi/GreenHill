/*
Copyright (C) 2022 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus-3D.

Platanus-3D is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus-3D is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus-3D; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "baseCommand.h"
#include "seqlib.h"
#include <vector>
#include <unordered_map>
#include <cmath>

using std::vector;
using std::cerr;
using std::endl;
using std::ifstream;


//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
//const double SeqLib::INS_DISTR_TRUNC = 0.01;
const double SeqLib::INS_DISTR_TRUNC = 0.025;
const double SeqLib::INS_DISTR_TRUNC_SD_RATE = 3.0;
const unsigned SeqLib::INS_DISTR_TRUNC_NUM_ITERATION = 1000;
const double SeqLib::INS_CUTOFF_RATE_TO_PEAK = 0.5;
const long SeqLib::INS_PEAK_WINDOW = 101;



//////////////////////////////////////////////////////////////////////////////////////
// cut outlier value in distribution of insert size
// default cut fraction is 0.01 in both-side
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SeqLib::truncateDistribution(vector<T> &distribution, const double edge) const
{
    if (distribution.size() <= 1)
        return;

    T totalSum = 0;
    for (unsigned long i = 1; i < distribution.size(); ++i)
        totalSum += distribution[i] * i;

    T edgeSum = 0;
    for (unsigned long i = 1; static_cast<double>(edgeSum) < totalSum * edge; ++i) {
        edgeSum += distribution[i] * i;
        distribution[i - 1] = 0;
    }
    edgeSum = (distribution.size() - 1) * distribution[distribution.size() - 1];
    for (unsigned long i = distribution.size() - 2; static_cast<double>(edgeSum) < totalSum * edge; --i) {
        edgeSum += distribution[i] * i;
        distribution[i + 1] = 0;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// cut outlier value in distribution of insert size
// default cut fraction is 0.01 in both-side
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SeqLib::truncateDistributionByNumber(vector<T> &distribution, const double edge) const
{
    if (distribution.size() <= 1)
        return;

    T totalSum = 0;
    for (unsigned long i = 1; i < distribution.size(); ++i)
        totalSum += distribution[i];
    T finalSum = edge / 2.0 * totalSum;

    unsigned long lowerThreshold = 0;
    T edgeSum = 0;
    for (; edgeSum < finalSum; ++lowerThreshold) {
        edgeSum += distribution[lowerThreshold];
        if (edgeSum > finalSum) {
            distribution[lowerThreshold] = edgeSum - finalSum;
            break;
        } else {
            distribution[lowerThreshold] = 0;
        }
    }
    unsigned long upperThreshold = distribution.size() - 1;
    edgeSum = 0;
    for (; edgeSum < finalSum; --upperThreshold) {
        edgeSum += distribution[upperThreshold];
        if (edgeSum + distribution[upperThreshold] > finalSum) {
            distribution[upperThreshold] = edgeSum - finalSum;
            break;
        } else {
            distribution[upperThreshold] = 0;
        }
    }
    cerr << "LOWER_THRESHOLD = " << lowerThreshold << ", UPPER_THRESHOLD = " << upperThreshold << endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// cut outlier value in distribution of insert size
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void SeqLib::truncateDistributionBySD(vector<T> &distribution, const double edge)
{
    if (distribution.size() <= 1)
        return;

    averageInsSize = calcDistributionAverage(distribution, distribution.size() - 1) + 0.5;
    sdInsSize = calcDistributionSD(distribution, distribution.size() - 1) + 0.5;
    long bestAverage = getAverageInsSize();
    long bestSD = getSDInsSize();

    long lowerThreshold;
    long upperThreshold;
    unsigned count = 0;
    while (true) {
        long preAverage = getAverageInsSize();
        long preSD = getSDInsSize();
        vector<T> tmpDistribution = distribution;
        lowerThreshold = preAverage - (edge * preSD - 0.5);
        upperThreshold = preAverage + (edge * preSD + 0.5);
        for (long idx = 0; idx < lowerThreshold; ++idx) {
            tmpDistribution[idx] = 0;
        }
        for (unsigned long idx = upperThreshold; idx < tmpDistribution.size(); ++idx) {
            tmpDistribution[idx] = 0;
        }

        averageInsSize = calcDistributionAverage(distribution, distribution.size() - 1) + 0.5;
        sdInsSize = calcDistributionSD(distribution, distribution.size() - 1) + 0.5;
        long tmpAverage = getAverageInsSize();
        long tmpSD = getSDInsSize();
        if (tmpAverage == preAverage) {
            break;
        }
        if (tmpSD <= preSD) {
            bestAverage = tmpAverage;
            bestSD = tmpSD;
        }
        if (++count > INS_DISTR_TRUNC_NUM_ITERATION) {
	  cerr << "Warning: fail to estimate AVE_INS and SD_INS!" << endl;
	  lowerThreshold = bestAverage - (edge * bestSD + 0.5);
	  upperThreshold = bestAverage + (edge * bestSD + 0.5);
	  break;
        }
    }
    cerr << "LOWER_THRESHOLD = " << lowerThreshold << ", UPPER_THRESHOLD = " << upperThreshold << endl;

    for (long idx = 0; idx < lowerThreshold; ++idx) {
        distribution[idx] = 0;
    }
    for (unsigned long idx = upperThreshold; idx < distribution.size(); ++idx) {
        distribution[idx] = 0;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate average of distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
double SeqLib::calcDistributionAverage(const vector<T> &distribution, const long lowerBound, const long upperBound) const
{
	if (lowerBound >= static_cast<long>(distribution.size()))
		return 0.0;

    T sumDistribution = 0;
    T numDistribution = 0;

    for (long i = lowerBound; i <= upperBound; ++i) {
        sumDistribution += distribution[i] * i;
        numDistribution += distribution[i];
    }

	if (numDistribution > 0)
		return static_cast<double>(sumDistribution) / numDistribution + 0.5;
	else
		return 0.0;
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate standard deviation of distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
double SeqLib::calcDistributionSD(const vector<T> &distribution, const long lowerBound, const long upperBound) const
{
    double numerator = 0;
    T denominator = 0;

    for (long i = lowerBound; i <= upperBound; ++i) {
        numerator += (static_cast<double>(i) - averageInsSize) * (static_cast<double>(i) - averageInsSize) * distribution[i];
        denominator += distribution[i];
    }
    if (denominator > 1)
        return sqrt(numerator / static_cast<double>(denominator - 1));
    else
        return 0.0;
}


void SeqLib::normalizeDistribution(const vector<long>& preDist, vector<double>& postDist, const vector<long>& seqLengths) const
{
    vector<long> probability(preDist.size(), 0);
    for (unsigned long idx = 0; idx < seqLengths.size(); ++idx) {
        unsigned long tmpLength = seqLengths[idx];
        for (unsigned long i = 1, end = std::min(tmpLength + 1, preDist.size()); i < end; ++i) {
            probability[i] += tmpLength + 1 - i;
        }
    }

    postDist.clear();
    postDist.resize(preDist.size(), 0);
    for (unsigned long idx = 1; idx < preDist.size(); ++idx) {
        postDist[idx] = static_cast<double>(preDist[idx]) / probability[idx];
    }

    double sumPreDist  = 0;
    double sumPostDist = 0;
    for (unsigned long idx = 1; idx < preDist.size(); ++idx) {
        sumPreDist  += preDist[idx];
        sumPostDist += postDist[idx];
    }
    double rate = static_cast<double>(sumPreDist) / sumPostDist;
    for (unsigned long idx = 1; idx < postDist.size(); ++idx) {
        postDist[idx] *= rate;
    }

}


//////////////////////////////////////////////////////////////////////////////////////
// estimate library insert size
//////////////////////////////////////////////////////////////////////////////////////
void SeqLib::estimateInsSize(const vector<long>& distribution, const long minPeakThreshold, const double lowerBoundFactor, const double upperBoundFactor)
{
    cerr << "estimating insert-size..." << endl;

	long peakInsSize = findDistributionPeak(distribution, INS_PEAK_WINDOW, minPeakThreshold);
	long upperBound = std::min(static_cast<long>(upperBoundFactor * peakInsSize + 0.5), static_cast<long>(distribution.size() - 1));
	long lowerBound = std::min(static_cast<long>(lowerBoundFactor * peakInsSize + 0.5), upperBound);

    averageInsSize = calcDistributionAverage(distribution, lowerBound, upperBound) + 0.5;

	if (averageInsSize != 0) {
		sdInsSize = calcDistributionSD(distribution, lowerBound, upperBound) + 0.5;
		cerr << "PEAK = " << peakInsSize << endl;
		cerr << "LOWER_LIMIT (permissible range to estimate AVE_INS)= " << lowerBound << endl;
		cerr << "UPPER_LIMIT (permissible range to estimate AVE_INS)= " << upperBound << endl;
		cerr << "AVE_INS = " << averageInsSize << endl;
		cerr << "SD_INS = " << sdInsSize << endl;
	}
	else {
		sdInsSize = 0;
	}
}

//////////////////////////////////////////////////////////////////////////////////////
// estimate library insert size
//////////////////////////////////////////////////////////////////////////////////////
void SeqLib::estimateInsSizeNormalized(const vector<long>& insSizeDistribution, vector<long>& seqLengths)
{
    vector<double> distribution;
    normalizeDistribution(insSizeDistribution, distribution, seqLengths);
    truncateDistributionByNumber(distribution, insCutoffRate);

    //long maxInsSize = std::min(2 * findDistributionPeak(normalizedDistribution, INS_PEAK_WINDOW), static_cast<long>(normalizedDistribution.size() - 1), 0);
    averageInsSize = calcDistributionAverage(distribution, 0, distribution.size() - 1) + 0.5;
    sdInsSize = calcDistributionSD(distribution, 0, distribution.size() - 1) + 0.5;
}


void SeqLib::readInsertSizeFile(vector<long>& distribution)
{
    distribution.clear();
    long insLen;
    rewind(insertLengthFP);
    while (fread(&insLen, sizeof(long), 1, insertLengthFP)) {
        if (insLen >= static_cast<long>(distribution.size())) {
            distribution.resize(insLen + 1, static_cast<long>(0));
        }
        ++distribution[insLen];
    }
    if (distribution.size() == 0) {
        throw platanus::MapError("No read mapped in the same contig!!");
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// Output insert size frequency file
//////////////////////////////////////////////////////////////////////////////////////
void SeqLib::printInsertSizeFreq(const std::string &outputFilename)
{
    vector<long> distribution;
    readInsertSizeFile(distribution);

    std::ofstream out(outputFilename.c_str());
//    long maxInsSize = std::min(2 * findDistributionPeak(distribution, INS_PEAK_WINDOW), static_cast<long>(distribution.size() - 1));
    long maxInsSize = distribution.size() - 1;
    long i = 1;
    for (; i <= maxInsSize; ++i) {
        out << i << "\t" << distribution[i] << std::endl;
    }
    unsigned long long tmp = 0;
    for (; static_cast<unsigned long long>(i) < distribution.size(); ++i) {
        tmp += distribution[i];
    }
    out << maxInsSize + 1 << "\t" << tmp << std::endl;

    out.close();
}


//////////////////////////////////////////////////////////////////////////////////////
// find distribution peak value
//////////////////////////////////////////////////////////////////////////////////////
template <typename T>
long SeqLib::findDistributionPeak(const vector<T> &distribution, const long windowSize, long minPeakThreshold) const
{
    if (distribution.size() <= static_cast<unsigned long long>(std::min(windowSize, minPeakThreshold))) {
        return distribution.size() / 2;
    }

	minPeakThreshold = std::max(minPeakThreshold, windowSize / 2);

    T pre = 0;
    for (long i = minPeakThreshold - windowSize / 2; i < windowSize; ++i) {
        pre += distribution[i];
    }

    T peak = pre;
    long peakI = minPeakThreshold;

    for (long i = minPeakThreshold - windowSize / 2 + 1, n = distribution.size() - windowSize; i <= n; ++i) {
        T cur = pre - distribution[i - 1] + distribution[i + windowSize - 1];
        if (cur > peak) {
            peak = cur;
            peakI = i + windowSize / 2;
        }
        pre = cur;
    }

    return peakI;
}


//////////////////////////////////////////////////////////////////////////////////////
// Read single FASTA(Q) file and write tempfile
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaSingleMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair)
{
	platanus::FILECOMPRESSION format = platanus::checkFileCompression(filename);

	if (format == platanus::FILECOMPRESSION::UNCOMPRESSED)
		ReadFastaSingleUncompressedMT(lib, filename, numThread, isMate, isFastq, notPair);
	else
		ReadFastaSingleCompressedMT(lib, filename, numThread, isMate, isFastq, notPair);
}


void ReadFastaSingleUncompressedMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair)
{
    long i;
    std::ifstream forwardFile(filename.c_str());
    if (!forwardFile) throw platanus::FILEError(filename);

    platanus::SEQ forwardSeq, reverseSeq;
    void (* const ReadFunction)(platanus::SEQ &seq, std::ifstream &ifs) = isFastq ? ReadFastqSeqUncompressed : ReadFastaSeqUncompressed;
    void (* const ReadHeaderFunction)(std::ifstream &ifs) = isFastq ? ReadFastqHeaderUncompressed : ReadFastaHeaderUncompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    if (!forwardFile.good()) {
        throw platanus::IOError();
    }
    ReadHeaderFunction(forwardFile);
    i = 0;
    while (forwardFile && !forwardFile.eof()) {
        ReadFunction(forwardSeq, forwardFile);
        if (!notPair && forwardFile.eof()) {
            throw platanus::FormatError("the number of read is odd in file.");
        }
        ReadFunction(reverseSeq, forwardFile);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }
        forwardSeq.writeTemporaryFile(lib[i].pairFP);
        reverseSeq.writeTemporaryFile(lib[i].pairFP);
        i = (i + 1) % numThread;
        lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
        lib[0].addNumPair(1);
    }

    forwardFile.close();
}


void ReadFastaSingleCompressedMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair)
{
    long i;
	FILE* fp = platanus::openFileAllowingCompression(filename, "r");

    platanus::SEQ forwardSeq, reverseSeq;
    void (* const ReadFunction)(platanus::SEQ &seq, FILE *fp) = isFastq ? ReadFastqSeqCompressed : ReadFastaSeqCompressed;
    void (* const ReadHeaderFunction)(FILE *fp) = isFastq ? ReadFastqHeaderCompressed : ReadFastaHeaderCompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    ReadHeaderFunction(fp);
    i = 0;
    while (!feof(fp)) {
        ReadFunction(forwardSeq, fp);
        if (!notPair && feof(fp)) {
            throw platanus::FormatError("the number of read is odd in file.");
        }
        ReadFunction(reverseSeq, fp);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }
        forwardSeq.writeTemporaryFile(lib[i].pairFP);
        reverseSeq.writeTemporaryFile(lib[i].pairFP);
        i = (i + 1) % numThread;
        lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
        lib[0].addNumPair(1);
    }

	platanus::closeFileAllowingCompression(fp, filename);
}


void ReadFastaSingleTaggedMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair, std::unordered_map<std::string, int> &tagStringConverter)
{
	platanus::FILECOMPRESSION format = platanus::checkFileCompression(filename);

	if (format == platanus::FILECOMPRESSION::UNCOMPRESSED)
		ReadFastaSingleTaggedUncompressedMT(lib, filename, numThread, isMate, isFastq, notPair, tagStringConverter);
	else
		ReadFastaSingleTaggedCompressedMT(lib, filename, numThread, isMate, isFastq, notPair, tagStringConverter);
}


void ReadFastaSingleTaggedUncompressedMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair, std::unordered_map<std::string, int> &tagStringConverter)
{
    long i;
    std::ifstream forwardFile(filename.c_str());
    if (!forwardFile) throw platanus::FILEError(filename);

    platanus::SEQ forwardSeq, reverseSeq;
	int forwardTagID, reverseTagID;
    void (* const ReadFunction)(platanus::SEQ &seq, int &tagID, std::ifstream &ifs, std::unordered_map<std::string, int> &tagStringConverter) = isFastq ? ReadFastqSeqTaggedUncompressed : ReadFastaSeqTaggedUncompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    if (!forwardFile.good()) {
        throw platanus::IOError();
    }

    i = 0;
    while (forwardFile && !forwardFile.eof()) {
        ReadFunction(forwardSeq, forwardTagID, forwardFile, tagStringConverter);
        if (!notPair && forwardFile.eof()) {
            throw platanus::FormatError("the number of read is odd in file.");
        }
        ReadFunction(reverseSeq, reverseTagID, forwardFile, tagStringConverter);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }

		if (forwardTagID != -1 && reverseTagID != -1) {
			forwardSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&forwardTagID, sizeof(int), 1, lib[i].pairFP);

			reverseSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&reverseTagID, sizeof(int), 1, lib[i].pairFP);

			i = (i + 1) % numThread;
			lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
			lib[0].addNumPair(1);
		}
    }

    forwardFile.close();
}


void ReadFastaSingleTaggedCompressedMT(vector<SeqLib> &lib, const std::string &filename, const int numThread, const bool isMate, const bool isFastq, const bool notPair, std::unordered_map<std::string, int> &tagStringConverter)
{
    long i;
	FILE* fp = platanus::openFileAllowingCompression(filename, "r");

    platanus::SEQ forwardSeq, reverseSeq;
	int forwardTagID, reverseTagID;
    void (* const ReadFunction)(platanus::SEQ &seq, int &tagID, FILE *fp, std::unordered_map<std::string, int> &tagStringConverter) = isFastq ? ReadFastqSeqTaggedCompressed : ReadFastaSeqTaggedCompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);

    i = 0;
    while (!feof(fp)) {
        ReadFunction(forwardSeq, forwardTagID, fp, tagStringConverter);
        if (!notPair && feof(fp)) {
            throw platanus::FormatError("the number of read is odd in file.");
        }
        ReadFunction(reverseSeq, reverseTagID, fp, tagStringConverter);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }

		if (forwardTagID != -1 && reverseTagID != -1) {
			forwardSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&forwardTagID, sizeof(int), 1, lib[i].pairFP);

			reverseSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&reverseTagID, sizeof(int), 1, lib[i].pairFP);

			i = (i + 1) % numThread;
			lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
			lib[0].addNumPair(1);
		}
    }

	platanus::closeFileAllowingCompression(fp, filename);
}


//////////////////////////////////////////////////////////////////////////////////////
// Read pair FASTA(Q) file and write tempfile
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaPairMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq)
{
	platanus::FILECOMPRESSION format1 = platanus::checkFileCompression(filename1);
	platanus::FILECOMPRESSION format2 = platanus::checkFileCompression(filename2);

	if (format1 == platanus::FILECOMPRESSION::UNCOMPRESSED && format2 == platanus::FILECOMPRESSION::UNCOMPRESSED)
		ReadFastaPairUncompressedMT(lib, filename1, filename2, numThread, isMate, isFastq);
	else
		ReadFastaPairCompressedMT(lib, filename1, filename2, numThread, isMate, isFastq);
}


void ReadFastaPairUncompressedMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq)
{
    long i;
    ifstream forwardFile(filename1.c_str());
    if (!forwardFile) throw platanus::FILEError(filename1);
    ifstream reverseFile(filename2.c_str());
    if (!reverseFile) throw platanus::FILEError(filename2);
    platanus::SEQ forwardSeq, reverseSeq;
    void (*ReadFunction)(platanus::SEQ &seq, ifstream &ifs) = isFastq ? ReadFastqSeqUncompressed : ReadFastaSeqUncompressed;
    void (*ReadHeaderFunction)(ifstream &ifs) = isFastq ? ReadFastqHeaderUncompressed : ReadFastaHeaderUncompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    if (!forwardFile.good() || !reverseFile.good()) {
        throw platanus::IOError();
    }
    ReadHeaderFunction(forwardFile);
    ReadHeaderFunction(reverseFile);

    i = 0;

    while (forwardFile && reverseFile && !forwardFile.eof() && !reverseFile.eof()) {
        ReadFunction(forwardSeq, forwardFile);
        ReadFunction(reverseSeq, reverseFile);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }
        forwardSeq.writeTemporaryFile(lib[i].pairFP);
        reverseSeq.writeTemporaryFile(lib[i].pairFP);
        i = (i + 1) % numThread;

        lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
        lib[0].addNumPair(1);
    }

    if (!forwardFile.eof() || !reverseFile.eof()) {
        throw platanus::FormatError("the number of read is different in paired-file.");
    }

    forwardFile.close();
    reverseFile.close();

}


void ReadFastaPairCompressedMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq)
{
    long i;
	FILE* forwardFP = platanus::openFileAllowingCompression(filename1, "r");
	FILE* reverseFP = platanus::openFileAllowingCompression(filename2, "r");
    platanus::SEQ forwardSeq, reverseSeq;
    void (*ReadFunction)(platanus::SEQ &seq, FILE *fp) = isFastq ? ReadFastqSeqCompressed : ReadFastaSeqCompressed;
    void (*ReadHeaderFunction)(FILE *fp) = isFastq ? ReadFastqHeaderCompressed : ReadFastaHeaderCompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    ReadHeaderFunction(forwardFP);
    ReadHeaderFunction(reverseFP);

    i = 0;
    while (!feof(forwardFP) && !feof(reverseFP)) {
        ReadFunction(forwardSeq, forwardFP);
        ReadFunction(reverseSeq, reverseFP);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }
        forwardSeq.writeTemporaryFile(lib[i].pairFP);
        reverseSeq.writeTemporaryFile(lib[i].pairFP);
        i = (i + 1) % numThread;

        lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
        lib[0].addNumPair(1);
    }

    if (!feof(forwardFP) || !feof(reverseFP)) {
        throw platanus::FormatError("the number of read is different in paired-file.");
    }

	platanus::closeFileAllowingCompression(forwardFP, filename1);
	platanus::closeFileAllowingCompression(reverseFP, filename2);
}


void ReadFastaPairTaggedMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq, std::unordered_map<std::string, int> &tagStringConverter)
{
	platanus::FILECOMPRESSION format1 = platanus::checkFileCompression(filename1);
	platanus::FILECOMPRESSION format2 = platanus::checkFileCompression(filename2);

	if (format1 == platanus::FILECOMPRESSION::UNCOMPRESSED && format2 == platanus::FILECOMPRESSION::UNCOMPRESSED)
		ReadFastaPairTaggedUncompressedMT(lib, filename1, filename2, numThread, isMate, isFastq, tagStringConverter);
	else
		ReadFastaPairTaggedCompressedMT(lib, filename1, filename2, numThread, isMate, isFastq, tagStringConverter);
}


void ReadFastaPairTaggedUncompressedMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq, std::unordered_map<std::string, int> &tagStringConverter)
{
    long i;
    ifstream forwardFile(filename1.c_str());
    if (!forwardFile) throw platanus::FILEError(filename1);
    ifstream reverseFile(filename2.c_str());
    if (!reverseFile) throw platanus::FILEError(filename2);
    platanus::SEQ forwardSeq, reverseSeq;
	int forwardTagID, reverseTagID;
    void (*ReadFunction)(platanus::SEQ &seq, int &tagID, ifstream &ifs, std::unordered_map<std::string, int> &tagStringConverter) = isFastq ? ReadFastqSeqTaggedUncompressed : ReadFastaSeqTaggedUncompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);
    if (!forwardFile.good() || !reverseFile.good()) {
        throw platanus::IOError();
    }

    i = 0;
    while (forwardFile && reverseFile && !forwardFile.eof() && !reverseFile.eof()) {
        ReadFunction(forwardSeq, forwardTagID, forwardFile, tagStringConverter);
        ReadFunction(reverseSeq, reverseTagID, reverseFile, tagStringConverter);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }

		if (forwardTagID == reverseTagID && forwardTagID != -1 && reverseTagID != -1) {
			forwardSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&forwardTagID, sizeof(int), 1, lib[i].pairFP);

			reverseSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&reverseTagID, sizeof(int), 1, lib[i].pairFP);

			i = (i + 1) % numThread;
			lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
			lib[0].addNumPair(1);
		}

    }
    if (!forwardFile.eof() || !reverseFile.eof()) {
        throw platanus::FormatError("the number of read is different in paired-file.");
    }

    forwardFile.close();
    reverseFile.close();
}


void ReadFastaPairTaggedCompressedMT(vector<SeqLib> &lib, const std::string &filename1, const std::string &filename2, const int numThread, const bool isMate, const bool isFastq, std::unordered_map<std::string, int> &tagStringConverter)
{
    long i;
	FILE* forwardFP = platanus::openFileAllowingCompression(filename1, "r");
	FILE* reverseFP = platanus::openFileAllowingCompression(filename2, "r");
    platanus::SEQ forwardSeq, reverseSeq;
	int forwardTagID, reverseTagID;
    void (*ReadFunction)(platanus::SEQ &seq, int &tagID, FILE *fp, std::unordered_map<std::string, int> &tagStringConverter) = isFastq ? ReadFastqSeqTaggedCompressed : ReadFastaSeqTaggedCompressed;
    for (i = 0; i < numThread; ++i)
        fseek(lib[i].pairFP, 0, SEEK_END);

    i = 0;
    while (!feof(forwardFP) && !feof(reverseFP)) {
        ReadFunction(forwardSeq, forwardTagID, forwardFP, tagStringConverter);
        ReadFunction(reverseSeq, reverseTagID, reverseFP, tagStringConverter);
        if (isMate) {
            forwardSeq.reverse();
            reverseSeq.reverse();
        }

		if (forwardTagID == reverseTagID && forwardTagID != -1 && reverseTagID != -1) {
			forwardSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&forwardTagID, sizeof(int), 1, lib[i].pairFP);

			reverseSeq.writeTemporaryFile(lib[i].pairFP);
            fwrite(&reverseTagID, sizeof(int), 1, lib[i].pairFP);

			i = (i + 1) % numThread;
			lib[0].addTotalLength(forwardSeq.length + reverseSeq.length);
			lib[0].addNumPair(1);
		}

    }
    if (!feof(forwardFP) || !feof(reverseFP)) {
        throw platanus::FormatError("the number of read is different in paired-file.");
    }

	platanus::closeFileAllowingCompression(forwardFP, filename1);
	platanus::closeFileAllowingCompression(reverseFP, filename2);
}


//////////////////////////////////////////////////////////////////////////////////////
// Jump fasta head header
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaHeaderUncompressed(ifstream &ifs)
{
    std::string line;
    while (ifs && getline(ifs, line)) {
        if (line[0] == '>') {
            break;
        }
    }
}


void ReadFastaHeaderCompressed(FILE *fp)
{
    std::string line;
    while (platanus::getlineFILE(line, fp) != NULL) {
        if (line[0] == '>') {
            break;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// Jump fastq head header
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastqHeaderUncompressed(ifstream &ifs)
{
    std::string line;
    while (ifs && getline(ifs, line)) {
        if (line[0] == '@') {
            break;
        }
    }
}


void ReadFastqHeaderCompressed(FILE *fp)
{
    std::string line;
    while (platanus::getlineFILE(line, fp) != NULL) {
        if (line[0] == '@') {
            break;
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// Read one fasta seq
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastaSeqUncompressed(platanus::SEQ &seq, std::ifstream &ifs)
{
    std::string line;
    seq.length = 0;
    seq.base = "";
    while (ifs && getline(ifs, line)) {
        if (line[0] == '>') {
            break;
        } else {
            seq.put(line);
        }
    }
}


void ReadFastaSeqTaggedUncompressed(platanus::SEQ &seq, int &tagID, std::ifstream &ifs, std::unordered_map<std::string, int> &tagStringConverter)
{
    std::string line;

    seq.length = 0;
    seq.base = "";

	size_t tagInfoStart = 0, tagInfoEnd = 0;
    getline(ifs, line);
	if (tagPositionInline(line, tagInfoStart, tagInfoEnd))
		tagID = tagStringConverter[line.substr(tagInfoStart, tagInfoEnd - tagInfoStart)];
	else
		tagID = -1;

    while (ifs && getline(ifs, line)) {
        if (line[0] == '>') {
            break;
        } else {
            seq.put(line);
        }
    }
}


void ReadFastaSeqCompressed(platanus::SEQ &seq, FILE *fp)
{
    std::string line;
    seq.length = 0;
    seq.base = "";
    while (platanus::getlineFILE(line, fp) != NULL) {
        if (line[0] == '>') {
            break;
        } else {
            seq.put(line);
        }
    }
}


void ReadFastaSeqTaggedCompressed(platanus::SEQ &seq, int &tagID, FILE *fp, std::unordered_map<std::string, int> &tagStringConverter)
{
    std::string line;

    seq.length = 0;
    seq.base = "";

	size_t tagInfoStart = 0, tagInfoEnd = 0;
    platanus::getlineFILE(line, fp);
	if (tagPositionInline(line, tagInfoStart, tagInfoEnd))
		tagID = tagStringConverter[line.substr(tagInfoStart, tagInfoEnd - tagInfoStart)];
	else
		tagID = -1;

    while (platanus::getlineFILE(line, fp) != NULL) {
        if (line[0] == '>') {
            break;
        } else {
            seq.put(line);
        }
    }
}


bool tagPositionInline(std::string line, size_t &start, size_t &end)
{
	static const std::string TAG_START_STRING = "BX:Z:";

	start = line.find(TAG_START_STRING, 0);
	if (start == std::string::npos)
		return false;

	start += TAG_START_STRING.size();
	for (end = start + 1; end < line.size(); ++end) {
		if (!isalnum(line[end]))
			return true;
	}
	return true;
}


//////////////////////////////////////////////////////////////////////////////////////
// Read one fastq seq
//////////////////////////////////////////////////////////////////////////////////////
void ReadFastqSeqUncompressed(platanus::SEQ &seq, std::ifstream &ifs)
{
    unsigned readLine = 0;
    std::string line;
    seq.length = 0;
    seq.base = "";
    while (ifs && getline(ifs, line)) {
        if (line[0] == '+') {
            break;
        } else {
            ++readLine;
            seq.put(line);
        }
    }
    for (unsigned i = 0; i < readLine + 1; ++i)
        getline(ifs, line);
}


void ReadFastqSeqTaggedUncompressed(platanus::SEQ &seq, int &tagID, std::ifstream &ifs, std::unordered_map<std::string, int> &tagStringConverter)
{
    unsigned readLine = 0;
    std::string line;

    seq.length = 0;
    seq.base = "";

	size_t tagInfoStart = 0, tagInfoEnd = 0;
    getline(ifs, line);
	if (tagPositionInline(line, tagInfoStart, tagInfoEnd))
		tagID = tagStringConverter[line.substr(tagInfoStart, tagInfoEnd - tagInfoStart)];
	else
		tagID = -1;

    while (ifs && getline(ifs, line)) {
        if (line[0] == '+') {
            break;
        } else {
            ++readLine;
            seq.put(line);
        }
    }
    for (unsigned i = 0; i < readLine; ++i)
        getline(ifs, line);
}


void ReadFastqSeqCompressed(platanus::SEQ &seq, FILE *fp)
{
    unsigned readLine = 0;
    std::string line;
    seq.length = 0;
    seq.base = "";
    while (platanus::getlineFILE(line, fp) != NULL) {
        if (line[0] == '+') {
            break;
        } else {
            ++readLine;
            seq.put(line);
        }
    }
    for (unsigned i = 0; i < readLine + 1; ++i)
		platanus::getlineFILE(line, fp);
}


void ReadFastqSeqTaggedCompressed(platanus::SEQ &seq, int &tagID, FILE *fp, std::unordered_map<std::string, int> &tagStringConverter)
{
    unsigned readLine = 0;
    std::string line;

    seq.length = 0;
    seq.base = "";

	size_t tagInfoStart = 0, tagInfoEnd = 0;
    platanus::getlineFILE(line, fp);
	if (tagPositionInline(line, tagInfoStart, tagInfoEnd))
		tagID = tagStringConverter[line.substr(tagInfoStart, tagInfoEnd - tagInfoStart)];
	else
		tagID = -1;

    while (platanus::getlineFILE(line, fp) != NULL) {
        if (line[0] == '+') {
            break;
        } else {
            ++readLine;
            seq.put(line);
        }
    }
    for (unsigned i = 0; i < readLine; ++i)
		platanus::getlineFILE(line, fp);
}


void setTagStringConverter(const vector<std::string> &filenames, std::unordered_map<std::string, int> &tagStringConverter)
{
	tagStringConverter.clear();

	for (auto fileItr = filenames.begin(); fileItr != filenames.end(); ++fileItr) {
		std::string line;
		char headerChar = '\0';

		if (platanus::checkFileCompression(*fileItr) == platanus::FILECOMPRESSION::UNCOMPRESSED) {
			ifstream ifs(*fileItr);
			if (!ifs)
				continue;

			getline(ifs, line);
			if (line.size() > 0)
				headerChar = line[0];
			if (headerChar != '@' && headerChar != '>')
				continue;

			ifs.close();

			ifs.open(*fileItr);
			if (!ifs)
				continue;

			size_t tagInfoStart = 0, tagInfoEnd = 0;
			while (ifs && getline(ifs, line)) {
				if (line.size() > 0 && line[0] == headerChar && tagPositionInline(line, tagInfoStart, tagInfoEnd))
					tagStringConverter[line.substr(tagInfoStart, tagInfoEnd - tagInfoStart)] = 1;
			}
			ifs.close();
		}
		else {
			FILE* fp = platanus::openFileAllowingCompression(*fileItr, "r");

			platanus::getlineFILE(line, fp);
			if (line.size() > 0)
				headerChar = line[0];
			platanus::closeFileAllowingCompression(fp, *fileItr);
			if (headerChar != '@' && headerChar != '>')
				continue;

			fp = platanus::openFileAllowingCompression(*fileItr, "r");

			size_t tagInfoStart = 0, tagInfoEnd = 0;
			while (platanus::getlineFILE(line, fp) != NULL) {
				if (line.size() > 0 && line[0] == headerChar && tagPositionInline(line, tagInfoStart, tagInfoEnd))
					tagStringConverter[line.substr(tagInfoStart, tagInfoEnd - tagInfoStart)] = 1;
			}
			platanus::closeFileAllowingCompression(fp, *fileItr);
		}
    }

	vector<std::string> stringBuffer;
	for (auto itr = tagStringConverter.begin(); itr != tagStringConverter.end(); ++itr)
		stringBuffer.push_back(itr->first);
	std::sort(stringBuffer.begin(), stringBuffer.end());

	int tagID = 0;
	for (auto itr = stringBuffer.begin(); itr != stringBuffer.end(); ++itr) {
		++tagID;
		tagStringConverter[*itr] = tagID;
	}
}


void SeqLib::deleteErrorneousMappedPair(const std::unordered_map<std::pair<int, int>, bool, platanus::PairHash, platanus::PairEqual> &errorLink)
{
    platanus::Position forwardResult;
    platanus::Position reverseResult;
    FILE *newMappedFP = platanus::makeTemporaryFile();

    rewind(mappedFP);
    while (fread(&forwardResult, sizeof(platanus::Position), 1, mappedFP)) {
        fread(&reverseResult, sizeof(platanus::Position), 1, mappedFP);

        std::pair<int, int> key;
        if (abs(forwardResult.id) < abs(reverseResult.id)) {
            key.first = forwardResult.id;
            key.second = -(reverseResult.id);
        }
        else {
            key.first = reverseResult.id;
            key.second = -(forwardResult.id);
        }

        if (errorLink.find(key) == errorLink.end()) {
            fwrite(&forwardResult, sizeof(platanus::Position), 1, newMappedFP);
            fwrite(&reverseResult, sizeof(platanus::Position), 1, newMappedFP);
        }
    }

    fclose(mappedFP);
    mappedFP = newMappedFP;
}
