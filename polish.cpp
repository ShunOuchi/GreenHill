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

#include "polish.h"
#include <utility>
#include <algorithm>
#include <sstream>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const unsigned long Polish::MIN_COVERAGE = 1;

//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
Polish::Polish()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-s"] = "0";
    optionSingleArgs["-e"] = "0.98";
    optionSingleArgs["-r"] = "1.0";
    optionSingleArgs["-l"] = "0";
    optionSingleArgs["-t"] = "1";
    optionMultiArgs["-c"] = vector<string>();
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionSingleArgs["-tmp"] = ".";
}


//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void Polish::usage(void) const
{

    std::cerr << "\nUsage: platanus_allee polish [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : scaffold_file (fasta format)\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -s INT                             : mapping seed length (default " << optionSingleArgs.at("-s") << " (auto))\n"
              << "    -e FLOAT                           : minimum sequence identity for mapping (default " << optionSingleArgs.at("-e") << ")\n"
              << "    -r FLOAT                           : minimum low-identity reads ratio (#low-identity-reads / #mapped-reads) (default " << optionSingleArgs.at("-r") << ")\n"
              << "    -l INT                             : minimum contig length for output (default " << optionSingleArgs.at("-l") << " (auto))\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Output:\n"
              << "    PREFIX_polished.fa\n"
              << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// exec polish
//////////////////////////////////////////////////////////////////////////////////////
void Polish::exec(void)
{
    unsigned i;
	const double minIdentity = atof(optionSingleArgs["-e"].c_str());
	const double minNumOddReadRatio = atof(optionSingleArgs["-r"].c_str());
    int numLibrary = 0;
    int numThread = atoi(optionSingleArgs["-t"].c_str());
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
    vector<int> numFilePerLibraryID(optionPairFile.size());
    vector<int> libraryIDList(optionPairFile.size());
    vector<vector<SeqLib> > libraryMT;
    unsigned long long readLength = platanus::Contig::getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
    if (readLength == 0)
        readLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;


	contigMaxK = contig.getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
	const long minContigLength = atoi(optionSingleArgs["-l"].c_str()) == 0 ? contigMaxK : atoi(optionSingleArgs["-l"].c_str());
	const unsigned int seedLength = atoi(optionSingleArgs["-s"].c_str()) == 0 ? contigMaxK : atoi(optionSingleArgs["-s"].c_str());
    const int keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);

    // nest {} for calling mapper destructor
    {
        Mapper mapper(seedLength, keyLength);
        sort(optionPairFile.begin(), optionPairFile.end());

        // count up the number of file in ecah library
        for (i = 0; i < optionPairFile.size(); ++i) {
            ++(numFilePerLibraryID[numLibrary]);
            libraryIDList[numLibrary] = optionPairFile[i].libraryID;
            if (optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
                ++numLibrary;
            }
        }

        libraryMT.resize(numLibrary);
        omp_set_num_threads(numThread);

        readLibrary(mapper, contig, libraryMT, numFilePerLibraryID, numLibrary, numThread);


        for (int i = 0; i < numLibrary; ++i) {
            std::cerr << "[PAIR_LIBRARY " << libraryIDList[i] << "]" << std::endl;
            libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));

            mapper.mapPairToCalculateCoverage(libraryMT[i], numThread);

            vector<long> insSizeDistribution;
            libraryMT[i][0].readInsertSizeFile(insSizeDistribution);
            libraryMT[i][0].estimateInsSize(insSizeDistribution, i > 0 ? libraryMT[i - 1][0].getAverageInsSize() : 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);
        }

        mapper.releaseContig(contig);
    }

    for (int i = 0; i < numLibrary; ++i) {
        std::cerr << "[PAIR_LIBRARY " << i + 1 << "]" << std::endl;
        pileupMappedReads(libraryMT[i], minIdentity, numThread);

        for (long j = 0; j < numThread; ++j) {
            fclose(libraryMT[i][j].mappedFP);
            libraryMT[i][j].mappedFP = NULL;
        }
        libraryMT[i].clear();
	}

	maskErrorBases(minNumOddReadRatio);
//	maskLowCoverageEdge(MIN_COVERAGE);

	for (unsigned long i = 0; i < pileup.size(); ++i) {
		maskShortContig(contig.seq[i], minContigLength);
		trimEdgeN(contig.seq[i]);
	}

	setFastaFileInfo(optionMultiArgs["-c"]);
	printSeq(optionSingleArgs["-o"], contigReadLength, contigMaxK);

    std::cerr << "polish completed!!" << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// Read library files
// Change read function which fasta or fastq, single file or pair one.
//////////////////////////////////////////////////////////////////////////////////////
void Polish::readLibrary(Mapper &mapper, platanus::Contig &contig, vector<vector<SeqLib> > &libraryMT, vector<int> &numFilePerLibraryID, const int numLibrary, const int numThread)
{
    unsigned nowFileNumber = 0;
    int j;
    omp_set_num_threads(numThread);
    #    pragma omp parallel for private(j, nowFileNumber) schedule(static, 1)
    for (int i = -1; i < numLibrary; ++i) {
        try {
            if (i == -1) {
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i) {
				   contig.readFastaCoverage(optionMultiArgs["-c"][i]);
                }

                mapper.setContig(contig);
                mapper.makeKmerTable();
            } else {
                nowFileNumber = 0;
                libraryMT[i].resize(numThread);
                for (j = 0; j < i; ++j)
                    nowFileNumber += numFilePerLibraryID[j];

                for (j = 0; j < numThread; ++j) {
                    libraryMT[i][j].makeTempPairFP();
                }
                for (j = 0; j < numFilePerLibraryID[i]; ++j) {
                    bool isFastq, isMate;
                    isMate = optionPairFile[nowFileNumber+j].libraryType == "-op"
                           || optionPairFile[nowFileNumber+j].libraryType == "-OP" ? true : false;
                    // if read file is single file (forward and reverse seq contain same file)
                    if (optionPairFile[nowFileNumber+j].libraryType == "-ip" || optionPairFile[nowFileNumber+j].libraryType == "-op") {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        switch (fileFormat) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaSingleMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, numThread, isMate, isFastq);
                    // if read file is pair file (forward and reverse seq contain differenct file)
                    } else {
                        platanus::FILETYPE fileFormat1 = checkFileFormat(optionPairFile[nowFileNumber+j].fileFirst);
                        platanus::FILETYPE fileFormat2 = checkFileFormat(optionPairFile[nowFileNumber+j].fileSecond);
                        if (fileFormat1 != fileFormat2) {
                            throw platanus::FormatError("Different file type in paired-file.");
                        }
                        switch (fileFormat1) {
                        case platanus::UNKNOWN:
                            throw platanus::FormatError();
                            break;
                        case platanus::FASTA:
                            isFastq = false;
                            break;
                        default:
                            isFastq = true;
                            break;
                        }
                        ReadFastaPairMT(libraryMT[i], optionPairFile[nowFileNumber+j].fileFirst, optionPairFile[nowFileNumber+j].fileSecond, numThread, isMate, isFastq);
                    }
                }
            }
        } catch (platanus::ErrorBase &e) {
            e.showErrorMessage();
            exit(e.getID());
        }
    }
}

void Polish::pileupMappedReads(const vector<SeqLib> &library, const double minIdentity, const long numThread)
{
    std::cerr << "pileup mapped reads..." << std::endl;

	if (pileup.size() != contig.numSeq) {
		pileup.resize(contig.numSeq);
		for (unsigned long j = 0; j < pileup.size(); ++j)
			pileup[j].assign(contig.seq[j].length, pileupRecord());
	}

    omp_set_num_threads(numThread);
    #pragma omp parallel for schedule(static, 1) ordered
    for (long i = 0; i < numThread; ++i) {
        platanus::Position forwardPosition;
        platanus::Position reversePosition;
        platanus::SEQ forward;
        platanus::SEQ reverse;
		double forwardIdentity;
		double reverseIdentity;
		vector<vector<pileupRecord> > pileupInThread(contig.numSeq);
		for (unsigned long j = 0; j < pileupInThread.size(); ++j)
			pileupInThread[j].assign(contig.seq[j].length, pileupRecord());

        rewind (library[i].mappedFP);
        while (fread(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP)) {
            fread(&(forward.length), sizeof(long), 1, library[i].mappedFP);
            std::unique_ptr<char[]> tmp(new char[forward.length]());
            forward.base.resize(forward.length);
            fread(tmp.get(), sizeof(char), forward.length, library[i].mappedFP);
            std::copy(tmp.get(), tmp.get() + forward.length, forward.base.begin());
            fread(&forwardIdentity, sizeof(double), 1, library[i].mappedFP);

            fread(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            fread(&(reverse.length), sizeof(long), 1, library[i].mappedFP);
            tmp.reset(new char[reverse.length]());
            fread(tmp.get(), sizeof(char), reverse.length, library[i].mappedFP);
            reverse.base.resize(reverse.length);
            std::copy(tmp.get(), tmp.get() + reverse.length, reverse.base.begin());
            fread(&reverseIdentity, sizeof(double), 1, library[i].mappedFP);

			if (forwardPosition.id != 0) {
				if (forwardIdentity >= minIdentity)
					pileupSingleRead(forward, forwardPosition, pileupInThread);
				else
					pileupSingleOddRead(forward, forwardPosition, pileupInThread);
			}

			if (reversePosition.id != 0) {
				if (reverseIdentity >= minIdentity)
					pileupSingleRead(reverse, reversePosition, pileupInThread);
				else
					pileupSingleOddRead(reverse, reversePosition, pileupInThread);
			}
        }

		#pragma omp ordered
		for (unsigned long j = 0; j < pileupInThread.size(); ++j) {
			for (unsigned long k = 0; k < pileupInThread[j].size(); ++k) {
				pileup[j][k].numRead += pileupInThread[j][k].numRead;
				pileup[j][k].numOddRead += pileupInThread[j][k].numOddRead;
			}
		}
    }

/*
//for (unsigned long j = 0; j < pileup.size(); ++j) {
for (unsigned long j = 4; j < 5; ++j) {
	std::cout << '>' << j << std::endl;
	for (unsigned long k = 0; k < pileup[j].size(); ++k) {
		std::cout << pileup[j][k].numRead << '\t' << pileup[j][k].numOddRead << std::endl;
	}
}
*/
    return;
}


void Polish::pileupSingleRead(const platanus::SEQ &read, const platanus::Position &position, vector<vector<pileupRecord> > &result)
{
	long i, st, ed;

	if (position.id > 0) {
		i = position.id - 1;
		st = std::max(position.offset, 0);
		ed = std::min(position.offset + read.length, static_cast<long>(result[i].size()));
	}
	else {
		i = -position.id - 1;
		st = std::max(position.offset - read.length + 1, 0L);
		ed = std::min(position.offset + 1, static_cast<int>(result[i].size()));
	}

	for (long j = st; j < ed; ++j)
		++result[i][j].numRead;
}


void Polish::pileupSingleOddRead(const platanus::SEQ &read, const platanus::Position &position, vector<vector<pileupRecord> > &result)
{
	long i, st, ed;

	if (position.id > 0) {
		i = position.id - 1;
		st = std::max(position.offset, 0);
		ed = std::min(position.offset + read.length, static_cast<long>(result[i].size()));
	}
	else {
		i = -position.id - 1;
		st = std::max(position.offset - read.length + 1, 0L);
		ed = std::min(position.offset + 1, static_cast<int>(result[i].size()));
	}

	for (long j = st; j < ed; ++j)
		++result[i][j].numOddRead;
}


void Polish::maskErrorBases(const double minNumOddReadRatio)
{

	for (unsigned long i = 0; i < pileup.size(); ++i) {
		for (unsigned long j = 0; j < pileup[i].size(); ++j) {
			if (static_cast<double>(pileup[i][j].numOddRead) / pileup[i][j].numRead > minNumOddReadRatio)
				contig.seq[i].base[j] = 4;
		}
	}
}


void Polish::maskLowCoverageEdge(const unsigned long minCoverage)
{

	for (unsigned long i = 0; i < pileup.size(); ++i) {
/*
		unsigned long sum = 0;
		for (unsigned long j = 0; j < pileup[i].size(); ++j)
			sum += pileup[i][j].numRead;
		double meanCoverage = static_cast<double>(sum) / pileup[i].size();
*/
		for (unsigned long j = 0; j < pileup[i].size(); ++j) {
			if (pileup[i][j].numRead >= minCoverage)
				break;
			contig.seq[i].base[j] = 4;
		}
		for (unsigned long j = 0; j < pileup[i].size(); ++j) {
			if (pileup[i][pileup[i].size() - j - 1].numRead >= minCoverage)
				break;
			contig.seq[i].base[pileup[i].size() - j - 1] = 4;
		}
	}
}


void Polish::maskShortContig(platanus::SEQ &seq, const long minContigLength)
{
	long st = 0;
	long ed = 0;
	bool gapFlag = true;

	for (long i = 0; i < seq.length; ++i) {
		if (seq.base[i] != 4) {
			if (gapFlag) {
				st = i;
				ed = i + 1;
			}
			else {
				++ed;
			}
			gapFlag = false;
		}
		else {
			if (!gapFlag) {
				if (ed - st  < minContigLength) {
					for (long j = st; j < ed; ++j)
						seq.base[j] = 4;
				}
			}
			gapFlag = true;
		}
	}

	if (!gapFlag && ed - st  < minContigLength) {
		for (long j = st; j < ed; ++j)
			seq.base[j] = 4;
	}
}


void Polish::trimEdgeN(platanus::SEQ &seq)
{
	long i = 0;

	for (i = 0; i < seq.length; ++i) {
		if (seq.base[i] != 4)
			break;
	}
	long st = i;

	for (i = 0; i < seq.length; ++i) {
		if (seq.base[seq.length - i - 1] != 4)
			break;
	}
	long ed = seq.length - i;

	if (ed <= st) {
		seq.base.clear();
		seq.length = 0;
	}
	else {
		std::copy(seq.base.begin() + st, seq.base.begin() + ed, seq.base.begin());
		seq.length = ed - st;
		seq.base.resize(seq.length);
	}
}

void Polish::setFastaFileInfo(const vector<string> fastaFileName)
{
	platanus::setFastaFileNameAndNumber(fastaFileName, contigFileName, numContigInFile);
}

void Polish::printSeq(const string &outPrefix, const unsigned long long readLength, const unsigned long long contigMaxK)
{
	long seqID = 0;
	long numContigInOneFile = 0;
	long fileIndex = 0;

	unsigned long prefixPosition = contigFileName[fileIndex].find(outPrefix);
	if (prefixPosition != std::string::npos)
		prefixPosition += outPrefix.size();
	else
		prefixPosition = 0;
	if (!isalnum(contigFileName[fileIndex][prefixPosition]))
		++prefixPosition;
	string outName = outPrefix + "_polished_" + contigFileName[fileIndex].substr(prefixPosition);
    std::ofstream out(outName);


	for (long i = 0; i < contig.numSeq; ++i) {
		if (numContigInOneFile >= numContigInFile[fileIndex]) {
			++fileIndex;
			numContigInOneFile = 0;

			out.close();
			unsigned long prefixPosition = contigFileName[fileIndex].find(outPrefix);
			if (prefixPosition != std::string::npos)
				prefixPosition += outPrefix.size();
			else
				prefixPosition = 0;
			if (!isalnum(contigFileName[fileIndex][prefixPosition]))
				++prefixPosition;
			outName = outPrefix + "_polished_" + contigFileName[fileIndex].substr(prefixPosition);
			out.open(outName);
		}

		++seqID;
		++numContigInOneFile;

		if (contig.seq[i].length == 0)
			continue;

		unsigned long j = 0;
		for (; contig.name[i].size(); ++j) {
			if (isdigit(contig.name[i][j])) {
				for (;contig.name[i].size(); ++j) {
					if (!isdigit(contig.name[i][j]))
						break;
				}
				break;
			}
		}
		if (contig.name[i].size() > 0)
			out << '>' << contig.name[i].substr(0, j);
		else
			out << ">seq" << seqID;

		out << "_len" << contig.seq[i].length << "_cov" << contig.coverage[i] << "_read" << readLength << "_maxK" << contigMaxK << " pre_name:" << contig.name[i] << '\n';


		for (long j = 0; j < contig.seq[i].length; ++j) {
			out.put(platanus::Bin2Char(contig.seq[i].base[j]));
			if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
				out.put('\n');
		}
		if (contig.seq[i].length % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
			out.put('\n');
	}

    out.close();
}
