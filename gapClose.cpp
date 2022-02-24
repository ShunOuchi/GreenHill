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

#include "gapClose.h"
#include "gapCloseOLC.h"
#include "gapCloseDBG.h"
#include <deque>
#include <utility>
#include <algorithm>
#include <sstream>
#include <array>

using std::vector;
using std::string;


//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const unsigned short GapClose::HEAD_TAIL_SEQ_LEN = 150;
const unsigned short GapClose::SD_RATIO_MAPPED_GAP = 3;
const unsigned short GapClose::Closer::MIN_NUM_READS_COVERING_SMALL_GAP = 3;
const int GapClose::MIN_OVERLAP_FOR_CIRCLE = 100;
const unsigned GapClose::Closer::BRUIJN_MIN_KMER = 24;
const unsigned GapClose::Closer::BRUIJN_MAX_KMER = 72;
const unsigned GapClose::Closer::MIN_COVERAGE = 3;
const unsigned short GapClose::Closer::SD_RATIO_TOO_MANY_READS_ON_GAP = 3;

// not change value below parameter!!
const unsigned short GapClose::GAP_CANDIDATE_SEQ_1ST = 2;

//////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////
GapClose::GapClose()
{
    optionSingleArgs["-o"] = "out";
    optionSingleArgs["-k"] = "32";
    optionSingleArgs["-vo"] = "32";
    optionSingleArgs["-vd"] = "32";
    optionSingleArgs["-d"] = "5000";
    optionSingleArgs["-eo"] = "1";
    optionSingleArgs["-ed"] = "0.05";
    optionSingleArgs["-ro"] = "0.66";
    optionSingleArgs["-rs"] = "0.9";
    optionSingleArgs["-t"] = "1";
    optionBool["-a"] = false;
    optionBool["-n"] = false;
    optionBool["-extend"] = false;
    optionBool["-reduce_redundancy"] = false;
    optionBool["-no_partial"] = false;

    optionMultiArgs["-s"] = vector<string>(3);
    optionMultiArgs["-s"][0] = "32";
    optionMultiArgs["-s"][1] = "64";
    optionMultiArgs["-s"][2] = "96";

    optionMultiArgs["-c"] = vector<string>();
    optionMultiArgs["-f"] = vector<string>();
    pairedEndSingleFileType.emplace_back("-ip");
    pairedEndSingleFileType.emplace_back("-op");
    pairedEndPairFileType.emplace_back("-IP");
    pairedEndPairFileType.emplace_back("-OP");
    optionSingleArgs["-tmp"] = ".";
}




//////////////////////////////////////////////////////////////////////////////////////
// usage
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::usage(void) const
{

    std::cerr << "\nUsage: platanus_allee gap_close [Options]\n"
              << "Options:\n"
              << "    -o STR                             : prefix of output file (default " << optionSingleArgs.at("-o") << ", length <= " << platanus::ConstParam::MAX_FILE_LEN << ")\n"
              << "    -c FILE1 [FILE2 ...]               : scaffold_file (fasta format)\n"
              << "    -f FILE1 [FILE2 ...]               : single end files (fasta or fastq, number <= " << platanus::ConstParam::MAX_FILE_NUM << ")\n"
              << "    -ip{INT} PAIR1 [PAIR2 ...]         : lib_id inward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -op{INT} PAIR1 [PAIR2 ...]         : lib_id outward_pair_file (reads in 1 file, fasta or fastq)\n"
              << "    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)\n"
              << "    -s INT1 [INT2 ...]                 : mapping seed length for short reads (default " << optionMultiArgs.at("-s")[0] << " " << optionMultiArgs.at("-s")[1] << " " << optionMultiArgs.at("-s")[2]  << ")\n"
              << "    -k INT                             : mapping seed length in overlpap-layout-consensus (default " << optionSingleArgs.at("-k") << ")\n"
              << "    -d INT                             : maximum the number of read in gaps exec closing gap (default " << optionSingleArgs.at("-d") << ")\n"
              << "    -vo INT                            : minimum overlap length among each read in OLC gap closing (default " << optionSingleArgs.at("-vo") << ")\n"
              << "    -vd INT                            : minimum overlap length between contig and edge seq in De Bruijn gap closing (default " << optionSingleArgs.at("-vd") << ")\n"
              << "    -eo INT                            : maximum edit distance of overlap in OLC gap closing (identity, default " << optionSingleArgs.at("-eo") << ")\n"
              << "    -ed FLOAT                          : maximum error rate among gap edge seq in De Bruijn gap closing (identity, default " << optionSingleArgs.at("-ed") << ")\n"
              << "    -ro FLOAT                          : minimum consensus rate in OLC gap closing (identity, default " << optionSingleArgs.at("-ro") << ")\n"
              << "    -rs FLOAT                          : minimum consensus rate in Single Read gap closing (identity, default " << optionSingleArgs.at("-rs") << ")\n"
              << "    -a                                 : do gap close only one time using all libraries\n"
              << "    -n                                 : not extend the ends of scaffolds\n"
              << "    -t INT                             : number of threads (default " << optionSingleArgs.at("-t") << ")\n"
              << "    -extend                            : extend the ends of scaffolds (default false)\n"
              << "    -reduce_redundancy                 : reduce redundant sequences that exactly matche others (default, off)\n"
              << "    -no_partial                        : not close gaps partially, i.e. only close ones completely (default, off)\n"
              << "    -tmp DIR                           : directory for temporary files (default " << optionSingleArgs.at("-tmp") << ")\n\n\n"

              << "Output:\n"
//              << "    PREFIX_gapClosed_tmp[LIBRARY_ID].fa\n"
              << "    PREFIX_gapClosed.fa\n"
              << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// exec gap close
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::exec(void)
{
    const unsigned seedLength = atoi(optionMultiArgs["-s"][0].c_str());
    vector<int> multiSeedLength;
    multiSeedLength.clear();
	for (auto itr = optionMultiArgs["-s"].begin(); itr != optionMultiArgs["-s"].end(); ++itr) {
		unsigned length = std::stoi(*itr);
		multiSeedLength.push_back(length);
	}

    unsigned i;
    const int keyLength = std::min(seedLength, platanus::ConstParam::SCAFFOLD_HASH_OVERLAP);
    int numLibrary = 0;
    int numThread = atoi(optionSingleArgs["-t"].c_str());
	platanus::setGlobalTmpFileDir(optionSingleArgs["-tmp"].c_str());
    const long kmerLength = atoi(optionSingleArgs["-k"].c_str());
    const unsigned long long olcThreshold = atoi(optionSingleArgs["-d"].c_str());
    const unsigned minOverlapOLC = atoi(optionSingleArgs["-vo"].c_str());
    const unsigned minOverlapDe = atoi(optionSingleArgs["-vd"].c_str());
    const unsigned maxEditDistance = atoi(optionSingleArgs["-eo"].c_str());
    const double maxMissRate = atof(optionSingleArgs["-ed"].c_str());
    const double consensusThreshold = atof(optionSingleArgs["-ro"].c_str());
    const double consensusThresholdSingle = atof(optionSingleArgs["-rs"].c_str());
    vector<int> numFilePerLibraryID(optionPairFile.size());
    vector<int> libraryIDList(optionPairFile.size());
    vector<vector<SeqLib> > libraryMT;
    vector<SeqLib> singleLibrary;
    platanus::Contig contig;
    FILE *gapSeqFP;
    FILE *unusedFP;
    unsigned long long readLength = platanus::Contig::getReadLengthFromFastaHeader(optionMultiArgs["-c"][0]);
	unsigned long long contigMaxK = platanus::Contig::getMaxKFromFastaHeader(optionMultiArgs["-c"][0]);
    unsigned long long totalReadLength = 0;
    unsigned long long totalNumRead = 0;
    if (readLength == 0)
        readLength = platanus::ConstParam::DEFAULT_CONTIG_READ_LEN;


    // nest {} for calling mapper destructor
    {
        Mapper mapper(seedLength, keyLength);
		mapper.setMultiSeedLength(multiSeedLength);

        sort(optionPairFile.begin(), optionPairFile.end());

        // count up the number of file in ecah library
        for (i = 0; i < optionPairFile.size(); ++i) {
            ++(numFilePerLibraryID[numLibrary]);
            libraryIDList[numLibrary] = optionPairFile[i].libraryID;
            if (i + 1 < optionPairFile.size() && optionPairFile[i].libraryID != optionPairFile[i + 1].libraryID) {
                ++numLibrary;
            }
        }
		if (i > 0)
			++numLibrary;

        libraryMT.resize(numLibrary);
        omp_set_num_threads(numThread);

        readLibrary(mapper, contig, libraryMT, singleLibrary, numFilePerLibraryID, numLibrary, numThread);

        // estimate library insert size
        gapSeqFP = platanus::makeTemporaryFile();


        if (optionMultiArgs["-f"].size() > 0) {
            std::cerr << "[SINGLE_LIBRARY]" << std::endl;
            mapper.mapSmallGap(singleLibrary, gapSeqFP, numThread);
            vector<SeqLib>().swap(singleLibrary);
        }


        for (int i = 0; i < numLibrary; ++i) {
            std::cerr << "[PAIR_LIBRARY " << libraryIDList[i] << "]" << std::endl;
            libraryMT[i][0].setAverageLength((long)((double)(libraryMT[i][0].getTotalLength()) / (2 * libraryMT[i][0].getNumPair()) + 0.5));

            mapper.gatherPairReadMappedSameContig(libraryMT[i], numThread);

            vector<long> insSizeDistribution;
            libraryMT[i][0].readInsertSizeFile(insSizeDistribution);
            libraryMT[i][0].estimateInsSize(insSizeDistribution, i > 0 ? libraryMT[i - 1][0].getAverageInsSize() : 0, platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR, platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR);

            mapper.mapSmallGap(libraryMT[i], gapSeqFP, numThread);
        }

//        sort(libraryMT.begin(), libraryMT.end());

        mapper.releaseContig(contig);
    }


    Closer closer(contig, olcThreshold, minOverlapOLC, minOverlapDe, maxEditDistance, maxMissRate, consensusThreshold);
    contig.clear();
    closer.makeGapTable();

    // gap close using single read
    closer.closeSmallGaps(consensusThresholdSingle, gapSeqFP);
    fclose(gapSeqFP);

    unusedFP = platanus::makeTemporaryFile();
    // gap close using pair read library
    for (int i = 0; i < numLibrary; ++i) {
        std::cerr << "[PAIR_LIBRARY " << i + 1 << "]" << std::endl;
        gapSeqFP = closer.saveGapCoveringReads(libraryMT[i], numThread);

        for (long j = 0; j < numThread; ++j) {
            fclose(libraryMT[i][j].mappedFP);
            libraryMT[i][j].mappedFP = NULL;
        }
        const long librarySD = libraryMT[i][0].getSDInsSize();
        libraryMT[i].clear();

        totalReadLength += libraryMT[i][0].getTotalLength();
        totalNumRead += 2 * libraryMT[i][0].getNumPair();

        closer.loadLocalReads(gapSeqFP);
        fclose(gapSeqFP);

        unsigned long long numBase = closer.calcMaxGapSeqBase(numThread);

        if (!optionBool["-a"]) {
            closer.gapCloseUsingPairReads(kmerLength, librarySD, numThread, optionBool["-no_partial"]);
//            std::ostringstream outStream;
//            outStream << optionSingleArgs["-o"] << "_gapClosed_tmp" << i + 1 << ".fa";
//            string outName = outStream.str();
//            closer.printSeq(outName);
        }

        //if (numLibrary > 1)
        closer.saveUnusedReads(numBase, unusedFP);
    }

//    double averageReadLength = static_cast<double>(totalReadLength) / totalNumRead;

    // gap close using all library
    if (numLibrary > 1) {
        fprintf(stderr, "[ALL LIBRARY]\n");
        closer.loadUnusedReads(unusedFP);
        closer.gapCloseUsingPairReads(kmerLength, 0, numThread, optionBool["-no_partial"], true);
    }

	closer.generateGapClosedSeq(!optionBool["-extend"], contigMaxK);
//	closer.findCircularGapClosedSeq(numThread);

	closer.setFastaFileInfo(optionMultiArgs["-c"]);

	if (optionBool["-reduce_redundancy"]) {
		closer.setHeteroInfoOfGapClosedSeq();
		closer.markRedundantSeq(contigMaxK, numThread);
	}

	closer.printGapClosedSeq(optionSingleArgs["-o"], readLength, contigMaxK, optionBool["-extend"]);

    if (numLibrary == 1)
		closer.loadUnusedReads(unusedFP);

    fclose(unusedFP);

//    FILE *contigFP = closer.localAssemble(numThread);
//    string outName = optionSingleArgs["-o"] + "_extraContig.fa";
//    closer.printContig(outName, contigFP, averageReadLength, readLength, contigMaxK);

    std::cerr << "gap_close completed!!" << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// Read library files
// Change read function which fasta or fastq, single file or pair one.
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::readLibrary(Mapper &mapper, platanus::Contig &contig, vector<vector<SeqLib> > &libraryMT, vector<SeqLib> &singleLibrary, vector<int> &numFilePerLibraryID, const int numLibrary, const int numThread)
{
    unsigned nowFileNumber = 0;
    int j;
    omp_set_num_threads(numThread);
    #    pragma omp parallel for private(j, nowFileNumber) schedule(static, 1)
    for (int i = -2; i < numLibrary; ++i) {
        try {
            if (i == -2) {
                // make contig(in fact, scaffold) kmer mapper table
                for (unsigned i = 0; i < optionMultiArgs["-c"].size(); ++i)
                    contig.readFastaCoverage(optionMultiArgs["-c"][i]);
//					contig.readFastaCoverageAddNEdge(optionMultiArgs["-c"][i]);
                mapper.setContig(contig);
                mapper.makeKmerTable();
            // read single read library
            } else if (i == -1) {
                if (optionMultiArgs["-f"].size() > 0) {
                    singleLibrary.resize(numThread);
                    for (j = 0; j < numThread; ++j) {
                        singleLibrary[j].pairFP = platanus::makeTemporaryFile();
                    }
                    for (j = 0; j < static_cast<long>(optionMultiArgs["-f"].size()); ++j) {
                        platanus::FILETYPE fileFormat = checkFileFormat(optionMultiArgs["-f"][j]);
                        if (fileFormat == platanus::FILETYPE::FASTA)
                            ReadFastaSingleMT(singleLibrary, optionMultiArgs["-f"][j], numThread, false, false, true);
                        else if (fileFormat == platanus::FILETYPE::FASTQ)
                            ReadFastaSingleMT(singleLibrary, optionMultiArgs["-f"][j], numThread, false, true, true);
                        else {
                            throw platanus::FormatError();
                        }
                    }
                }

            // read pair read library
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


//////////////////////////////////////////////////////////////////////////////////////
// make gap list
// gap list contains each gaps below information
// - scaffoldID, and position in scaffold
// - candidate seqs which fill gap
// - gap's terminal seq
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::makeGapTable(void)
{
    long totalGapLength = 0;
    long numScaffold = this->scaffold.size();
    long numGap = 0;
    long i, j;


    std::cerr << "making hash table of gaps..." << std::endl;

    // calculate total number of gaps
    for (i = 0; i < numScaffold; ++i) {
        j = 0;
        while (j < scaffold[i].length) {
            if (scaffold[i].base[j] == 4) {
                while (scaffold[i].base[j] == 4 && j < scaffold[i].length) {
                    ++totalGapLength;
                    ++j;
                }
                ++numGap;
            }
            ++j;
        }
    }
    gap.resize(numGap);
    numGap = 0;
    for (i = 0; i < numScaffold; ++i) {
        j = 0;
        while (j < scaffold[i].length) {
            if (scaffold[i].base[j] == 4) {
                gap[numGap].scaffoldID = i + 1;
                gap[numGap].start = j;
                gap[numGap].state = GAP::UNCLOSED;

                // head seq
                int length = std::min(j, static_cast<long>(HEAD_TAIL_SEQ_LEN));
                gap[numGap].headSeq.resize(length);
                int cutLength = 0;
                for (int k = 0; k < length; ++k) {
                    if (scaffold[i].base[j-length+k] == 4) {
                        cutLength = k + 1;
                    } else {
                        gap[numGap].headSeq[k-cutLength] = scaffold[i].base[j-length+k];
                    }
                }

                gap[numGap].headSeq.resize(length - cutLength);

                while (scaffold[i].base[j] == 4 && j < scaffold[i].length) {
                    this->insertGap(gap[numGap].scaffoldID, j, numGap);
                    ++j;
                }

                gap[numGap].end = j;
                // tail seq
                length = std::min(scaffold[i].length - j, static_cast<long>(HEAD_TAIL_SEQ_LEN));
                gap[numGap].tailSeq.resize(length);
                for (int k = 0; k < length; ++k) {
                    if (scaffold[i].base[j+k] == 4) {
                        gap[numGap].tailSeq.resize(k);
                        break;
                    } else {
                        gap[numGap].tailSeq[k] = scaffold[i].base[j+k];
                    }
                }

                gap[numGap].remainedGap = gap[numGap].end - gap[numGap].start;

                ++numGap;
            }
            ++j;
        }
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// insert gaptable
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::insertGap(const int id, const int offset, const long value)
{
    unsigned long long key = (static_cast<unsigned long long>(id) << 32) | static_cast<unsigned long long>(offset);
    gapTable[key] = value;
}


//////////////////////////////////////////////////////////////////////////////////////
// find gaptable
//////////////////////////////////////////////////////////////////////////////////////
long GapClose::Closer::findGapID(const int id, const int offset)
{
    unsigned long long key = (static_cast<unsigned long long>(id) << 32) | static_cast<unsigned long long>(offset);
    auto tableIterator = gapTable.find(key);
    if (tableIterator != gapTable.end()) {
        return tableIterator->second;
    } else {
        return -1;
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// gather candidate seq which mapped near gap on library
//////////////////////////////////////////////////////////////////////////////////////
FILE *GapClose::Closer::saveGapCoveringReads(const vector<SeqLib> &library, const long numThread)
{
    long tolerance = library[0].getSDInsSize() * SD_RATIO_MAPPED_GAP;
    FILE *tmpFP[numThread];

    std::cerr << "saving reads covering gaps..." << std::endl;
    omp_set_num_threads(numThread);
    #    pragma omp parallel for  schedule(static, 1)
    for (long i = 0; i < numThread; ++i) {
        platanus::Position forwardPosition;
        platanus::Position reversePosition;
        platanus::SEQ forward;
        platanus::SEQ reverse;
        const long averageInsSize = library[0].getAverageInsSize();

        tmpFP[i] = platanus::makeTemporaryFile();

        rewind (library[i].mappedFP);
        while (fread(&forwardPosition, sizeof(platanus::Position), 1, library[i].mappedFP)) {
            // load forward and reverse seq
            fread(&(forward.length), sizeof(long), 1, library[i].mappedFP);
            std::unique_ptr<char[]> tmp(new char[forward.length]());
            fread(tmp.get(), sizeof(char), forward.length, library[i].mappedFP);
            forward.base.resize(forward.length);
            std::copy(tmp.get(), tmp.get() + forward.length, forward.base.begin());
            fread(&reversePosition, sizeof(platanus::Position), 1, library[i].mappedFP);
            fread(&(reverse.length), sizeof(long), 1, library[i].mappedFP);
            tmp.reset(new char[reverse.length]());
            fread(tmp.get(), sizeof(char), reverse.length, library[i].mappedFP);
            reverse.base.resize(reverse.length);
            std::copy(tmp.get(), tmp.get() + reverse.length, reverse.base.begin());

            judgePairReadMappedNearGap(forwardPosition, reverse, averageInsSize, tolerance, tmpFP[i]);
            judgePairReadMappedNearGap(reversePosition, forward, averageInsSize, tolerance, tmpFP[i]);
        }
    }
    for (long i = 1; i < numThread; ++i) {
        char c;
        rewind(tmpFP[i]);
        while (fread(&c, sizeof(char), 1, tmpFP[i]))
            fwrite(&c, sizeof(char), 1, tmpFP[0]);
        fclose(tmpFP[i]);
    }
    return tmpFP[0];
}




//////////////////////////////////////////////////////////////////////////////////////
// judge whether seq mapped near gap
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::judgePairReadMappedNearGap(platanus::Position &mappedPosition, platanus::SEQ &pairSeq, const long averageInsSize, const long tolerance, FILE *tmpFP)
{
    long start, end;
    if (mappedPosition.id != 0) {
        if (mappedPosition.id > 0) {
            start = std::max(static_cast<long>(mappedPosition.offset), mappedPosition.offset + averageInsSize - tolerance - pairSeq.length);
            start = std::min(start, scaffold[mappedPosition.id - 1].length - 1);
            start = std::max(start, static_cast<long>(0));
            end = std::min(mappedPosition.offset + averageInsSize + tolerance, scaffold[mappedPosition.id - 1].length);
            pairSeq.reverse();
        }
        else {
            mappedPosition.id *= -1;
            start = std::min(static_cast<long>(mappedPosition.offset), mappedPosition.offset - averageInsSize - tolerance);
            start = std::max(start, static_cast<long>(0));
            end = std::min(mappedPosition.offset - averageInsSize + tolerance + pairSeq.length, scaffold[mappedPosition.id - 1].length);
        }
        while (start < end) {
            if (scaffold[mappedPosition.id - 1].base[start] == 4) {
                long gapID = this->findGapID(mappedPosition.id, start);
                if (gap.at(gapID).length == 0) {
                    fwrite(&gapID, sizeof(long), 1, tmpFP);
                    fwrite(&(pairSeq.length), sizeof(long), 1, tmpFP);
                    fwrite(pairSeq.base.c_str(), sizeof(char), pairSeq.length, tmpFP);
                }
                while (scaffold[mappedPosition.id - 1].base[start] == 4 && start < end)
                ++start;
            }
            ++start;
        }
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// judge whether seq mapped near gap
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::loadLocalReads(FILE *gapSeqFP)
{
    long gapID;
    long length;
    unsigned numGap = gap.size();

    std::cerr << "loading reads covering gaps..." << std::endl;
    rewind(gapSeqFP);
    // calculate gap length
    while (fread(&gapID, sizeof(long), 1, gapSeqFP)) {
        fread(&length, sizeof(long), 1, gapSeqFP);
        gap[gapID].length += length + 1;
        fseek(gapSeqFP, length * sizeof(char), SEEK_CUR);
    }

    // push head and tail seqs
    for (unsigned i = 0; i < numGap; ++i) {
        if (gap[i].state == GAP::CLOSED) continue;
        gap[i].seq.emplace_back(gap[i].headSeq);
        gap[i].seq.emplace_back(gap[i].tailSeq);
    }

    // push candidate seqs
    rewind(gapSeqFP);
    while (fread(&gapID, sizeof(long), 1, gapSeqFP)) {
        fread(&length, sizeof(long), 1, gapSeqFP);
        gap[gapID].totalSeqBase += length;
        std::unique_ptr<char[]> tmp(new char[length]());
        fread(tmp.get(), sizeof(char), length, gapSeqFP);
        vector<char> str(length);
        std::copy(tmp.get(), tmp.get() + length, str.begin());
        gap[gapID].seq.emplace_back(std::move(str));
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// exec gap close actually using overlap layout consensus
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::gapCloseUsingPairReads(const long keyLength, const long sd, const long numThread, const bool noPartialFlag, const bool final)
{

    typedef unsigned long long u64_t;
    u64_t numGap = 0;
    u64_t numGapDe = 0;
    u64_t numGapOLC = 0;
    u64_t numClosedDe = 0;
    u64_t numClosedDePartial = 0;
    u64_t numClosedOLC = 0;
    u64_t too = 0;
    u64_t size = this->gap.size();

    vector<GAP *> buffer(size);
    std::deque<bool> tooManyReadsThisLibrary(size, false);

    std::cerr << "assembling localized reads..." << std::endl;

    omp_set_num_threads(numThread);

    for (unsigned long long i = 0; i < size; ++i)
        buffer[i] = &(this->gap[i]);

    sort(buffer.begin(), buffer.end(), GAP());

    // ignore too many reads mappend on gap
    if (!final) {
        vector<double> coverageRatio(size);
        double average = 0;
        double standarddev = 0;
        unsigned long long num = 0;
        // # pragma omp parallel for schedule(static, 1) reduction(+: average)
        for (unsigned long long i = 0; i < size; ++i) {
            if (buffer[i]->state == GAP::UNCLOSED) {
                coverageRatio[i] = this->tooManyReadsGathered(buffer[i], sd);
                average += coverageRatio[i];
                ++num;
            }
        }
        average /= num;
        # pragma omp parallel for schedule(static, 1) reduction(+: standarddev)
        for (unsigned long long i = 0; i < size; ++i) {
            if (buffer[i]->state == GAP::UNCLOSED) {
                standarddev += (coverageRatio[i] - average) * (coverageRatio[i] - average);
            }
        }
        standarddev = sqrt(standarddev / (num - 1));
        average += standarddev * SD_RATIO_TOO_MANY_READS_ON_GAP;
        // # pragma omp parallel for schedule(static, 1)
        /*
        for (unsigned long long i = 0; i < size; ++i) {
            if (coverageRatio[i] > average) {
                buffer[i]->tooManyReads = true;
                tooManyReadsThisLibrary[i] = true;
            }
        }
        */
    }


    #    pragma omp parallel for schedule(dynamic) reduction(+: numGap, numGapDe, numGapOLC, numClosedDe, numClosedDePartial, numClosedOLC, too)
    for (unsigned long long i = 0; i < size; ++i) {
        bool closed = false;
        if (buffer[i]->state == GAP::CLOSED) continue;
        ++numGap;
        /*
        if ((final && buffer[i]->tooManyReads) || (!final && tooManyReadsThisLibrary[i])) {
            ++too;
            continue;
        }
        */

        // if the num of read exceed olcThreshold, use De Bruijn Graph

        ++numGapDe;
        buffer[i]->sortSeq();
        GapCloseDBG<BRUIJN_MIN_KMER, BRUIJN_MAX_KMER> ele(BRUIJN_MIN_KMER, BRUIJN_MAX_KMER, MIN_COVERAGE, minOverlapDe, maxMissRate, buffer[i]);
        if (buffer[i]->start != 0 && buffer[i]->end != this->scaffold[buffer[i]->scaffoldID - 1].length) {
			if (ele.enableGapClose()) {
				ele.gapAssemble();
				closed = ele.isClosedGap();
				if (closed) {
					++numClosedDe;
					buffer[i]->closingSeq = ele.getSeq();
					buffer[i]->closingLength = ele.getSeqLength();
					buffer[i]->state = GAP::CLOSED;
				}
				else if (buffer[i]->seq.size() <= olcThreshold) {
					++numGapOLC;
					OverlapLayoutConsensus olc(keyLength, buffer[i]);
					olc.makeKmerTable();
					olc.makeOverlapGraph(minOverlapOLC, maxMissRate);
					if (olc.hasHeadSeqLink()) {
						closed = olc.greedyExtension(minConsensusRate);
						if (closed) {
							++numClosedOLC;
							buffer[i]->closingSeq = olc.getLayoutSeq();
							buffer[i]->closingLength = olc.getLayoutLength();
							buffer[i]->state = GAP::CLOSED;
						}
					}
				}

				if (!noPartialFlag && buffer[i]->state != GAP::CLOSED) {
					closed = ele.isClosedGapPartial();
					if (closed && (buffer[i]->state == GAP::UNCLOSED || (buffer[i]->state == GAP::PART_CLOSED && ele.getRemainedGap() < buffer[i]->remainedGap))) {
						++numClosedDePartial;
						buffer[i]->closingSeq = ele.getSeq();
						buffer[i]->closingLength = ele.getSeqLength();
						buffer[i]->remainedGap = ele.getRemainedGap();
						buffer[i]->state = GAP::PART_CLOSED;
					}
				}

			}
        } else {
            ele.gapAssemble();
            if (ele.extendEdge()) {
                buffer[i]->closingSeq = ele.getSeq();
                buffer[i]->closingLength = ele.getSeqLength();
                buffer[i]->state = GAP::CLOSED;
            }
        }
    }

    std::cerr << "NUM_GAPS = " << numGap << std::endl;
//    std::cerr << "NUM_GAPS_IN_DE_BRUIJN_GAP_CLOSE = " << numGapDe << std::endl;
//    std::cerr << "NUM_GAPS_IN_OVERLAP_LAYOUT_CONSENSUS = " << numGapOLC << std::endl;
    std::cerr << "NUM_NOT_CLOSED_GAPS_BECAUSE_OF_TOO_MANY_READS = " << too << std::endl;
    std::cerr << "NUM_CLOSED_GAPS_IN_DE_BRUIJN = " << numClosedDe << std::endl;
    std::cerr << "NUM_CLOSED_GAPS_PARTIAL = " << numClosedDePartial << std::endl;
    std::cerr << "NUM_CLOSED_GAPS_IN_OVERLAP_LAYOUT_CONSENSUS = " << numClosedOLC << std::endl;
}


//////////////////////////////////////////////////////////////////////////////////////
// calc rate between average coverage in gap and scaffold which has this gap coverage
//////////////////////////////////////////////////////////////////////////////////////
double GapClose::Closer::tooManyReadsGathered(const GAP * const gap, const long sd) const
{
    return static_cast<double>(gap->totalSeqBase) / ((std::abs(gap->end - gap->start) + 1) + (2 * SD_RATIO_MAPPED_GAP * sd)) / this->scaffoldCoverage[gap->scaffoldID - 1];
}




//////////////////////////////////////////////////////////////////////////////////////
// save candidate seq
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::saveUnusedReads(const unsigned long long maxBase, FILE *unusedFP)
{
    fseek(unusedFP, 0, SEEK_END);
    std::unique_ptr<char[]> tmp(new char[maxBase]());
    for (long gapID = 0; gapID < static_cast<long>(gap.size()); ++gapID) {
        if (gap[gapID].state == GAP::CLOSED) continue;
        fwrite(&gapID, sizeof(long), 1, unusedFP);
        unsigned total = 0;
        for (unsigned seqID = GapClose::GAP_CANDIDATE_SEQ_1ST; seqID < gap[gapID].seq.size(); ++seqID) {
            std::move(gap[gapID].seq[seqID].begin(), gap[gapID].seq[seqID].end(), &tmp[total]);
            total += gap[gapID].seq[seqID].size() + 1;
            tmp[total - 1] = 'N';
        }
        fwrite(&total, sizeof(int), 1, unusedFP);
        fwrite(tmp.get(), sizeof(char), total, unusedFP);
    }

    # pragma omp parallel for schedule(dynamic)
    for (long gapID = 0; gapID < static_cast<long>(gap.size()); ++gapID) {
        if (gap[gapID].state != GAP::CLOSED) {
			gap[gapID].length = 0;
			gap[gapID].totalSeqBase = 0;
			gap[gapID].seq.clear();
		}
    }
}

//////////////////////////////////////////////////////////////////////////////////////
// clear gap
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::clearGap(unsigned gapID)
{
    gap[gapID].length = 0;
    gap[gapID].totalSeqBase = 0;
    gap[gapID].seq.clear();
    gap[gapID].state = GAP::UNCLOSED;
}



//////////////////////////////////////////////////////////////////////////////////////
// load candidate seq
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::loadUnusedReads(FILE *unusedFP)
{

    for (unsigned gapID = 0; gapID < gap.size(); ++gapID) {
        if (gap[gapID].state == GAP::CLOSED) continue;
        gap[gapID].seq.emplace_back(gap[gapID].headSeq);
        gap[gapID].seq.emplace_back(gap[gapID].tailSeq);
    }

    rewind(unusedFP);
    long id = 0;
    int length = 0;
    int it = 0;
    while (fread(&id, sizeof(long), 1, unusedFP)) {
        fread(&length, sizeof(int), 1, unusedFP);
        if (gap[id].state == GAP::CLOSED) {
            fseek(unusedFP, length * sizeof(char), SEEK_CUR);
            continue;
        }
        gap[id].length += length + 1;
        std::unique_ptr<char[]> tmpSeq(new char[length]());
        fread(tmpSeq.get(), sizeof(char), length, unusedFP);
        vector<char> str(length);
        for (int i = 0; i < length; ++i) {
            if (tmpSeq[i] != 'N') {
                str[it] = tmpSeq[i];
                ++it;
            } else {
                str.resize(it);
                gap[id].seq.emplace_back(str);
                str.clear();
                str.resize(length);
                it = 0;
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// single read gap closer
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::closeSmallGaps(const double consensusRate, FILE *gapSeqFP)
{
    platanus::Position position;
    long length;
    long numClosed = 0;
    const unsigned long long size = this->gap.size();
    vector<std::multiset<std::pair<vector<char>, int> > > tmpGap(size);
    typedef std::pair<vector<char>, int>  READ;
    vector<char> READ::* readSeq = &READ::first;
    int READ::* readLength = &READ::second;
    READ consensus;

    std::cerr << "making consensus sequences to close small gaps..." << std::endl;


    rewind(gapSeqFP);
    while (fread(&position, sizeof(platanus::Position), 1, gapSeqFP)) {
        long gapID = this->findGapID(position.id, position.offset);
        fread(&length, sizeof(long), 1, gapSeqFP);
        if (length > 0) {
            vector<char> tmp(length);
            fread(&tmp[0], sizeof(char), length, gapSeqFP);
            tmpGap[gapID].insert(std::make_pair(tmp, static_cast<int>(length)));
        } else {
            vector<char> tmp;
            tmpGap[gapID].insert(std::make_pair(tmp, static_cast<int>(length)));
        }
    }

    for (unsigned long long gapID = 0; gapID < size; ++gapID) {
        if (tmpGap[gapID].size() < MIN_NUM_READS_COVERING_SMALL_GAP) continue;
        if (decideConsensusFromReads(tmpGap[gapID], consensusRate, consensus, gapID)) {
            if (consensus.*readLength > 0) {
                this->gap[gapID].closingLength = consensus.*readLength;
                this->gap[gapID].closingSeq = consensus.*readSeq;
            } else {
                int i = 0;
                for (; i < -(consensus.*readLength); ++i) {
                    if (this->scaffold[this->gap[gapID].scaffoldID - 1].base[this->gap[gapID].start - i - 1] != this->scaffold[this->gap[gapID].scaffoldID - 1].base[this->gap[gapID].end - consensus.*readLength - i - 1]
                        || this->scaffold[this->gap[gapID].scaffoldID - 1].base[this->gap[gapID].start - i - 1] == 4
                        || this->scaffold[this->gap[gapID].scaffoldID - 1].base[this->gap[gapID].end - consensus.*readLength - i - 1] == 4)
                        break;
                }
                if (i != -(consensus.*readLength)) continue;

                this->gap[gapID].closingLength = consensus.*readLength - 1;
                this->gap[gapID].closingSeq.clear();
            }
            ++numClosed;
            this->gap[gapID].state = GAP::CLOSED;
        }

    }

    std::cerr << "NUM_GAP=" << this->gap.size() << std::endl;
    std::cerr << "NUM_CLOSED_GAP=" << numClosed << std::endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// exec overlap layout consensus like gap closing actually in single read mode
//////////////////////////////////////////////////////////////////////////////////////
bool GapClose::Closer::decideConsensusFromReads(const std::multiset<std::pair<vector<char>, int> > &seq, const double threshold, std::pair<vector<char>, int> &consensus, const unsigned long long id)
{
    vector<int> positiveLengthCount(platanus::ConstParam::MAX_READ_LEN + 1, 0);
    vector<int> negativeLengthCount(platanus::ConstParam::MAX_READ_LEN + 1, 0);
    int baseCount[5];
    int most = 0;
    int mostLength = 0;

    consensus.first.clear();
    for (auto it = seq.begin(), end = seq.end(); it != end; ++it) {
        if (it->second >= 0) {
            ++positiveLengthCount[it->second];
        } else {
            ++negativeLengthCount[-(it->second)];
        }
    }
    /*
    for (long seqID = 0, n = seq.size(); seqID < n; ++seqID) {
        if (seq[seqID].second >= 0) {
            ++positiveLengthCount[seq[seqID].second];
        } else {
            ++negativeLengthCount[-seq[seqID].second];
        }
    }
    */
    for (long i = 0, n = positiveLengthCount.size(); i < n; ++i) {
        if (positiveLengthCount[i] > most) {
            most = positiveLengthCount[i];
            mostLength = i;
        }
        if (negativeLengthCount[i] > most) {
            most = negativeLengthCount[i];
            mostLength = -i;
        }
    }

    if (static_cast<double>(most) / seq.size() < threshold) return false;

    if (mostLength <= 0) {
        consensus.second = mostLength;
        return true;
    }

    most = 0;
    for (long i = 0; i < mostLength; ++i) {
        for (unsigned base = 0; base < 4; ++base) {
            baseCount[base] = 0;
        }
        /*
        for (unsigned seqID = 0, n = seq.size(); seqID < n; ++seqID) {
            if (seq[seqID].second == mostLength) {
                ++baseCount[static_cast<int>(seq[seqID].first[i])];
            }
        }
        */
        for (auto it = seq.begin(), end = seq.end(); it != end; ++it) {
            if (it->second == mostLength) {
                ++baseCount[static_cast<int>(it->first[i])];
            }
        }
        int tmp = 0;
        char tmpBase = 0;
        for (unsigned base = 0; base < 4; ++base) {
            if (baseCount[base] > tmp) {
                tmp = baseCount[base];
                tmpBase = base;
            }
        }
        consensus.first.emplace_back(tmpBase);
        most += tmp;
    }

    if (static_cast<double>(most) / (mostLength * positiveLengthCount[mostLength]) < threshold)
        return false;
    consensus.second = consensus.first.size();
    return true;
}







//////////////////////////////////////////////////////////////////////////////////////
// output
//////////////////////////////////////////////////////////////////////////////////////
void GapClose::Closer::printSeq(const string &outName, const unsigned long long readLength, const bool addN)
{
    long gapID = 0;
    long numClosed = 0;
    std::ofstream out(outName.c_str());


    for (unsigned long long scafID = 0 ; static_cast<long>(scafID) < numScaffold; ++scafID) {
        out << ">scaffold" << scafID + 1 << "_cov" << scaffoldCoverage[scafID] << "_read" << readLength << "\n";
        long position = addN ? 1 : 0;
        long end = addN ? scaffold[scafID].length - 1 : scaffold[scafID].length;
        position = 0;
        end = scaffold[scafID].length;
        unsigned mod = 0;
        while (position < end) {
            if (scaffold[scafID].base[position] != 4) {
                out.put(platanus::Bin2Char(scaffold[scafID].base[position]));
                ++mod;
                if (mod % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                    out.put('\n');
                ++position;
                continue;
            } else if ((gap[gapID].start == 0 || gap[gapID].end == scaffold[gap[gapID].scaffoldID - 1].length) && gap[gapID].state == GAP::UNCLOSED) {
                ++position;
                for (; position < gap[gapID].end; ++position) {
                    out.put('N');
                    ++mod;
                    if (mod % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                        out.put('\n');
                }
                gapID = std::min(gapID + 1, static_cast<long>(gap.size() - 1));
                continue;
            }
            if (gap[gapID].length > 0) {
                if (gap[gapID].state != GAP::UNCLOSED) {
                    for (long i = 0; i < gap[gapID].length; ++i) {
                        out.put(platanus::Bin2Char(gap[gapID].seq[0][i]));
                        ++mod;
                        if (mod % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                            out.put('\n');
                    }
                    position += gap[gapID].end - gap[gapID].start;
                    ++numClosed;
                } else {
                    for (; position < gap[gapID].end; ++position) {
                        out.put('N');
                        ++mod;
                        if (mod % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                            out.put('\n');
                    }
                }
            } else if (gap[gapID].length < 0) {
                position += gap[gapID].end - gap[gapID].start - (gap[gapID].length + 1);
                ++numClosed;
            } else {
                for (; position < gap[gapID].end; ++position) {
                    out.put('N');
                    ++mod;
                    if (mod % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                        out.put('\n');
                }
            }
            gapID = std::min(gapID + 1, static_cast<long>(gap.size() - 1));
        }
        if (mod % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
            out.put('\n');
    }
    out.close();

    std::cerr << "TOTAL_NUM_CLOSED_GAPS = " << numClosed << std::endl;
}


void GapClose::Closer::generateGapClosedSeq(const bool addN, const long contigMaxK)
{
    long gapID = 0;
    long numClosed = 0;

	gapClosedSeq.resize(numScaffold);

    for (unsigned long long scafID = 0 ; static_cast<long>(scafID) < numScaffold; ++scafID) {
		gapClosedSeq[scafID].coverage = scaffoldCoverage[scafID];
		gapClosedSeq[scafID].baseName = scaffoldName[scafID];

        long position = addN ? 1 : 0;
        long end = addN ? scaffold[scafID].length - 1 : scaffold[scafID].length;
        position = 0;
        end = scaffold[scafID].length;
        while (position < end) {
            if (scaffold[scafID].base[position] != 4) {
				gapClosedSeq[scafID].seq.push_back(scaffold[scafID].base[position]);
                ++position;
                continue;
            } else if ((gap[gapID].start == 0 || gap[gapID].end == scaffold[gap[gapID].scaffoldID - 1].length) && gap[gapID].state == GAP::UNCLOSED) {
                ++position;
                for (; position < gap[gapID].end; ++position) {
					gapClosedSeq[scafID].seq.push_back(4);
                }
                gapID = std::min(gapID + 1, static_cast<long>(gap.size() - 1));
                continue;
            }
            if (gap[gapID].closingLength > 0) {
                if (gap[gapID].state != GAP::UNCLOSED) {
					if (addN) {
						if (gap[gapID].start == 0) {
							gapClosedSeq[scafID].leftExtendedLen = gap[gapID].length;
						}
						if (gap[gapID].end == static_cast<int>(scaffold[scafID].base.size())) {
							gapClosedSeq[scafID].rightExtendedLen = gap[gapID].closingLength;
						}
					}

                    for (long i = 0; i < gap[gapID].closingLength; ++i) {
						gapClosedSeq[scafID].seq.push_back(gap[gapID].closingSeq[i]);
                    }
                    position += gap[gapID].end - gap[gapID].start;
                    ++numClosed;
                } else {
                    for (; position < gap[gapID].end; ++position) {
						gapClosedSeq[scafID].seq.push_back(4);
                    }
                }
            } else if (gap[gapID].closingLength < 0) {
				if ((gap[gapID].start + gap[gapID].closingLength <= contigMaxK) || (scaffold[scafID].length - gap[gapID].end + gap[gapID].closingLength) <= contigMaxK) {
					++position;
					for (; position < gap[gapID].end; ++position) {
						gapClosedSeq[scafID].seq.push_back(4);
					}
					gapID = std::min(gapID + 1, static_cast<long>(gap.size() - 1));
					continue;
				}

				position += gap[gapID].end - gap[gapID].start - (gap[gapID].closingLength + 1);
				++numClosed;
            } else {
                for (; position < gap[gapID].end; ++position) {
					gapClosedSeq[scafID].seq.push_back(4);
                }
            }
            gapID = std::min(gapID + 1, static_cast<long>(gap.size() - 1));
        }
    }

    std::cerr << "TOTAL_NUM_CLOSED_GAPS = " << numClosed << std::endl;
}


void GapClose::Closer::findCircularGapClosedSeq(const long numThread)
{
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(dynamic)
	for (unsigned long i = 0; i < gapClosedSeq.size(); ++i) {
		int overlap = this->selfOverlap(gapClosedSeq[i].seq, MIN_OVERLAP_FOR_CIRCLE);
		if (overlap >= MIN_OVERLAP_FOR_CIRCLE) {
			gapClosedSeq[i].circularFlag = true;
			gapClosedSeq[i].seq.resize(gapClosedSeq[i].seq.size() - overlap);
			gapClosedSeq[i].rightExtendedLen = std::min(gapClosedSeq[i].rightExtendedLen - overlap, 0);
		}
	}
}

long GapClose::Closer::selfOverlap(string &seq, const long minOverlap)
{
	long i;
	for (i = seq.size() - 1; i >= minOverlap; --i) {
		long j;
		for (j = i; j >= 1; --j) {
			if (seq[j - 1] != seq[seq.size() - i + j - 1])
				break;
		}
		if (j < 1)
			return i;
	}

	return 0;
}

void GapClose::Closer::printGapClosedSeq(const string &outPrefix, const unsigned long long readLength, const unsigned long long contigMaxK, const bool extentionFlag)
{
	long seqID = 0;
	long numScaffoldInOneFile = 0;
	long fileIndex = 0;

	unsigned long prefixPosition = scaffoldFileName[fileIndex].find(outPrefix);
	if (prefixPosition != std::string::npos)
		prefixPosition += outPrefix.size();
	else
		prefixPosition = 0;
	if (!isalnum(scaffoldFileName[fileIndex][prefixPosition]))
		++prefixPosition;
	string outName = outPrefix + "_gapClosed_" + scaffoldFileName[fileIndex].substr(prefixPosition);
    std::ofstream out(outName);

	for (auto it = gapClosedSeq.begin(); it != gapClosedSeq.end(); ++it) {
		if (numScaffoldInOneFile >= numScaffoldInFile[fileIndex]) {
			numScaffoldInOneFile = 0;
			++fileIndex;
			out.close();
			unsigned long prefixPosition = scaffoldFileName[fileIndex].find(outPrefix);
			if (prefixPosition != std::string::npos)
				prefixPosition += outPrefix.size();
			else
				prefixPosition = 0;
			if (!isalnum(scaffoldFileName[fileIndex][prefixPosition]))
				++prefixPosition;
			outName = outPrefix + "_gapClosed_" + scaffoldFileName[fileIndex].substr(prefixPosition);
			out.open(outName);

			while (numScaffoldInFile[fileIndex] == 0) {
				++fileIndex;
				out.close();
				unsigned long prefixPosition = scaffoldFileName[fileIndex].find(outPrefix);
				if (prefixPosition != std::string::npos)
					prefixPosition += outPrefix.size();
				else
					prefixPosition = 0;
				if (!isalnum(scaffoldFileName[fileIndex][prefixPosition]))
					++prefixPosition;
				outName = outPrefix + "_gapClosed_" + scaffoldFileName[fileIndex].substr(prefixPosition);
				out.open(outName);
			}
		}

		++seqID;
		++numScaffoldInOneFile;

		if (it->heteroFlag) {
			if (it->redundantFlag && gapClosedSeq[it->heteroCounterPartIndex].redundantFlag)
				continue;
		}
		else {
			if (it->redundantFlag)
				continue;
		}

		unsigned long i = 0;
		for (; it->baseName.size(); ++i) {
			if (isdigit(it->baseName[i])) {
				for (;it->baseName.size(); ++i) {
					if (!isdigit(it->baseName[i]))
						break;
				}
				break;
			}
		}
		if (it->baseName.size() > 0)
			out << '>' << it->baseName.substr(0, i);
		else
			out << ">scaffold" << seqID;

		unsigned long outputSeqLength = it->seq.size();
		if (!extentionFlag)
			outputSeqLength -= (it->leftExtendedLen + it->rightExtendedLen);
		out << "_len" << outputSeqLength << "_cov" << it->coverage << "_read" << readLength << "_maxK" << contigMaxK;

		if (it->circularFlag)
			out << "_circular_candidate";
//		out << " pre_name:" << it->baseName;
		out << "\n";

		if (extentionFlag) {
			for (unsigned long j = 0; j < it->seq.size(); ++j) {
				out.put(platanus::Bin2Char(it->seq[j]));
				if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
					out.put('\n');
			}
			if (it->seq.size() % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
				out.put('\n');
		}
		else {
			for (unsigned long j = 0; j < it->seq.size() - it->leftExtendedLen - it->rightExtendedLen; ++j) {
				out.put(platanus::Bin2Char(it->seq[j + it->leftExtendedLen]));
				if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
					out.put('\n');
			}
			if ((it->seq.size() - it->leftExtendedLen - it->rightExtendedLen) % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
				out.put('\n');
		}
	}

	if (fileIndex + 1 < static_cast<long>(numScaffoldInFile.size())) {
		out.close();
		++fileIndex;
		unsigned long prefixPosition = scaffoldFileName[fileIndex].find(outPrefix);
		if (prefixPosition != std::string::npos)
			prefixPosition += outPrefix.size();
		else
			prefixPosition = 0;
		if (!isalnum(scaffoldFileName[fileIndex][prefixPosition]))
			++prefixPosition;
		outName = outPrefix + "_gapClosed_" + scaffoldFileName[fileIndex].substr(prefixPosition);
		out.open(outName);
	}
    out.close();
}

//////////////////////////////////////////////////////////////////////////////////////
// calculate max seq length because of decreasing new operator
//////////////////////////////////////////////////////////////////////////////////////
unsigned long long GapClose::Closer::calcMaxGapSeqBase(const long numThread) const
{
    vector<unsigned long long> base(numThread, 0);
    omp_set_num_threads(numThread);
    # pragma omp parallel for schedule(dynamic)
    for (long threadID = 0; threadID < numThread; ++threadID) {
        for (unsigned long long i = threadID; i < gap.size(); i += numThread) {
            unsigned long long tmp = 0;
            for (auto it = gap[i].seq.begin(), end = gap[i].seq.end(); it != end; ++it)
                tmp += it->size() + 1;
            if (tmp > base[threadID]) base[threadID] = tmp;
        }
    }
    return *(std::max_element(base.begin(), base.end()));

}


FILE *GapClose::Closer::localAssemble(const unsigned long long numThread)
{
    FILE *temporaryFP[numThread];
    const unsigned long long size = this->gap.size();
    std::vector<GAP *> buffer(size);

    std::cerr << "assembling localized reads..." << std::endl;

    for (unsigned long long i = 0; i < size; ++i)
        buffer[i] = &(this->gap[i]);

    std::sort(buffer.begin(), buffer.end(), GAP());

    for (unsigned long long i = 0; i < numThread; ++i)
        temporaryFP[i] =  platanus::makeTemporaryFile();


    # pragma omp parallel for schedule(dynamic)
    for (unsigned long long i = 0; i < size; ++i) {
        if (buffer[i]->state != GAP::UNCLOSED || buffer[i]->tooManyReads) continue;
        GapCloseDBG<BRUIJN_MIN_KMER, BRUIJN_MAX_KMER> ele(BRUIJN_MIN_KMER, BRUIJN_MAX_KMER, MIN_COVERAGE, minOverlapDe, maxMissRate, buffer[i]);

        ele.gapAssemble();
        ele.saveContig(temporaryFP[omp_get_thread_num()]);
    }

    for (unsigned long long i = 1; i < numThread; ++i) {
        rewind(temporaryFP[i]);
        char c;
        while (fread(&c, sizeof(char), 1, temporaryFP[i])) {
            putc(c, temporaryFP[0]);
        }
        fclose(temporaryFP[i]);
    }

    return temporaryFP[0];
}

void GapClose::Closer::setFastaFileInfo(const vector<string> fastaFileName)
{
	platanus::setFastaFileNameAndNumber(fastaFileName, scaffoldFileName, numScaffoldInFile);
}

void GapClose::Closer::setHeteroInfoOfGapClosedSeq()
{
	const vector<string> heteroPrefix = {"primary_bubble", "secondary_bubble"};

	unsigned long seqIndex;
	vector<long> seqIndexToID(gapClosedSeq.size(), 0);
	std::unordered_map<long, std::array<long, 2> > bubbleIDToIndex;

	for (seqIndex = 0; seqIndex < gapClosedSeq.size(); ++seqIndex) {
		for (unsigned prefixIndex = 0; prefixIndex < heteroPrefix.size(); ++prefixIndex) {
			if (gapClosedSeq[seqIndex].baseName.size() <= heteroPrefix[prefixIndex].size() || gapClosedSeq[seqIndex].baseName.substr(0, heteroPrefix[prefixIndex].size()) != heteroPrefix[prefixIndex])
				continue;

			unsigned long bubbleIDEnd = heteroPrefix[prefixIndex].size();
			for (; bubbleIDEnd < gapClosedSeq[seqIndex].baseName.size(); ++bubbleIDEnd) {
				if (!isdigit(gapClosedSeq[seqIndex].baseName[bubbleIDEnd]))
					break;
			}
			if (bubbleIDEnd == heteroPrefix[prefixIndex].size())
				continue;

			gapClosedSeq[seqIndex].heteroFlag = true;

			long bubbleID = stoi(gapClosedSeq[seqIndex].baseName.substr(heteroPrefix[prefixIndex].size(), bubbleIDEnd - heteroPrefix[prefixIndex].size()));

			bubbleIDToIndex[bubbleID][prefixIndex] = seqIndex;
			auto mapItr = bubbleIDToIndex.find(bubbleID);
			if (mapItr == bubbleIDToIndex.end()) {
				bubbleIDToIndex[bubbleID] = std::array<long, 2>();
				bubbleIDToIndex[bubbleID][prefixIndex] = seqIndex;
			}
			else {
				mapItr->second[prefixIndex] = seqIndex;
			}
		}
	}

	for (auto mapItr = bubbleIDToIndex.begin(); mapItr != bubbleIDToIndex.end(); ++mapItr) {
		for (unsigned prefixIndex = 0; prefixIndex < heteroPrefix.size(); ++prefixIndex) {
			seqIndex = mapItr->second[prefixIndex];
			if (gapClosedSeq[seqIndex].heteroFlag)
				gapClosedSeq[seqIndex].heteroCounterPartIndex = mapItr->second[prefixIndex^1];
		}
	}
}

void GapClose::Closer::markRedundantSeq(const long maxPrefixLength, const unsigned long long numThread)
{
	unsigned long prefixLength = maxPrefixLength;

	for (unsigned long seqIndex = 0; seqIndex < gapClosedSeq.size(); ++seqIndex) {
		if (prefixLength > gapClosedSeq[seqIndex].seq.size())
			prefixLength = gapClosedSeq[seqIndex].seq.size();
	}

	std::unordered_map<string, vector<long> > prefixToSeqIndex;

	for (unsigned long seqIndex = 0; seqIndex < gapClosedSeq.size(); ++seqIndex) {
		auto mapItr = prefixToSeqIndex.find(gapClosedSeq[seqIndex].seq.substr(0, prefixLength));
		if (mapItr == prefixToSeqIndex.end())
			prefixToSeqIndex[gapClosedSeq[seqIndex].seq.substr(0, prefixLength)] = std::vector<long>(1, seqIndex);
		else
		 	mapItr->second.push_back(seqIndex);
	}


    omp_set_num_threads(numThread);
	vector<vector<char> > redundantFlag(numThread);

	# pragma omp parallel for schedule(static, 1)
	for (unsigned threadID = 0; threadID < numThread; ++threadID) {
		redundantFlag[threadID].assign(gapClosedSeq.size(), false);

		for (unsigned long seqIndex = threadID; seqIndex < gapClosedSeq.size(); seqIndex += numThread) {
			for (unsigned strandIndex = 0; strandIndex < 2; ++strandIndex) {
				string targetSeq;
				if (strandIndex == 0) {
					targetSeq = gapClosedSeq[seqIndex].seq;
				}
				else {
					targetSeq.resize(gapClosedSeq[seqIndex].seq.size());
					for (unsigned long baseIndex = 0; baseIndex < gapClosedSeq[seqIndex].seq.size(); ++baseIndex) {
						if (gapClosedSeq[seqIndex].seq[baseIndex] < 4)
							targetSeq[targetSeq.size() - baseIndex - 1] = 3^(gapClosedSeq[seqIndex].seq[baseIndex]);
						else
							targetSeq[targetSeq.size() - baseIndex - 1] = 4;
					}
				}

				for (unsigned long baseIndex = 0; baseIndex < targetSeq.size() - prefixLength + 1; ++baseIndex) {
					auto mapItr = prefixToSeqIndex.find(targetSeq.substr(baseIndex, prefixLength));
					if (mapItr == prefixToSeqIndex.end())
						continue;

					for (auto vecItr = mapItr->second.begin(); vecItr != mapItr->second.end(); ++vecItr) {
						if (*vecItr != seqIndex &&
							(targetSeq.size() > gapClosedSeq[*vecItr].seq.size() || (targetSeq.size() == gapClosedSeq[*vecItr].seq.size() && seqIndex < *vecItr)) &&
							targetSeq.size() - baseIndex >= gapClosedSeq[*vecItr].seq.size() &&
							std::search(gapClosedSeq[*vecItr].seq.begin()+prefixLength, gapClosedSeq[*vecItr].seq.end(), targetSeq.begin()+baseIndex+prefixLength, targetSeq.begin()+baseIndex+gapClosedSeq[*vecItr].seq.size()) == gapClosedSeq[*vecItr].seq.begin()+prefixLength) {

							redundantFlag[threadID][*vecItr] = true;
						}
					}
				}
			}
		}
	}

	for (unsigned threadID = 0; threadID < numThread; ++threadID) {
		for (unsigned long seqIndex = 0; seqIndex < gapClosedSeq.size(); ++seqIndex)
			if (gapClosedSeq[seqIndex].redundantFlag == false && redundantFlag[threadID][seqIndex])
				gapClosedSeq[seqIndex].redundantFlag = true;
	}
}
