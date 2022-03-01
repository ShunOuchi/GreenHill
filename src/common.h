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

#ifndef PLATANUS_COMMON_H
#define PLATANUS_COMMON_H


#include <memory>
#include <string>
#include <string.h>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <cmath>
#include <unistd.h>
#include <sys/resource.h>
#include "binstr.h"


namespace platanus
{
    using std::string;
    using std::vector;
    typedef unsigned long long u64_t;

    enum FILECOMPRESSION {UNCOMPRESSED, GZIP, BZIP2};
    enum FILETYPE {UNKNOWN, FASTA, FASTQ};
    enum ERROR
    {IO=0, FOPEN, TMP, FORMAT, READ, KMERDIST, MAP, BASE, COVERAGE, DOUBLEHASH, CLINK, CDIR, SCAF, SOLVE, GAP, MERGE};

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// ConstParam Class
// set constatnt parameter
    struct ConstParam
    {
        static const unsigned long long MAX_READ_LEN;
        static const unsigned SCAFFOLD_HASH_OVERLAP;
        static const unsigned OUTPUT_LINE_LENGTH;
        static const unsigned MAX_FILE_NUM;
        static const unsigned MAX_FILE_LEN;
        static const unsigned MAX_THREAD;
        static const unsigned long long DEFAULT_CONTIG_READ_LEN;
        static const std::string VERSION;
        static const double DOUBLE_HASH_MAX_LOAD_FACTOR;
		static const double SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR;
		static const double SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR;
		static const double LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR;
		static const double LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR;
    };
//////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// Error Class
// throwing object in exception
    class ErrorBase
    {
    private:
        const ERROR errorID;
        const std::string message;

    public:
        ErrorBase() = delete;
        ErrorBase(const ERROR err, const std::string &mes): errorID(err), message(mes) {}
        ~ErrorBase() {}

        void showErrorMessage(void) const
        {
            std::cerr << "Error(" << errorID << "): " << message << std::endl;
        }

        ERROR getID(void) const { return errorID; }
    };


    struct IOError : public ErrorBase
    {
        IOError(): ErrorBase(IO, "Error, IO exception!!") {}
        IOError(const std::string &mes): ErrorBase(IO, "Error, IO exception!!\n" + mes) {}
        ~IOError() {}
    };

    struct FILEError : public ErrorBase
    {
        FILEError(): ErrorBase(FOPEN, "Error, File open exception!!") {}
        FILEError(const std::string &mes): ErrorBase(FOPEN, "Error, File not open exception!!\n" + mes + " is not open!!") {}
        ~FILEError() {}
    };

    struct TMPError : public ErrorBase
    {
        TMPError(): ErrorBase(TMP, "Error, Temporary file exception!!\nCannot open temporary file.") {}
        TMPError(const std::string &mes): ErrorBase(TMP, "Error, Temporary file exception!!\n" + mes) {}
        ~TMPError() {}
    };

    struct FormatError : public ErrorBase
    {
        FormatError(): ErrorBase(FORMAT, "Error, File fomat exception!!\nRead file is unknown format.") {}
        FormatError(const std::string &mes): ErrorBase(TMP, "Error, File format exception!!\n" + mes) {}
        ~FormatError() {}
    };

    struct ReadError : public ErrorBase
    {
        ReadError(): ErrorBase(READ, "Error, Read file exception!!\nDNA read seq length is too long.") {}
        ReadError(const std::string &mes): ErrorBase(TMP, "Error, Read file exception!!\n" + mes) {}
        ~ReadError() {}
    };

    struct KmerDistError : public ErrorBase
    {
        KmerDistError(): ErrorBase(KMERDIST, "Error, Kmer distribution exception!!\nKmer distribution cannot be caluculated correctly.") {}
        KmerDistError(const std::string &mes): ErrorBase(TMP, "Error, Kmer distribution exception!!\n" + mes) {}
        ~KmerDistError() {}
    };

    struct MapError : public ErrorBase
    {
        MapError() = delete;
        MapError(const std::string &mes): ErrorBase(MAP, "Error, Kmer mapping exception!!\n" + mes) {}
        ~MapError() {}
    };

    struct BaseError : public ErrorBase
    {
        BaseError(): ErrorBase(KMERDIST, "Error, Seq base exception!!\nfind unknwon base (without \"ATGCN\").") {}
        BaseError(const std::string &mes): ErrorBase(MAP, "Error, Seq base exception!!\n" + mes) {}
        ~BaseError() {}
    };

    struct CoverageError : public ErrorBase
    {
        CoverageError(): ErrorBase(COVERAGE, "Error, Coverage exception!!\ncontig/scaffold fasta file's header doesn't contain \"_cov\" word.") {}
        ~CoverageError() {}
    };

    struct DoubleHashError : public ErrorBase
    {
        DoubleHashError(): ErrorBase(DOUBLEHASH, "Error, DoubleHash exception!!\ndouble hash's space should be 2^n.") {}
        ~DoubleHashError() {}
    };

    struct CreateLinkError : public ErrorBase
    {
        CreateLinkError(): ErrorBase(CLINK, "Error, Create link exception!!\nln, cp, mv or cat command failed.") {}
        ~CreateLinkError() {}
    };

    struct CreateDirError : public ErrorBase
    {
        CreateDirError(const std::string &mes): ErrorBase(CDIR, "Error, Create directory exception!!\n" + mes + "is created") {}
        ~CreateDirError() {}
    };

    struct ScaffoldError : public ErrorBase
    {
        ScaffoldError(): ErrorBase(SCAF, "Error, Scaffold exception!!\nplatanus_allee scaffold command failed.") {}
        ~ScaffoldError() {}
    };

    struct SolveDBGError : public ErrorBase
    {
        SolveDBGError(): ErrorBase(SOLVE, "Error, SolveDBG exception!!\nplatanus_allee solve_DBG command failed.") {}
        ~SolveDBGError() {}
    };

    struct PolishError : public ErrorBase
    {
        PolishError(): ErrorBase(SCAF, "Error, Polish exception!!\nplatanus_allee polish command failed.") {}
        ~PolishError() {}
    };

    struct GapError : public ErrorBase
    {
        GapError(): ErrorBase(GAP, "Error, Gap_close exception!!\nplatanus_allee gap_close command failed.") {}
        ~GapError() {}
    };

    struct MergeError : public ErrorBase
    {
        MergeError(): ErrorBase(MERGE, "Error, Merge exception!!\nplatanus_allee merge command failed.") {}
        ~MergeError() {}
    };

    struct AlignerError : public ErrorBase
    {
        AlignerError(): ErrorBase(MERGE, "Error, Aligner exception!!\naligner (minimap2 or minialign) returned an error.") {}
        ~AlignerError() {}
    };
//////////////////////////////////////////////////////////////////////////////////////




    //////////////////////////////////////////////////////////////////////////////////////
    // Convert from DNA seq (A,T,G,C) to binary code and reverse convert
    //////////////////////////////////////////////////////////////////////////////////////
    inline unsigned char Char2Bin(const char c) { return ".\x0.\x1\x3..\x2......\x4"[(c&0xF)]; }
    inline char Bin2Char(const unsigned char b) { return "ACGTN"[(int)b]; }

    inline void MemoryAlert(void)
    {
      std::cerr << "WARNING:: Sorry, memory exceeds specified value!!" << std::endl;
    }




    //////////////////////////////////////////////////////////////////////////////////////
    // make temporary file
    //////////////////////////////////////////////////////////////////////////////////////
	extern std::string globalTmpFileDir;
	inline void setGlobalTmpFileDir(std::string name)
	{
		globalTmpFileDir = name;
	}

    inline FILE *makeTemporaryFile(void)
    {
        int fd;
        FILE *fp;

		char templat[globalTmpFileDir.size() + 7 + 1];
		strcpy(templat, globalTmpFileDir.c_str());
		strcat(templat,"/XXXXXX");

        fd = mkstemp(templat);
        if (fd == -1) {
            throw TMPError();
        }

        fp = fdopen(fd, "wb+");
        unlink(templat);
        return fp;
    }




    //////////////////////////////////////////////////////////////////////////////////////
    // Get de brujin graph connectining condition functions
    // whether the node only one node connected, multiple node connected, or no connection
    //////////////////////////////////////////////////////////////////////////////////////
    inline unsigned long long Flag2Base(const unsigned char flag)
    {
        static const unsigned long long array[] = {4,0,1,4,2,4,4,4,3,4,4,4,4,4,4,4};
        return array[flag];
    }
    inline unsigned long long FlagCount(const unsigned char flag)
    {
        return (flag & 1) + ((flag >> 1) & 1) + ((flag >> 2) & 1) + ((flag >> 3) & 1);
    }


    inline void printContig(const std::string &outputFilename, FILE *contigFP, const double coverageRatio, const double averageLength, const long maxKmerLength, const char *seqNamePrefix, const bool appendFlag=false)
    {
        std::ofstream out;

		if (appendFlag)
			out.open(outputFilename, std::ios::app);
		else
			out.open(outputFilename, std::ios::out);

        unsigned long long length;
        unsigned short cov;
        unsigned long long numSeq = 0;

        rewind(contigFP);

        while (fread(&length, sizeof(unsigned long long), 1, contigFP)) {
            fread(&cov, sizeof(unsigned short), 1, contigFP);
            binstr_t seq(length);
            seq.readTemporaryFile(contigFP);
            if (cov == UINT16_MAX) continue;
            ++numSeq;
            out << ">" << seqNamePrefix << numSeq << "_len" << length << "_cov" << static_cast<unsigned short>(cov * coverageRatio + 0.5) << "_read" << averageLength << "_maxK" << maxKmerLength << "\n";
            unsigned long long j = 0;
            for (; j < length; ++j) {
                out.put(platanus::Bin2Char(seq.get(length - 1 - j)));
                if ((j + 1) % platanus::ConstParam::OUTPUT_LINE_LENGTH == 0)
                    out.put('\n');
            }
            if (j % platanus::ConstParam::OUTPUT_LINE_LENGTH != 0)
                out.put('\n');
        }
        out.close();
    }

	void setLimitNumFileOpenMax();
	void setFastaFileNameAndNumber(const vector<string> fastaFileName, std::vector<std::string> &scaffoldFileName, std::vector<long> &numScaffoldInFile);
	FILECOMPRESSION checkFileCompression(const std::string &filename);
	FILE *openFileAllowingCompression(const std::string &str, const char *mode);
	int closeFileAllowingCompression(FILE *fp, const std::string &filename);
	FILE *getlineFILE(std::string &str, FILE *fp);
	void printNumOpenFiles();


//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// inner class Position
// this class has conting (scaffold) ID and position which read is mapped
    struct Position
    {
        int offset;
        int id;

        Position(): offset(0), id(0) {}
        Position(const int pos, const int i): offset(pos), id(i) {}
        Position(const Position &pos) = default;
        Position &operator=(const Position &pos) = default;
        ~Position() {}

        bool operator<(const Position &a) const
        {
            if (id != a.id)
                return id - a.id < 0;
            else
                return offset - a.offset < 0;
        }

        bool operator==(const Position &a) const
        {
			return (id == a.id && offset == a.offset);
        }

        operator bool() const
        {
            return id != 0;
        }

        void swap(Position &a)
        {
            std::swap(offset, a.offset);
            std::swap(id, a.id);
        }
    };
//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// inner class SEQ
// this class has read information from file.
    struct SEQ
    {
        long length;
        string base;
        int numUnknown;
        vector<int> positionUnknown;

        SEQ(): length(0), base(), numUnknown(0), positionUnknown() {}
        SEQ(const SEQ &seq) = default;
        SEQ &operator=(const SEQ &seq) = default;
        ~SEQ() {};


        bool operator<(const SEQ &seq) const { return length < seq.length; }

        void reverse(void)
        {
            char reverse[length];
            for (int i = 0; i < length; ++i)
                reverse[i] = base[i];
            for (int i = 0; i < length; ++i) {
                base[i] = reverse[length - i - 1] != 4 ? reverse[length - i - 1] ^ 0x3 : 4;
            }
        }

        inline void writeTemporaryFile(FILE *seqFile) const
        {
            fwrite(&numUnknown, sizeof(int), 1, seqFile);
            fwrite(&positionUnknown[0], sizeof(int), numUnknown, seqFile);
            fwrite(&length, sizeof(int), 1, seqFile);
            fwrite(&base[0], sizeof(char), length, seqFile);
        }


        inline int readTemporaryFile(FILE *seqFile)
        {
            if (fread(&numUnknown, sizeof(int), 1, seqFile) != 1)
                return 0;
            positionUnknown.resize(numUnknown);
            fread(&positionUnknown[0], sizeof(int), numUnknown, seqFile);

            fread(&length, sizeof(int), 1, seqFile);
            base.resize(length);
            fread(&base[0], sizeof(char), length, seqFile);
            return 1;
        }


        void put(const string &line)
        {
            int tmpLength = line.length();
            length += tmpLength;
            for (int i = 0; i < tmpLength; ++i) {
                base += Char2Bin((char)line[i]);
            }
        }

        void convertFromString(const string &str)
        {
            length = str.length();
            base.resize(length);
            numUnknown = 0;
            if (str.length() >= ConstParam::MAX_READ_LEN) {
                throw ReadError();
            }
            for (unsigned i = 0; i < length; ++i) {
                if (Char2Bin(str[i]) == 4) {
                    positionUnknown.resize(numUnknown + 1);
                    positionUnknown[numUnknown] = i;
                    ++numUnknown;
                } else {
                    base[i] = Char2Bin(str[i]);
                }
            }
        }

        void show(void) const
        {
            for (unsigned l = 0; l < length; ++l)
                std::cout << platanus::Bin2Char(base[l]);
                std::cout << std::endl;
        }

        //added by ouchi
        SEQ divide(const int pos)
        {
            SEQ tmp;
            if (pos < 0 || pos >= length)
                return tmp;
            tmp.length = length - pos;
            length = pos + 1;
            tmp.numUnknown = 0;
            for (int i = 0; i < numUnknown; ++i) {
                if (positionUnknown[i] > pos) {
                    tmp.positionUnknown.resize(tmp.numUnknown + 1);
                    tmp.positionUnknown[numUnknown] = i;
                    ++tmp.numUnknown;
                }
            }
            numUnknown += -tmp.numUnknown;
            positionUnknown.resize(numUnknown);
            tmp.base = base.substr(pos, tmp.length);
            base.erase(base.begin() + pos, base.end());
            return tmp;
        }
        //added by ouchi

    };
//////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// inner class Contig
// this class has contig information from contig file.
    struct Contig
    {
        FILE *seqFP;
        FILE *coverageFP;
        vector<SEQ> seq;
        vector<std::string> name;
		std::unordered_map<std::string, unsigned> nameIndex;
        unsigned numSeq;
        vector<char> seqPool;
        vector<unsigned short> coverage;

        Contig(): seqFP(NULL), coverageFP(NULL), seq(), name(0), numSeq(0), seqPool(), coverage() {}
        Contig(const Contig &) = delete;
        Contig &operator=(const Contig &) = delete;
        ~Contig() {}


        void clear(void)
        {
            seq.clear();
            seqPool.clear();
            coverage.clear();
            numSeq = 0; //added by ouchi
            name.clear(); //added by ouchi
            nameIndex.clear(); //added by ouchi
        }

        enum READMODE {NORMAL=0, CUTN, ADD_N_EDGE};

        void readFastaCoverageAddNEdge(const string &filename)
        { readFasta(filename, "cov", true, READMODE::ADD_N_EDGE); }
        void readFastaCoverageCutN(const string &filename, const unsigned long long lengthThreshold)
        { readFasta(filename, "cov", true, READMODE::CUTN, lengthThreshold); }
        void readFastaCoverage(const string &filename)
        { readFasta(filename, "cov", true, READMODE::NORMAL); }

        void setSeq(const std::string &contig, const std::string &seqName, const unsigned short cov, const unsigned long long lengthThreshold)
        {
            long size = contig.length();
            SEQ sseq;
            sseq.length = size;
            sseq.base.resize(size);
            for (int i = 0; i < size; ++i) {
                sseq.base[i] = Char2Bin(contig[i]);
            }
            seq.emplace_back(std::move(sseq));

			std::istringstream ss(seqName);
			string newSeqName;
			ss >> newSeqName;
			name.push_back(newSeqName);
            coverage.emplace_back(cov);
            ++numSeq;
        }

        void setSeqAndCutNNNN(const std::string &contig, const std::string &seqName, const unsigned short cov, const unsigned long long lengthThreshold)
        {
            long size = contig.length();
            SEQ sseq;
            if (size < static_cast<long>(lengthThreshold)) return;
            for (int i = 0; i < size; ++i) {
                if (contig[i] != 'N') {
                    sseq.base.push_back(Char2Bin(contig[i]));
                    ++sseq.length;
                } else {
                    if (static_cast<unsigned long long>(sseq.length) > 0) {
                        sseq.base.resize(sseq.length);
                        seq.push_back(sseq);
                        coverage.push_back(cov);
                        ++numSeq;
                    }
                    sseq.length = 0;
                    sseq.base.clear();
                }
            }
            if (static_cast<unsigned long long>(sseq.length) > 0) {
                sseq.base.resize(sseq.length);
                seq.emplace_back(std::move(sseq));
				name.push_back(seqName);
                coverage.emplace_back(cov);
                ++numSeq;
            }
        }


        void setSeqAndAddN(const std::string &contig, const std::string &seqName, const unsigned short cov, const unsigned long long lengthThreshold)
        {
            long size = contig.length() + 2;
            SEQ sseq;
            sseq.length = size;
            sseq.base.resize(size);
            sseq.base[0] = sseq.base[size-1] = 4;
            for (unsigned long long i = 0; i < contig.length(); ++i) {
                sseq.base[i + 1] = Char2Bin(contig[i]);
            }
            seq.emplace_back(std::move(sseq));
			name.push_back(seqName);
            coverage.emplace_back(cov);
            ++numSeq;
        }

        void readFasta(const string &filename, const string &covHeader="cov", const bool isReadCoverage=false, const READMODE mode=NORMAL, const unsigned long long lengthThreshold=0)
        {
			FILE* fp = platanus::openFileAllowingCompression(filename, "r");

            string oneLine, contig, seqName;
            u64_t readCoverage = 0;
            unsigned short readCoverage16Bit = 0;
            void (Contig::*function)(const std::string &contig, const std::string &seqName, const unsigned short coverage, const unsigned long long lengthThreshold) = NULL;
            switch (mode) {
            case NORMAL:
                function = &Contig::setSeq;
                break;
            case CUTN:
                function = &Contig::setSeqAndCutNNNN;
                break;
            case ADD_N_EDGE:
                function = &Contig::setSeqAndAddN;
                break;
            }
            if (seqFP == NULL && (seqFP = makeTemporaryFile()) == NULL) throw TMPError();
            fseek(seqFP, (long)0, SEEK_END);

            if (isReadCoverage) {
                if (coverageFP == NULL && (coverageFP = makeTemporaryFile()) == NULL) throw TMPError();
                fseek(coverageFP, 0, SEEK_END);
            }

			while (platanus::getlineFILE(oneLine, fp) != NULL) {
                if (oneLine[0] == '>') {
					seqName = oneLine.substr(1);
                    break;
                }
            }
            if (feof(fp)) {
				platanus::closeFileAllowingCompression(fp, filename);
				return;
			}
            if (isReadCoverage) {
                readCoverage = findCoverageFromHeader(oneLine.c_str(), covHeader);
                readCoverage16Bit = static_cast<unsigned short>(std::min(readCoverage, static_cast<unsigned long long >(UINT16_MAX)));
                readCoverage = 0;
            }

			while (platanus::getlineFILE(oneLine, fp) != NULL) {
                if (oneLine[0] != '>') {
                    contig += oneLine;
                } else {
                    (this->*function)(contig, seqName, readCoverage16Bit, lengthThreshold);
                    contig = "";
					seqName = oneLine.substr(1);

                    if (isReadCoverage) {
                        readCoverage = findCoverageFromHeader(oneLine.c_str(), covHeader);
                        readCoverage16Bit = std::min((unsigned short)readCoverage, (unsigned short)UINT16_MAX);
                        readCoverage = 0;
                    }


                }
            }
            (this->*function)(contig, seqName, readCoverage16Bit, lengthThreshold);
            contig = "";

			platanus::closeFileAllowingCompression(fp, filename);
        }

        static u64_t findCoverageFromHeader(const char *line, const string &covHeader)
        {
            u64_t cov = 0;
            const char *covHeaderPointer = strstr(line, covHeader.c_str());
            if (covHeaderPointer != NULL) {
				covHeaderPointer += covHeader.length();
				while (*covHeaderPointer != '\0' && isdigit(*covHeaderPointer)) {
					cov = (cov * 10) + (*covHeaderPointer - 48);
					++covHeaderPointer;
				}
			}
			else {
				cov = 1;
			}
			return cov;
        }

        static unsigned long long getReadLengthFromFastaHeader(const std::string &filename)
        {
            return Contig::getReadLength(filename, "read");
        }

        static unsigned long long getMaxKFromFastaHeader(const std::string &filename)
        {
            return Contig::getReadLength(filename, "maxK");
        }

        static unsigned long long getReadLength(const std::string &filename, const std::string &readHeader)
        {
			FILE* fp = platanus::openFileAllowingCompression(filename, "r");
            string oneLine;

			while (platanus::getlineFILE(oneLine, fp) != NULL) {
                if (oneLine[0] == '>') {
                    break;
                }
            }

			platanus::closeFileAllowingCompression(fp, filename);

            unsigned long long readLength = findCoverageFromHeader(oneLine.c_str(), readHeader);
			if (readLength > 1)
				return readLength;
			else
				return 100;
        }

        void loadSeq(void)
        {
            if (seqFP != NULL) {
                rewind(seqFP);
                seq.clear();
                seq.resize(numSeq);
                for (long i = 0; i < numSeq; ++i) {
                    SEQ tmpSEQ;
                    fread(&(tmpSEQ.length), sizeof(long), 1, seqFP);
                    tmpSEQ.base.resize(tmpSEQ.length);
                    fread(&(tmpSEQ.base[0]), sizeof(char), tmpSEQ.length, seqFP);
                    seq[i] = tmpSEQ;
                }
            }
        }


        void loadCoverage(void)
        {
            if (coverageFP == NULL) throw TMPError();
            coverage.resize(numSeq, static_cast<unsigned short>(0));
            rewind(coverageFP);
            for (unsigned i = 0; i < numSeq; ++i) {
                fread(&(coverage[i]), sizeof(unsigned short), 1, coverageFP);
            }
        }

        double calculateAverageCoverage(const long minLength) const
        {
            long sum = 0;
            long num = 0;
            for (long i = 0; i < numSeq; ++i) {
                if (seq[i].length >= minLength) {
                sum += coverage[i] * seq[i].length;
                num += seq[i].length;
                }
            }
            return static_cast<double>(sum) / num;
        }


        double calculateAverageCoverageExcludingOutlier(const long minLength) const
        {
			const double EXCLUSION_FACTOR = 100.0;
            long sum = 0;
            long num = 0;
            for (long i = 0; i < numSeq; ++i) {
                if (seq[i].length >= minLength) {
					sum += coverage[i] * seq[i].length;
					num += seq[i].length;
                }
            }
            double temporaryMean = static_cast<double>(sum) / num;

			sum = num = 0;
            for (long i = 0; i < numSeq; ++i) {
                if (seq[i].length >= minLength) {
					if (coverage[i] >= temporaryMean / EXCLUSION_FACTOR && coverage[i] <= temporaryMean * EXCLUSION_FACTOR) {
						sum += coverage[i] * seq[i].length;
						num += seq[i].length;
					}
                }
            }

            return static_cast<double>(sum) / num;
        }


		void readContigFP(FILE *contigFP)
		{
			unsigned long long length;
			unsigned short cov;

			seq.clear();
			coverage.clear();
			numSeq = 0;

			rewind(contigFP);
			while (fread(&length, sizeof(unsigned long long), 1, contigFP)) {
				fread(&cov, sizeof(unsigned short), 1, contigFP);
				binstr_t bseq(length);
				bseq.readTemporaryFile(contigFP);

				++numSeq;

				seq.resize(numSeq);
				seq.back().length = length;
				seq.back().base.resize(length);
				for (unsigned long long j = 0; j < length; ++j)
					seq.back().base[j] = bseq.get(length - 1 - j);

				coverage.resize(numSeq);
				coverage.back() = cov;
			}
		}


        long getSeqLengthMedian(void) const
        {
            std::vector<SEQ> tmp = seq;
            sort(tmp.begin(), tmp.end());
            return tmp[tmp.size()/2].length;

        }

        void setNameIndex(void)
        {
			for (unsigned i = 0; i < name.size(); ++i) {
				nameIndex[name[i]] = i;
			}
        }

        void printNameIndexCSV(const std::string &outputFilename)
        {
			std::ofstream out(outputFilename.c_str());

			out << "Index,Name" << std::endl;
			for (unsigned i = 0; i < name.size(); ++i) {
				out << i + 1 << ',' << name[i] << std::endl;
			}
			out.close();
        }

    };
//////////////////////////////////////////////////////////////////////////////////////
// Normal distribution class for calculation
    template <typename T>
    class NormalDistribution
    {
    private:
        T average;
        double sd;
        double k1;
        double k2;
    public:
        NormalDistribution(T ave, double sdv)
        : average(ave), sd(sdv), k1(static_cast<double>(1) / (sqrt(M_PI) * 2 * sd)), k2(-2 * sd * sd) {}
        double calc(const T val) { return k1 * exp(static_cast<double>((val - average) * (val - average)) / k2); }
        double operator() (const T val) { return calc(val); }
    };
//////////////////////////////////////////////////////////////////////////////////////
// PairEqual and Pairhash class
// for unordered map keys
    struct PairEqual{
        bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const
        {
            return (a.first == b.first) && (a.second == b.second);
        }
    };

    struct PairHash{
        size_t operator()(const std::pair<int, int> &a) const
        {
            return (size_t)(a.first * a.second);
        }
    };

    struct PairFirstLess{
        bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const
        {
            return (a.first < b.first);
        }
    };

    struct PairFirstGreater{
        bool operator()(const std::pair<int, int> &a, const std::pair<int, int> &b) const
        {
            return (a.first > b.first);
        }
    };
//////////////////////////////////////////////////////////////////////////////////////

}





#endif
