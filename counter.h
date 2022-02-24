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

#ifndef COUNTER_H
#define COUNTER_H

#include "common.h"
#include "kmer.h"
#include "gapClose.h"
#include "doubleHash.h"
#include <vector>
#include <unordered_map>
#include <cmath>
#include <memory>

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// kmer occurrence counter class
template <typename KMER>
class Counter
{
    typedef unsigned long long u64_t;
private:
    // constant parameter
    static const unsigned HASH_DIVID;
    static const unsigned long long HASH_SHIFT;


    std::vector<u64_t> lengthDistribution;
    std::vector<u64_t> occurrenceDistribution;
    DoubleHash<typename KMER::keyType, unsigned short> occurrenceTable;
    std::unordered_map<typename KMER::keyType, unsigned short, typename KMER::hasher> occurrenceGapTable;
    u64_t kmerLength;
    u64_t maxOccurrence;



    // swap memory delete function
    void deleteLengthDistribution(void)
    {std::vector<u64_t>().swap(lengthDistribution); }
    void deleteOccurrenceTable(void)
    {DoubleHash<typename KMER::keyType, unsigned short>().swap(occurrenceTable); }
    void deleteOccurrenceDistribution(void)
    {std::vector<unsigned long long>().swap(occurrenceDistribution); }



    double calcDistributionAverage(const std::vector<u64_t> &distribution, const u64_t start, const u64_t end) const;
    void countKmerPerThreadFirst(const unsigned long long kLength, bool &loopFlag, DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *readFP, FILE *unmappedFP, omp_lock_t lock[], const unsigned long long numThread);
    void countKmerPerThreadSecond(bool &loopFlag, DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *readFP, FILE *unmappedFP, omp_lock_t lock[], const unsigned iterateTimes, const unsigned long long numThread);
    void countKmerOrWriteTemporary(bool &loopFlag, const typename KMER::keyType &key, DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *unmappedFP, omp_lock_t lock[], const KMER &kmer, const unsigned iterateTimes=32);
    void divideKmerUsedMakingPreviousContig(const unsigned kmerLength, KMER &kmer, typename KMER::keyType &key, bool &loopFlag, FILE *readFP, FILE *unmappedFP);
    unsigned writeKmerDistribution(DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *kmerFP);

    void countKmerPerThreadSecond(bool &loopFlag, FILE *readFP, FILE *unmappedFP, omp_lock_t lock[], const unsigned iterateTimes, const unsigned long long numThread);
    void countKmerOrWriteTemporary(bool &loopFlag, const typename KMER::keyType &key, FILE *unmappedFP, omp_lock_t lock[], const KMER &kmer, const unsigned iterateTimes);
    unsigned writeKmerDistribution(FILE *kmerFP);

public:
    FILE *kmerFP;

    // constructor and destructor
    Counter(): lengthDistribution(), occurrenceDistribution(), occurrenceTable(), occurrenceGapTable(), kmerLength(0), maxOccurrence(0), kmerFP(NULL){}
    ~Counter() {}


    // getter, setter and similar function
    u64_t getMaxOccurrence(void) const {return maxOccurrence; }
    u64_t getLengthDistributionI(const u64_t position) const {return lengthDistribution[position]; }
    u64_t getKmerLength(void) const {return kmerLength; }
    u64_t getTableSize(void) const {return occurrenceTable.size(); }
    void setOccurrenceValue(const typename KMER::keyType &key, const unsigned short value) {occurrenceTable[key] = value; }
    void setKmerLength(const unsigned long long kmer) {kmerLength = kmer; }
    void setOccurrenceTableSize(const unsigned long long size) {occurrenceTable.resize(size); }

    unsigned short findValue(const typename KMER::keyType &key)
    {
        return occurrenceTable.find_any(key)->second;
    }

    unsigned short findValueAndOverwrite(const typename KMER::keyType &key, const unsigned short overwrite)
    {
        auto tableIterator = occurrenceTable.find_any(key);
        short tmp = tableIterator->second;
        tableIterator->second = overwrite;
        return tmp;
    }

    std::pair<typename KMER::keyType, unsigned short> *getOccurrenceIterator(const typename KMER::keyType &key)
    {
        return occurrenceTable.find_any(key);
    }

    std::pair<typename KMER::keyType, unsigned short> *getOccurrenceEnd(void)
    {
        return occurrenceTable.end();
    }

    // clear function
    void deleteAllTable(void)
    {
        deleteLengthDistribution();
        deleteOccurrenceTable();
        deleteOccurrenceDistribution();
    }
    void cleanAllTable(void)
    {
        deleteLengthDistribution();
        occurrenceTable.clean();
        deleteOccurrenceDistribution();
    }

    void clearAllTable(void)
    {
        lengthDistribution.clear();
        occurrenceTable.clear();
        occurrenceDistribution.clear();
    }

    void clearOccurrenceTable(void)
    {
        occurrenceTable.clear();
    }



    // calculate distribution average
    double calcLengthDistributionAverage(const u64_t start, const u64_t end) const
    {return calcDistributionAverage(lengthDistribution, start, end); }

    double calcOccurrenceDistributionAverage(const u64_t start, const u64_t end) const
    {return calcDistributionAverage(occurrenceDistribution, start, end); }


    unsigned long long getLeftLocalMinimalValue(const unsigned long long windowSize) const;
    unsigned long long makeKmerReadDistributionMT(const unsigned long long  kLength, FILE **readFP, const unsigned long long memory, const unsigned long long numThread);
    unsigned long long loadKmer(const unsigned long long minOccurrence, const unsigned long long doubleHashSize);
    void swapOccurrenceTable(const unsigned k, DoubleHash<typename KMER::keyType, unsigned short> &table);
    unsigned long long makeKmerReadDistributionConsideringPreviousGraph(const unsigned k, FILE **readFP, const unsigned long long memory, const unsigned long long numThread);
    unsigned long long makeKmerReadDistributionFromContig(platanus::Contig &contig, const unsigned long long minContigLength, const unsigned long long minOccurrence, const unsigned long long memory);
    void pickupReadMatchedEdgeKmer(FILE **readFP, const bool readFPCloseFlag);
    FILE *sortedKeyFromKmerFile(const unsigned long long minOccurrence);
    void outputOccurrenceDistribution(const std::string &filename);
    void countKmerForGapClose(const std::vector<std::vector<char> > &seq);


    // using only gap close
    void assignOccurrenceKeysGap(std::vector<typename KMER::keyType> &sorter) const
    {
        for (auto it = occurrenceGapTable.begin(), end = occurrenceGapTable.end(); it != end; ++it) {
            if (it->second != 0)
            sorter.push_back(it->first);
        }
        std::sort(sorter.begin(), sorter.end());
    }

    unsigned short findValueGap(const typename KMER::keyType &key)
    {
        auto tableIterator = occurrenceGapTable.find(key);
        if (tableIterator != occurrenceGapTable.end()) {
            return tableIterator->second;
        } else {
            return 0;
        }
    }

    unsigned short findValueAndOverwriteGap(const typename KMER::keyType &key, const unsigned short overwrite)
    {
        auto tableIterator = occurrenceGapTable.find(key);
        short tmp = tableIterator->second;
        tableIterator->second = overwrite;
        return tmp;
    }
    void setOccurrenceValueGap(const typename KMER::keyType &key, const unsigned short value) {occurrenceGapTable[key] = value; }
    void clearOccurrenceTableGap(void)
    {
        occurrenceGapTable.clear();
    }


};
//////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
const unsigned long long Counter<KMER>::HASH_SHIFT = 10;
template <typename KMER>
const unsigned Counter<KMER>::HASH_DIVID = 1 << HASH_SHIFT;




//////////////////////////////////////////////////////////////////////////////////////
// calculate distribution average
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
double Counter<KMER>::calcDistributionAverage(const std::vector<unsigned long long> &distribution, const unsigned long long start, const unsigned long long end) const
{
    unsigned long long sumDistribution = 0;
    unsigned long long numDistribution = 0;
    if (end > distribution.size() || start > end) {
        throw platanus::KmerDistError();
    }

    for (unsigned long long i = start; i <= end; ++i) {
        sumDistribution += i * distribution[i];
        numDistribution += distribution[i];
    }
    if (numDistribution != 0)
        return static_cast<double>(sumDistribution) / numDistribution;
    else {
        throw platanus::KmerDistError();
    }
}


//////////////////////////////////////////////////////////////////////////////////////
// calculate kmer distribution local minimal value
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long Counter<KMER>::getLeftLocalMinimalValue(const unsigned long long windowSize) const
{
    unsigned long long i;
    unsigned long long windowVectorSize = maxOccurrence - windowSize + 2;
    std::vector<unsigned long long> window;
    if (maxOccurrence <= windowSize)
        return 0;

    window.resize(windowVectorSize, 0);

    for (i = 0; i < windowSize; ++i)
        window[1] += occurrenceDistribution[1 + i];

    for (i = 2; i < windowVectorSize; ++i) {
        window[i] = window[i - 1] - occurrenceDistribution[i - 1] + occurrenceDistribution[i + windowSize - 1];
        if (window[i] >= window[i - 1])
            break;
    }
    if (i <= maxOccurrence)
        return (i - 1 + windowSize / 2);
    else
        return (1 + windowSize / 2);
}




//////////////////////////////////////////////////////////////////////////////////////
// read file and count kmer occurrence in multi thread
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long Counter<KMER>::makeKmerReadDistributionMT(const unsigned long long  kLength, FILE **readFP, const unsigned long long memory, const unsigned long long numThread)
{

    std::vector<Counter> counterMT(numThread);
    DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[HASH_DIVID];
    omp_lock_t lock[HASH_DIVID];
    FILE *unmappedFP[numThread];
    bool loopFlag[numThread];
    FILE *tmpFP[numThread];
    this->kmerLength = kLength;

    omp_set_num_threads(numThread);

    for (unsigned long long i = 0; i < numThread; ++i) {
        unmappedFP[i] = platanus::makeTemporaryFile();
        counterMT[i].deleteAllTable();
        loopFlag[i] = false;
    }
    if (kmerFP != NULL)
        fclose(kmerFP);
    kmerFP = platanus::makeTemporaryFile();


    // if DoubleHash is filled keys less than half (considering MAX_LOAD_FACTOR), add find iterate times
    long base = sizeof(std::pair<typename KMER::keyType, unsigned short>);
    unsigned long long tmpMemory = memory / base;
    base = log(tmpMemory) / log(2);
    tmpMemory = std::pow(2, base);
    if (tmpMemory > memory) {
        while (tmpMemory > memory) {
            tmpMemory >>= 1;
        }
    }
    const unsigned long long doubleHashSize = tmpMemory;
    tmpMemory /= HASH_DIVID;
    for (unsigned long long i = 0; i < HASH_DIVID; ++i) {
        tmpOccurrenceTable[i].resize(tmpMemory);
        omp_init_lock(&lock[i]);
    }

    if (kmerFP != NULL)
        fclose(kmerFP);
    kmerFP = platanus::makeTemporaryFile();


    // count kmer
    # pragma omp parallel for schedule(static, 1)
    for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
        counterMT[threadID].countKmerPerThreadFirst(kLength, loopFlag[threadID], tmpOccurrenceTable, readFP[threadID], unmappedFP[threadID], lock, numThread);
    }

    // collect read length distribution
    lengthDistribution.resize(platanus::ConstParam::MAX_READ_LEN + 1, 0);
    for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
        for (unsigned long long i = 0; i < platanus::ConstParam::MAX_READ_LEN + 1; ++i)
            lengthDistribution[i] += counterMT[threadID].getLengthDistributionI(i);
        counterMT[threadID].deleteLengthDistribution();
    }

    occurrenceDistribution.clear();
    occurrenceDistribution.resize(UINT16_MAX);


    // if kmer not be filled occurrenceTable exist, iterate count
    unsigned iterateTimes = 32;
    while (1) {
        iterateTimes += writeKmerDistribution(tmpOccurrenceTable, kmerFP);

        // check whether finish count kmer
        // if loopFlag on, part of kmer is written in temporary file because of used large memory
        unsigned long long i = 0;
        while (i < numThread) {
            if (loopFlag[i]) break;
            ++i;
        }
        if (i == numThread)
            break;

	  // std::clog<<"redo count"<<std::endl;

        #   pragma omp parallel for schedule(static, 1)
        for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
            loopFlag[threadID] = false;
            tmpFP[threadID] = platanus::makeTemporaryFile();
            countKmerPerThreadSecond(loopFlag[threadID], tmpOccurrenceTable, unmappedFP[threadID], tmpFP[threadID], lock, iterateTimes, numThread);
            fclose(unmappedFP[threadID]);
            unmappedFP[threadID] = tmpFP[threadID];
        }
    }

    // delete temporary file
    for (unsigned long long i = 0; i < numThread; ++i)
        fclose(unmappedFP[i]);

    // search max occurrence in kmer distribution
    for (unsigned short i = UINT16_MAX - 1; i > 0; --i) {
        if (occurrenceDistribution[i] > 0) {
            this->maxOccurrence = i;
            break;
        }
    }

    for (unsigned long long i = 0; i < HASH_DIVID; ++i) {
        omp_destroy_lock(&lock[i]);
    }

    return doubleHashSize;
}



//////////////////////////////////////////////////////////////////////////////////////
// count kmer multi thread for initial
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void Counter<KMER>::countKmerPerThreadFirst(const unsigned long long kLength, bool &loopFlag, DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *readFP, FILE *unmappedFP, omp_lock_t lock[], const unsigned long long numThread)
{
    platanus::SEQ seq;
    lengthDistribution.resize(platanus::ConstParam::MAX_READ_LEN + 1, 0);
    kmerLength = kLength;
    typename KMER::keyType key;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));


    rewind(readFP);
    KMER kmer(kmerLength);

    // set kmer count in occurrenceTable.
    // if occurrenceTable is too large, write temporary file unstoredFP.
    while (seq.readTemporaryFile(readFP)) {
        ++lengthDistribution[seq.length];
        if (static_cast<unsigned>(seq.length) < kmerLength)
            continue;
        ++seq.numUnknown;
        seq.positionUnknown.resize(seq.numUnknown);
        seq.positionUnknown[seq.numUnknown - 1] = platanus::ConstParam::MAX_READ_LEN + 1;
        seq.numUnknown = 0;
        for (unsigned i = 0; i < kmerLength - 1; ++i) {
            kmer.setForward(kmerLength - i - 2, seq.base[i]);
            kmer.setReverse(i + 1, 0x3 ^ seq.base[i]);
        }
        for (unsigned i = 0; i < seq.length - kmerLength + 1; ++i) {
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.setForward(0, seq.base[i + kmerLength - 1]);
            kmer.reverse >>= 2;
            kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[i + kmerLength - 1]);
            if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) < i + kmerLength) {
                if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) <= i) {
                    ++seq.numUnknown;
                }
                continue;
            }
            key = std::min(kmer.forward, kmer.reverse);
            countKmerOrWriteTemporary(loopFlag, key, tmpOccurrenceTable, unmappedFP, lock, kmer);
        }
    }

}



//////////////////////////////////////////////////////////////////////////////////////
// count kmer multi thread
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void Counter<KMER>::countKmerPerThreadSecond(bool &loopFlag, DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *readFP, FILE *unmappedFP, omp_lock_t lock[], const unsigned iterateTimes, const unsigned long long numThread)
{
    KMER kmer(kmerLength);
    rewind(readFP);
    while (kmer.readTemporaryFileForward(readFP)) {
        countKmerOrWriteTemporary(loopFlag, kmer.forward, tmpOccurrenceTable, unmappedFP, lock, kmer, iterateTimes);
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// if this kmer(key) already exist, add
// else write temporary or add new key in table
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void Counter<KMER>::countKmerOrWriteTemporary(bool &loopFlag, const typename KMER::keyType &key, DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *unmappedFP, omp_lock_t lock[], const KMER &kmer, const unsigned iterateTimes)
{
    const unsigned mod = key & static_cast<unsigned long long>(HASH_DIVID - 1);
    omp_set_lock(&lock[mod]);
    auto it = tmpOccurrenceTable[mod].find_times_any(key, iterateTimes);
    if (it->second == 0) {
        it->first = key;
        it->second = 1;
    } else if (it->first == key) {
        if (it->second != UINT16_MAX - 1) {
            ++(it->second);
        }
    } else {
            loopFlag = true;
            kmer.writeKey(unmappedFP, key);
    }
    omp_unset_lock(&lock[mod]);
}


//////////////////////////////////////////////////////////////////////////////////////
// write kmer distribution in temporary file
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline unsigned Counter<KMER>::writeKmerDistribution(DoubleHash<typename KMER::keyType, unsigned short> tmpOccurrenceTable[], FILE *kmerFP)
{
    KMER kmer(kmerLength);
    unsigned long long num = 0;
    // write kmer int temporary file
    for (unsigned hashID = 0; hashID < HASH_DIVID; ++hashID) {
        auto tableIterator = tmpOccurrenceTable[hashID].begin();
        auto tableEnd = tmpOccurrenceTable[hashID].end();
        for (; tableIterator != tableEnd; ++tableIterator) {
            if (tableIterator->second != 0) {
                ++num;
                kmer.writeKey(kmerFP, tableIterator->first);
                fwrite(&(tableIterator->second), sizeof(unsigned short), 1, kmerFP);
                ++occurrenceDistribution[tableIterator->second];
                tableIterator->second = 0;
            }
        }
    }

    // if DoubleHash is filled keys less than half (considering MAX_LOAD_FACTOR), add find iterate times
    if (num / HASH_DIVID < tmpOccurrenceTable[0].size() * platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR / 2)
        return 1;
    else 
        return 0;
}


template <typename KMER>
unsigned long long Counter<KMER>::makeKmerReadDistributionFromContig(platanus::Contig &contig, const unsigned long long minContigLength, const unsigned long long minOccurrence, const unsigned long long memory)
{
    KMER kmer(this->kmerLength);
    typename KMER::keyType key;
    const unsigned long long minLength = this->kmerLength;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));

    std::cerr << "K = " << this->kmerLength << ", loading kmers from contigs...\n";

    long base = sizeof(std::pair<typename KMER::keyType, unsigned short>);
    unsigned long long tmpMemory = memory / base;
    base = log(tmpMemory) / log(2);
    tmpMemory = std::pow(2, base);
    if (tmpMemory > memory) {
        while (tmpMemory > memory) {
            tmpMemory >>= 1;
        }
    }

    unsigned long long totalKmer = 0;
    for (unsigned long long i = 0; i < contig.numSeq; ++i) {
        totalKmer += contig.seq[i].length - this->kmerLength + 1;
    }
    unsigned long long size = log(totalKmer / platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR) / log(2);
    size = std::pow(2, size + 1);
    unsigned long long upper = std::max(size, tmpMemory);
    if (upper == size) platanus::MemoryAlert();

    occurrenceTable.resize(upper);

    for (unsigned long long contigID = 0; contigID < contig.numSeq; ++contigID) {
        platanus::SEQ &seq = contig.seq[contigID];
        if (static_cast<unsigned>(seq.length) < minLength)
            continue;
        ++seq.numUnknown;
        seq.positionUnknown.resize(seq.numUnknown);
        seq.positionUnknown[seq.numUnknown - 1] = platanus::ConstParam::MAX_READ_LEN + 1;
        seq.numUnknown = 0;
        for (unsigned i = 0; i < kmerLength - 1; ++i) {
            kmer.setForward(kmerLength - i - 2, seq.base[i]);
            kmer.setReverse(i + 1, 0x3 ^ seq.base[i]);
        }
        for (unsigned i = 0; i < seq.length - kmerLength + 1; ++i) {
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.setForward(0, seq.base[i + kmerLength - 1]);
            kmer.reverse >>= 2;
            kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[i + kmerLength - 1]);
            /*
            if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) < i + kmerLength) {
                if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) <= i) {
                    ++seq.numUnknown;
                }
                continue;
            }
            */
            key = std::min(kmer.forward, kmer.reverse);
            //auto ptr = occurrenceTable.find_any(key);
            unsigned short &ptr = occurrenceTable[key];
            //ptr->first = key;
            //if (ptr->second != 0) ++b;
            //if (ptr->second < std::max(static_cast<unsigned long long >(contig.coverage[contigID]), minOccurrence)) {
            //    ptr->second = std::max(static_cast<unsigned long long >(contig.coverage[contigID]), minOccurrence);
            //}
            if (ptr < std::max(static_cast<unsigned long long >(contig.coverage[contigID]), minOccurrence)) {
                ptr = std::max(static_cast<unsigned long long >(contig.coverage[contigID]), minOccurrence);
            }
        }
    }
    occurrenceDistribution.clear();
    occurrenceDistribution.resize(UINT16_MAX);
    kmerFP = platanus::makeTemporaryFile();
    this->writeKmerDistribution(kmerFP);

    for (unsigned short i = UINT16_MAX - 1; i > 0; --i) {
        if (occurrenceDistribution[i] > 0) {
            this->maxOccurrence = i;
            break;
        }
    }

    return upper;
}


//////////////////////////////////////////////////////////////////////////////////////
// load kmer from temporary file
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long Counter<KMER>::loadKmer(const unsigned long long minOccurrence, const unsigned long long doubleHashSize)
{

    unsigned short occurrence;
    unsigned long long total = 0;
    FILE *newKmerFP = platanus::makeTemporaryFile();

    std::cerr << "loading kmers..." << std::endl;
    KMER key(kmerLength);
    deleteOccurrenceTable();
    rewind(kmerFP);

    while (key.readTemporaryFileForward(kmerFP)) {
        fread(&occurrence, sizeof(unsigned short), 1, kmerFP);
        if (occurrence >= minOccurrence) {
            ++total;
        }
    }
    rewind(kmerFP);

    // define DoubelHash memory size (memory size must be 2^N)
    unsigned long long size = log(total / platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR) / log(2);
    size = std::pow(2, size + 1);
    if (size > doubleHashSize) {
        platanus::MemoryAlert();
    }

    occurrenceTable.resize(std::max(size, doubleHashSize));


    while (key.readTemporaryFileForward(kmerFP)) {
        fread(&occurrence, sizeof(unsigned short), 1, kmerFP);
        if (occurrence < minOccurrence) {
            continue;
        }
        occurrenceTable[key.forward] = occurrence;
    }
    fclose(kmerFP);
    kmerFP = newKmerFP;
    return size;
}





//////////////////////////////////////////////////////////////////////////////////////
// load contig
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Counter<KMER>::swapOccurrenceTable(const unsigned k, DoubleHash<typename KMER::keyType, unsigned short> &table)
{
    kmerLength = k;
    std::cerr << "K = " << kmerLength << ", loading kmers from contigs..." << std::endl;
     deleteOccurrenceTable();
    table.swap(occurrenceTable);
}


//////////////////////////////////////////////////////////////////////////////////////
// for multi thread function makeKmerReadMatchedEdgeKmer
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
unsigned long long Counter<KMER>::makeKmerReadDistributionConsideringPreviousGraph(const unsigned k, FILE **readFP, const unsigned long long memory, const unsigned long long numThread)
{
    platanus::SEQ seq;
    FILE *unmappedFP[numThread];
    FILE *tmpFP[numThread];
    bool loopFlag[numThread];
    KMER kmer(k);
    typename KMER::keyType key(k);
    std::cerr << "K = " << k << ", saving additional kmers(not found in contigs) from reads..." << std::endl;

    omp_lock_t lock[HASH_DIVID];
    if (kmerFP != NULL)
        fclose(kmerFP);
    kmerFP = platanus::makeTemporaryFile();

    for (unsigned long long i = 0; i < HASH_DIVID; ++i) {
        omp_init_lock(&lock[i]);
    }


    // check kmer used contig
    omp_set_num_threads(numThread);
    #    pragma omp parallel for private(seq) firstprivate(kmer, key) schedule(static, 1)
    for (unsigned long long i = 0; i < numThread; ++i) {
        unmappedFP[i] = platanus::makeTemporaryFile();
        divideKmerUsedMakingPreviousContig(k, kmer, key, loopFlag[i], readFP[i], unmappedFP[i]);
    }


    // write occurrence distribution
    occurrenceDistribution.clear();
    occurrenceDistribution.resize(UINT16_MAX);

    auto tableIterator = occurrenceTable.begin();
    auto tableEnd = occurrenceTable.end();
    for (; tableIterator != tableEnd; ++tableIterator) {
        if (tableIterator->second != 0) {
            kmer.writeKey(kmerFP, tableIterator->first);
            ++occurrenceDistribution[tableIterator->second];
            fwrite(&(tableIterator->second), sizeof(unsigned short), 1, kmerFP);
            tableIterator->second = 0;
        }
    }

    bool isFirst = true;
    unsigned iterateTimes = 32;
    while (1) {
        if (!isFirst) {
            iterateTimes += writeKmerDistribution(kmerFP);
        }
        isFirst = false;
        // check whether finish count kmer
        // if loopFlag on, part of kmer is written in temporary file because of used large memory
        unsigned long long i = 0;
        while (i < numThread) {
            if (loopFlag[i]) break;
            ++i;
        }
        if (i == numThread)
            break;

        #   pragma omp parallel for schedule(static, 1)
        for (unsigned long long threadID = 0; threadID < numThread; ++threadID) {
            loopFlag[threadID] = false;
            tmpFP[threadID] = platanus::makeTemporaryFile();
            countKmerPerThreadSecond(loopFlag[threadID], unmappedFP[threadID], tmpFP[threadID], lock, iterateTimes, numThread);
            fclose(unmappedFP[threadID]);
            unmappedFP[threadID] = tmpFP[threadID];
        }
    }

    // delete temporary file
    for (unsigned long long i = 0; i < numThread; ++i)
        fclose(unmappedFP[i]);

    // search max occurrence in kmer distribution
    for (unsigned i = UINT16_MAX - 1; i > 0; --i) {
        if (occurrenceDistribution[i] > 0) {
            this->maxOccurrence = i;
            break;
        }
    }

    for (unsigned long long i = 0; i < HASH_DIVID; ++i) {
        omp_destroy_lock(&lock[i]);
    }
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////
// count kmer multi thread
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void Counter<KMER>::countKmerPerThreadSecond(bool &loopFlag, FILE *readFP, FILE *unmappedFP, omp_lock_t lock[], const unsigned iterateTimes, const unsigned long long numThread)
{
    KMER kmer(kmerLength);
    rewind(readFP);
    while (kmer.readTemporaryFileForward(readFP)) {
        countKmerOrWriteTemporary(loopFlag, kmer.forward, unmappedFP, lock, kmer, iterateTimes);
    }
}




//////////////////////////////////////////////////////////////////////////////////////
// if this kmer(key) already exist, add
// else write temporary or add new key in table
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline void Counter<KMER>::countKmerOrWriteTemporary(bool &loopFlag, const typename KMER::keyType &key, FILE *unmappedFP, omp_lock_t lock[], const KMER &kmer, const unsigned iterateTimes)
{
    const unsigned mod = key & static_cast<unsigned long long>(HASH_DIVID - 1);
    omp_set_lock(&lock[mod]);
    auto it = occurrenceTable.find_times_any(key, iterateTimes);
    if (it->second == 0) {
        it->first = key;
        it->second = 1;
    } else if (it->first == key) {
        if (it->second != UINT16_MAX - 1) {
            ++(it->second);
        }
    } else {
            loopFlag = true;
            kmer.writeKey(unmappedFP, key);
    }
    omp_unset_lock(&lock[mod]);
}


//////////////////////////////////////////////////////////////////////////////////////
// write kmer distribution in temporary file
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
inline unsigned Counter<KMER>::writeKmerDistribution(FILE *kmerFP)
{
    KMER kmer(kmerLength);
    unsigned long long num = 0;
    // write kmer int temporary file
        auto tableIterator = occurrenceTable.begin();
        auto tableEnd = occurrenceTable.end();
        for (; tableIterator != tableEnd; ++tableIterator) {
            if (tableIterator->second != 0) {
                ++num;
                kmer.writeKey(kmerFP, tableIterator->first);
                fwrite(&(tableIterator->second), sizeof(unsigned short), 1, kmerFP);
                ++occurrenceDistribution[tableIterator->second];
                tableIterator->second = 0;
            }
        }

    // if DoubleHash is filled keys less than half (considering MAX_LOAD_FACTOR), add find iterate times
    if (num / HASH_DIVID < occurrenceTable.size() * platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR / 2)
        return 1;
    else 
        return 0;
}



//////////////////////////////////////////////////////////////////////////////////////
// divide kmer which contains contig
// if contain, use contig's coverage value, else, write temporary file
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Counter<KMER>::divideKmerUsedMakingPreviousContig(const unsigned kmerLength, KMER &kmer, typename KMER::keyType &key, bool &loopFlag, FILE *readFP, FILE *unmappedFP)
{
    platanus::SEQ seq;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    loopFlag = false;
    rewind(readFP);
    while(seq.readTemporaryFile(readFP)) {
        if (seq.length < kmerLength) continue;
        seq.positionUnknown.resize(seq.numUnknown + 1);
        seq.positionUnknown[seq.numUnknown] = platanus::ConstParam::MAX_READ_LEN + 1;
        seq.numUnknown = 0;
        for (unsigned j = 0; j < kmerLength - 1; ++j) {
            kmer.setForward(kmerLength - 2 - j, seq.base[j]);
            kmer.setReverse(j + 1, 0x3 ^ seq.base[j]);
        }
        for (unsigned j = 0; j < seq.length - kmerLength + 1; ++j) {
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.setForward(0, seq.base[j + kmerLength - 1]);
            kmer.reverse >>= 2;
            kmer.setReverse(kmerLength - 1, 0x3 ^ seq.base[j + kmerLength - 1]);
            if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) < j + kmerLength) {
                if (static_cast<unsigned>(seq.positionUnknown[seq.numUnknown]) <= j)
                    ++seq.numUnknown;
                continue;
            }
            key = std::min(kmer.forward, kmer.reverse);
            if (occurrenceTable.find_any(key)->second == 0) {
                loopFlag = true;
                kmer.writeKey(unmappedFP, key);
            }
        }
    }
}



//////////////////////////////////////////////////////////////////////////////////////
// pick up read which is matched edge kmer
// these read are only used in next step (only these reads make de brujin graph)
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Counter<KMER>::pickupReadMatchedEdgeKmer(FILE **readFP, const bool readFPCloseFlag)
{
    platanus::SEQ seq, tmpSeq;
    FILE *newReadFP = platanus::makeTemporaryFile();
    rewind(*readFP);
    KMER key(kmerLength);
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));


    while (seq.readTemporaryFile(*readFP)) {
        if (static_cast<unsigned>(seq.length) < kmerLength) continue;

        tmpSeq = seq;
        tmpSeq.positionUnknown.resize(tmpSeq.numUnknown + 1);
        tmpSeq.positionUnknown[tmpSeq.numUnknown] = platanus::ConstParam::MAX_READ_LEN + 1;
        tmpSeq.numUnknown = 0;
        for (unsigned i = 0; i < kmerLength - 1; ++i) {
            key.setForward(kmerLength - i - 2, tmpSeq.base[i]);
            key.setReverse(i + 1, 0x3 ^ tmpSeq.base[i]);
        }
        for (unsigned i = 0; i < tmpSeq.length - kmerLength + 1; ++i) {
            key.forward <<= 2;
            key.maskForward(mask);
            key.setForward(0, tmpSeq.base[i + kmerLength - 1]);
            key.reverse >>= 2;
            key.setReverse(kmerLength - 1, 0x3 ^ tmpSeq.base[i + kmerLength - 1]);

            if (static_cast<unsigned>(tmpSeq.positionUnknown[tmpSeq.numUnknown]) < i + kmerLength) {
                if (static_cast<unsigned>(tmpSeq.positionUnknown[tmpSeq.numUnknown]) <= i)
                    ++tmpSeq.numUnknown;
                continue;
            }
            if (findValue(std::min(key.forward, key.reverse)) > 0) {
                seq.writeTemporaryFile(newReadFP);
                break;
            }
        }
    }
	if (readFPCloseFlag)
		fclose(*readFP);
    *readFP = newReadFP;
}


//////////////////////////////////////////////////////////////////////////////////////
// sort kmer
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
FILE *Counter<KMER>::sortedKeyFromKmerFile(const unsigned long long minOccurrence)
{
    std::vector<typename KMER::keyType> kmerBuffer;
    unsigned short occurrence;
    unsigned long long totalKmer = 0;
    KMER kmer(kmerLength);
    FILE *sortedKeyFP;

    // reserve memory for kmer buffer
    for (unsigned i = minOccurrence; i <= maxOccurrence; ++i)
        totalKmer += occurrenceDistribution[i];
    deleteOccurrenceDistribution();
    kmerBuffer.resize(totalKmer, kmerLength);

    rewind(kmerFP);
    totalKmer = 0;

    // insert kmer from temporary file
    while(kmer.readTemporaryFileForward(kmerFP)) {
        fread(&occurrence, sizeof(unsigned short), 1, kmerFP);
        if (occurrence < minOccurrence) continue;
            kmerBuffer[totalKmer] = kmer.forward;
            ++totalKmer;
    }

    sort(kmerBuffer.begin(), kmerBuffer.end());

    // write temporary file
    sortedKeyFP = platanus::makeTemporaryFile();
    unsigned long long kmerBufferSize = kmerBuffer.size();
    for (unsigned long long i = 0; i < kmerBufferSize; ++i) {
        kmer.writeKey(sortedKeyFP, kmerBuffer[i]);
    }
    return sortedKeyFP;
}


//////////////////////////////////////////////////////////////////////////////////////
// output kmer occurrence distribution
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Counter<KMER>::outputOccurrenceDistribution(const std::string &filename)
{
    std::ofstream ofs(filename.c_str());
    for (unsigned long long i = 1; i <= maxOccurrence; ++i) {
        ofs << i << "\t" << occurrenceDistribution[i] << std::endl;
    }
    ofs.close();
}




//////////////////////////////////////////////////////////////////////////////////////
// count kmer for gap close
//////////////////////////////////////////////////////////////////////////////////////
template <typename KMER>
void Counter<KMER>::countKmerForGapClose(const std::vector<std::vector<char> > &seq)
{
	if (seq.empty())
		return;

    KMER kmer(kmerLength);
    unsigned start = 0;
    bool init;
    const unsigned long long mask = kmerLength >= 32 ? ~static_cast<unsigned long long>(0) : ~(~static_cast<unsigned long long>(0) << (2 * kmerLength));
    // +2 means to remove head and tail seq
    for (auto it = seq.begin() + GapClose::GAP_CANDIDATE_SEQ_1ST, end = seq.end(); it != end; ++it) {
        start = 0;
        init = true;
        if (it->size() < kmerLength) continue;
        while (start < it->size() - kmerLength + 1) {
            if (init) {
                unsigned i = 0;
                for (; i < kmerLength - 1; ++i) {
                    if ((*it)[start + i] == 4) break;
                    kmer.forward <<= 2;
                    kmer.maskForward(mask);
                    kmer.setForward(0, (*it)[start + i]);
                }

                if (i != kmerLength - 1) {
                    start += i + 1;
                    continue;
                }
                init = false;
            }
            if ((*it)[start + kmerLength - 1] == 4) {
                start += kmerLength;
                init = true;
                continue;
            }
            kmer.forward <<= 2;
            kmer.maskForward(mask);
            kmer.setForward(0, (*it)[start + kmerLength - 1]);
            if (occurrenceGapTable[kmer.forward] < UINT16_MAX - 1)
                ++occurrenceGapTable[kmer.forward];
            ++start;
        }
    }
}



#endif


