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

#ifndef KMER_H
#define KMER_H

#include "common.h"
#include "binstr.h"
#include <vector>
#include <iostream>
#include <utility>

typedef unsigned long long Kmer31Key;
typedef binstr_t KmerNKey;


template <typename K, typename H>
class KmerBase
{
public:
    typedef K keyType;
    typedef H hasher;


    K forward, reverse;
    unsigned kmerLength;
    KmerBase(): forward(0), reverse(0), kmerLength(0) {}
    KmerBase(long length): forward(0), reverse(0), kmerLength(length) {}
    KmerBase(long constArg, long length): forward(constArg), reverse(constArg), kmerLength(length) {}
    KmerBase(const KmerBase &) = default;
    KmerBase &operator=(const KmerBase &) = default;
    virtual ~KmerBase() = default;



    static void revCom_u64(unsigned long long &t) {
    t ^= UINT64_MAX;
    t = ((t&0x3333333333333333ull)<<2) | ((t&0xCCCCCCCCCCCCCCCCull)>>2);
    t = ((t&0x0F0F0F0F0F0F0F0Full)<<4) | ((t&0xF0F0F0F0F0F0F0F0ull)>>4);
    t = ((t&0x00FF00FF00FF00FFull)<<8) | ((t&0xFF00FF00FF00FF00ull)>>8);
    t = ((t&0x0000FFFF0000FFFFull)<<16) | ((t&0xFFFF0000FFFF0000ull)>>16);
    t = (t<<32) | (t>>32);
  }
    virtual void writeTemporaryFileForward(FILE *fp) const = 0;
    virtual void writeTemporaryFileReverse(FILE *fp) const = 0;

    virtual size_t readTemporaryFileForward(FILE *fp) = 0;
    virtual size_t readTemporaryFileReverse(FILE *fp) = 0;

    virtual size_t readKey(FILE *fp, K &key) const = 0;
    virtual size_t writeKey(FILE *fp, const K &key) const = 0;
    virtual void convertKeyToBinstr(binstr_t &b, const K &key) const = 0;

    virtual void setForward(const unsigned position, const unsigned char value) = 0;
    virtual void setReverse(const unsigned position, const unsigned char value) = 0;

    virtual void clearForward(void) = 0;
    virtual void maskForward(const unsigned long long mask) = 0;
    virtual void maskReverse(const unsigned long long mask) = 0;
};






class Kmer31 : public KmerBase<Kmer31Key, std::hash<Kmer31Key> >
{
public:
    typedef KmerBase<Kmer31Key, std::hash<Kmer31Key> > Base;
    Kmer31():Base() {}
    Kmer31(long length): KmerBase<Kmer31Key, std::hash<Kmer31Key> >(length)
    {
    }
    ~Kmer31() = default;

    void writeTemporaryFileForward(FILE *fp) const
    {
        fwrite(&this->forward, sizeof(Kmer31Key), 1, fp);
    }

    void writeTemporaryFileReverse(FILE *fp) const
    {
        fwrite(&this->reverse, sizeof(Kmer31Key), 1, fp);
    }

    size_t readTemporaryFileForward(FILE *fp)
    {
        return fread(&this->forward, sizeof(Kmer31Key), 1, fp);
    }

    size_t readTemporaryFileReverse(FILE *fp)
    {
        return fread(&this->reverse, sizeof(Kmer31Key), 1, fp);
    }

    size_t readKey(FILE *fp, Kmer31Key &key) const
    {
        return fread(&key, sizeof(Kmer31Key), 1, fp);
    }

    size_t writeKey(FILE *fp, const Kmer31Key &key) const
    {
        return fwrite(&key, sizeof(Kmer31Key), 1, fp);
    }

    void convertKeyToBinstr(binstr_t &b, const Kmer31Key &key) const
	{
		b.value[0] = key;
	}

    void setForward(const unsigned position, const unsigned char value)
    {
        forward = (forward & ~(0x3ull << ((position)*2))) | (static_cast<unsigned long long>(value) << ((position)*2));
    }

    void setReverse(const unsigned position, const unsigned char value)
    {
        reverse = (reverse & ~(0x3ull << ((position)*2))) | (static_cast<unsigned long long>(value) << ((position)*2));
    }

    void reverseComplement(void)
    {
        reverse = forward;
        Base::revCom_u64(reverse);
        if (kmerLength != 32)
            reverse = (reverse >> (64 - (kmerLength) * 2)) & ((0x1ull << (kmerLength * 2)) - 0x1ull);
    }

    bool operator<(const Kmer31 &a) const
    {
        if (this->forward == a.forward)
            return this->reverse < a.reverse;
        else
            return this->forward < a.forward;
    }
    static unsigned long long assign(const Kmer31Key &key)
    {
        return 0;
    }

    void clearForward(void)
    {
        forward = 0;
    }

    int setKmer(const platanus::SEQ &seq, const unsigned start, const long keyLength, const long seedLength)
    {
        return setKmer(seq.base, start, keyLength, seedLength);
    }

    int setKmer(const std::string &seq, const unsigned start, const long keyLength, const long seedLength)
    {
        forward = reverse = 0;
        for (long i = 0; i < keyLength; ++i) {
            if (seq[start + i] == 4)
                return 1;
            forward = forward << 2 | (unsigned long long)seq[start + i];
            reverse = reverse << 2 | (unsigned long long)(0x3 ^ seq[start + keyLength - i - 1]);
        }
        return 0;
    }
    void maskForward(const unsigned long long mask)
    {
        forward &= mask;
    }

    void maskReverse(const unsigned long long mask)
    {
        reverse &= mask;
    }
    void show(Kmer31Key &a)
    {}
};

template <typename BASE>
class KmerN : public KmerBase<BASE, typename BASE::hasher>
{
public:
    typedef KmerBase<BASE, typename BASE::hasher> Base;
    KmerN(long length): KmerBase<BASE, typename BASE::hasher>(length, length)
    {
    }
    ~KmerN() = default;

    void writeTemporaryFileForward(FILE *fp) const
    {
        const unsigned long long *v = this->forward.value;
        fwrite(v, sizeof(Kmer31Key), (this->forward.len + 31) / 32, fp);
    }

    void writeTemporaryFileReverse(FILE *fp) const
    {
        const unsigned long long *v = this->reverse.value;
        fwrite(v, sizeof(Kmer31Key), (this->reverse.len + 31) / 32, fp);
    }

    size_t readTemporaryFileForward(FILE *fp)
    {
        unsigned long long *v = this->forward.value;
        return fread(v, sizeof(Kmer31Key), (this->forward.len + 31) / 32, fp);
    }

    size_t readTemporaryFileReverse(FILE *fp)
    {
        unsigned long long *v = this->reverse.value;
        return fread(v, sizeof(Kmer31Key), (this->reverse.len + 31) / 32, fp);
    }

    size_t readKey(FILE *fp, BASE &key) const
    {
        Kmer31Key *v = key.value;

        for (unsigned long i = 0; i < (key.len + 31) / 32; ++i) {
            fread(v++, sizeof(Kmer31Key), 1, fp);
        }
        return (key.len + 31) / 32;
    }

    size_t writeKey(FILE *fp, const BASE &key) const
    {
        const Kmer31Key *v = key.value;

        for (unsigned long i = 0; i < (key.len + 31) / 32; ++i) {
            fwrite(v++, sizeof(Kmer31Key), 1, fp);
        }
        return (key.len + 31) / 32;
    }

    void convertKeyToBinstr(binstr_t &b, const BASE &key) const
	{
		unsigned long long n = (b.len + 31) / 32;
		for (unsigned long long i = 0; i < n; ++i)
		  b.value[i] = key.value[i];
	}

    void setForward(const unsigned position, const unsigned char value)
    {
        this->forward.set(position, value);
    }

    void setReverse(const unsigned position, const unsigned char value)
    {
        this->reverse.set(position, value);
    }

    void reverseComplement(void)
    {
        this->reverse = this->forward;
        int i, j, n = (this->reverse.len + 31) / 32;
        if (n > 0) {
            Kmer31Key *v = this->reverse.value;
            for (i=0; i < n / 2; i++) {
                std::swap(*(v + i), *(v + n - 1 - i));
            }
            for (i=0; i<n; i++) {
                Base::revCom_u64(*v++);
            }
            j = (this->reverse.len%32) * 2;
            if (j>0) {
                v = this->reverse.value;
                for (i=0; i<n-1; i++) {
                    *v = ((*v)>>(64-j)) | ((*(v+1))<<j);
                    v++;
                }
                (*v)>>=(64-j);
            }
        }
    }

    bool operator<(const KmerN &a) const
    {
        return this->forward < a.forward;
    }

    static unsigned long long assign(const BASE &key)
    {
        return (key.value[1] & 1023);
    }

    void clearForward(void)
    {
        this->forward.clear();
    }

    void maskForward(const unsigned long long mask)
    {
    }

    void maskReverse(const unsigned long long mask)
    {
    }
};



#endif
