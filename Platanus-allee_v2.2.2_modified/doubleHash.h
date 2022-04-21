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

#ifndef DOUBLEHASH_H
#define DOUBLEHASH_H

#include "common.h"
#include <memory>

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// DoubleHash class
// this class is hash_like table
// for accelerate, all member variable are assigned public
// VAL type need bool cast operator bool() and swap function
// hash table size must be 2 ^ N
template <typename KEY, typename VAL>
struct DoubleHash
{
    typedef unsigned long long u64_t;
private:
    u64_t indexSize;
    u64_t indexLength;
    u64_t shifter;
    std::unique_ptr<std::pair<KEY, VAL>[]> table;


    //////////////////////////////////////////////////////////////////////////////////////
    // check whether table size is 2 ^ N
    //////////////////////////////////////////////////////////////////////////////////////
    void checkSizeIspow2(const size_t n)
    {
        if (n < 0) throw platanus::DoubleHashError();
        if (n > 0) {
            if (((n) & (n -1)) != 0) {
                throw platanus::DoubleHashError();
            }
        }
    }

public:


    DoubleHash(): indexSize(0), indexLength(0), shifter(0), table() {}
    DoubleHash(const u64_t len): indexSize(len - 1), indexLength(calcLength(len)), shifter(0), table(new std::pair<KEY, VAL>[len]()) {
        shifter = indexLength >= 32 ? 0 : 2 * indexLength;
        checkSizeIspow2(len);
    }
    DoubleHash(const DoubleHash &a): indexSize(a.indexSize), indexLength(a.indexLength), shifter(a.shifter), table() {}
    ~DoubleHash() = default;


    // define inde key seed
    u64_t calcLength(const u64_t len)
    {
        for (u64_t i = 1; i < 64; ++i) {
            if (len >> i == 0) {
                return i;
            }
        }
        return 64;
    }

    // make hash key
    template <typename ANY>
    u64_t makeHashKey(const ANY &key) const
    {
        u64_t value = 0;
        for (u64_t i = 0, n = (key.len + 31) / 32; i < n; ++i)
            value += (key.value[i] + (key.value[i] >> indexLength) + (key.value[i] >> shifter));
        return value & indexSize;
    }

    u64_t makeHashKey(const u64_t &key) const
    {
        return (key + (key >> indexLength) + (key >> shifter)) & indexSize;
    }


    // remake hash key when hash value are crashed
    template <typename ANY>
    u64_t reHashKey(const ANY &key) const
    {
        u64_t value = 0;
        for (u64_t i = 0, n = (key.len + 31) / 32; i < n; ++i)
            value += (~key.value[i] ^ (key.value[i] >> indexLength) ^ (key.value[i] >> shifter));
        return value | 1;
    }

    u64_t reHashKey(const u64_t &key) const
    {
        return (~key ^ (key >> indexLength) ^ (key >> shifter)) | 1;
    }

    std::pair<KEY, VAL> *begin(void)
    {
        return table.get();
    }

    const std::pair<KEY, VAL> *begin(void) const
    {
        return table.get();
    }

    std::pair<KEY, VAL> *end(void)
    {
        return table.get() + indexSize + 1;
    }

    const std::pair<KEY, VAL> *end(void) const
    {
        return table.get() + indexSize + 1;
    }

    // find function
    // when key is not found, return unoccupied region anywhere
    std::pair<KEY, VAL> *find_any(const KEY &key)
    {
        auto it = begin();
        u64_t value = makeHashKey(key);
        if (!table[value].second || key == table[value].first)
            return it + value;
        u64_t step = reHashKey(key);
        value = (value + step) & indexSize;
        while (table[value].second) {
            if (key == table[value].first)
                return it + value;
            value = (value + step) & indexSize;
        }
        return it + value;
    }

    const std::pair<KEY, VAL> *find_any(const KEY &key) const
    {
        auto it = begin();
        u64_t value = makeHashKey(key);
        if (!table[value].second || key == table[value].first)
            return it + value;
        u64_t step = reHashKey(key);
        value = (value + step) & indexSize;
        while (table[value].second) {
            if (key == table[value].first)
                return it + value;
            value = (value + step) & indexSize;
        }
        return it + value;
    }

    std::pair<KEY, VAL> *find_times_any(const KEY &key, const unsigned times)
    {
        auto it = begin();
        unsigned i = 1;
        u64_t value = makeHashKey(key);
        if (!table[value].second || key == table[value].first)
            return it + value;
        u64_t step = reHashKey(key);
        value = (value + step) & indexSize;
        while (table[value].second && i < times) {
            ++i;
            if (key == table[value].first)
                return it + value;
            value = (value + step) & indexSize;
        }
        return it + value;
    }


    VAL &operator[](const KEY &key)
    {
        auto it = find_any(key);
        it->first = key;
        return it->second;
    }


    void resize(const size_t n)
    {
        checkSizeIspow2(n);
        table.reset(new std::pair<KEY, VAL>[n]());
        indexSize = n - 1;
        indexLength = calcLength(n);
        shifter = indexLength >= 32 ? 0 : 2 * indexLength;
    }

    void clear(void)
    {
        this->clean();
        indexSize = 0;
        indexLength = 0;
        shifter = 0;
    }

    void clean(void)
    {
        for (auto it = table.begin(), end = table.end(); it != end; ++it) {
            it->first = it->second = 0;
        }
    }

    size_t size(void) const
    {
        return indexSize + 1;
    }

    void swap(DoubleHash<KEY, VAL> &sw)
    {
        sw.table.swap(table);
        std::swap(indexSize, sw.indexSize);
        std::swap(indexLength, sw.indexLength);
        std::swap(shifter, sw.shifter);
    }

};




#endif


