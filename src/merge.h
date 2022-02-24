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

#ifndef MERGE_H
#define MERGE_H

#include "common.h"
#include "counter.h"
#include "graph.h"


class ContigMerger : public BaseCommand
{
private:
    platanus::Contig contig;
    unsigned long long readLength;
    unsigned long long contigmaxK;

public:
    ContigMerger();
    ContigMerger(const ContigMerger &) = delete;
    ~ContigMerger() = default;

    void usage(void) const;
    void exec(void);
    template <typename KMER> void exec2_ForIntegrateKmer(Counter<KMER> &counter, BruijnGraph<KMER> &graph, const unsigned long long kmerLength);


    bool checkFileEnough(void)
    {
        if (optionMultiArgs["-f"].size() == 0) {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        return true;
    }

    int checkOtherOption(char *argv) const
    {
        return 0;
    }

};





#endif

