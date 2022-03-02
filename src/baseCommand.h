/*
Copyright (C) 2022 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GreenHill.

GreenHill is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GreenHill is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GreenHill; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef BASECOMMAND_H
#define BASECOMMAND_H

#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <omp.h>
#include <cstdint>
#include "common.h"




//////////////////////////////////////////////////////////////////////////////////////
// Superclass for each command (assemble, scaffold, gap_close)
//////////////////////////////////////////////////////////////////////////////////////
class BaseCommand
{
    typedef unsigned long long u64_t;
private:



//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// PairFile Class
// reserve paired-file name and library type(inner pair or outer one)
    struct PairFile
    {
        int libraryID;
        std::string libraryType;
        std::string fileFirst;
        std::string fileSecond;

        PairFile(const int id, const std::string &type, const std::string &name1, const std::string &name2=""): libraryID(id), libraryType(type), fileFirst(name1), fileSecond(name2) {}
        PairFile(const PairFile &pair) = default;
        PairFile &operator=(const PairFile &pair) = default;
        ~PairFile() {}

        bool operator<(const PairFile &a) const
        {
            return (libraryID - a.libraryID) < 0;
        }
    };
//////////////////////////////////////////////////////////////////////////////////////


    virtual bool checkFileEnough(void) = 0;
    virtual int checkOtherOption(char *argv) const = 0;




public:
    // constructor and destructor
    BaseCommand() {}
    virtual ~BaseCommand() {}
    BaseCommand(const BaseCommand &base) = delete;
    BaseCommand *operator=(const BaseCommand &base) = delete;


    // pure virtual functions
    virtual void exec() = 0;
    virtual void usage() const = 0;



    platanus::FILETYPE checkFileFormat(const std::string &filename) const;
    virtual bool parseArgs(int argc, char **argv);

protected:
    std::unordered_map<std::string, bool> optionBool;
    std::unordered_map<std::string, std::string> optionSingleArgs;
    std::unordered_map<std::string, std::vector<std::string> > optionMultiArgs;
    std::vector<PairFile> optionPairFile;
    std::vector<std::string> pairedEndSingleFileType;
    std::vector<std::string> pairedEndPairFileType;



    int divideArgvInt(char *args) const;

};




#endif
