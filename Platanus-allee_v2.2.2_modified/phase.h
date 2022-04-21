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

#ifndef PHASE_H
#define PHASE_H

#include "baseCommand.h"
#include <memory>
#include <vector>
#include <string>
#include <sstream>

class Phase : public BaseCommand
{
private:
    std::string directoryName;
    std::string previousDirectoryName;
    std::string platanusBinary;


    virtual bool checkFileEnough(void)
    {
        if (optionMultiArgs["-c"].empty()) {
            std::cerr << "Error: No contig fasta files (-c) were specified!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    Phase();
    Phase(const Phase &) = delete;
    Phase &operator=(const Phase &) = delete;
    ~Phase() = default;

    void usage(void) const;
    void exec(void);

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            this->platanusBinary = argv[0];
            return true;
        }
        return false;
    }


	void setDirectoryName(const std::string newDirectoryName);
    void createDirectory(void) const;
	void execSolveDBG(const unsigned long times, const unsigned long numIterate);
    void execGapClose(const unsigned long times, const unsigned long numIterate);
	void moveAndConcatenateFinalRoundResult(const std::string intermediateDirectoryName, const unsigned long times);
};



#endif
