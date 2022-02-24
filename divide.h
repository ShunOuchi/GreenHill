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

#ifndef DIVIDE_H
#define DIVIDE_H

#include "baseCommand.h"
#include <memory>
#include <vector>
#include <string>
#include <sstream>

class Divide : public BaseCommand
{
private:
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
    Divide();
    Divide(const Divide &) = delete;
    Divide &operator=(const Divide &) = delete;
    ~Divide() = default;

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

	void execSolveDBG(void);
};



#endif
