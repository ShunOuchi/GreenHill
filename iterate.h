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

#ifndef ITERATE_H
#define ITERATE_H

#include "baseCommand.h"
#include "scaffold.h"
#include "gapClose.h"
#include "merge.h"
#include "polish.h"
#include <memory>
#include <sstream>

class IterateScaffold : public BaseCommand
{
private:
    static const std::string CONTIG_FOOTER;
    static const std::string SCAF_FOOTER;
    static const std::string POLISH_FOOTER;
    static const std::string GAP_FOOTER;
    static const std::string EX_FOOTER;
    static const std::string DIV_FOOTER;
    static const std::string MERGE_FOOTER;
    static const std::string ITERATION_FOOTER;
    std::string directoryName;
    std::string previousDirectoryName;
    std::string platanusBinary;


    virtual bool checkFileEnough(void)
    {
        if (optionSingleArgs["-c"] == "") {
            std::cerr << "Error: not specified contig file!!" << std::endl;
            return false;
        }
        return true;
    }

    virtual int checkOtherOption(char *argv) const
    {
        return 0;
    }


public:
    IterateScaffold();
    IterateScaffold(const IterateScaffold &) = delete;
    IterateScaffold &operator=(const IterateScaffold &) = delete;
    ~IterateScaffold() = default;

    void usage(void) const;
    void exec(void);

    void setDirectoryName(const unsigned long long times)
    {
        std::ostringstream oss;
        oss << optionSingleArgs["-o"] << times;
        previousDirectoryName = directoryName;
        directoryName = oss.str();
    }

    virtual bool parseArgs(int argc, char **argv)
    {
        if (BaseCommand::parseArgs(argc, argv)) {
            this->platanusBinary = argv[0];
            return true;
        }
        return false;
    }


    void createDirectory(void) const;
    void createContig(const unsigned long long times);
    void execDivideMode(void);
    void execMergeMode(const unsigned long long times);
    void execScaffoldMode(void);
    void execPolishMode(void);
	void execFinalPolishMode(void);
    void execGapCloseMode(const std::string &mode);
	void execClusterMode(void);
	void execPostClusteringScaffoldMode(void);

};



#endif
