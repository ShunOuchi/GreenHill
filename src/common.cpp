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

#include "common.h"

//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const std::string platanus::ConstParam::VERSION = "1.1.0";

const unsigned long long platanus::ConstParam::MAX_READ_LEN = 500000;
const unsigned platanus::ConstParam::SCAFFOLD_HASH_OVERLAP = 32;
const unsigned platanus::ConstParam::OUTPUT_LINE_LENGTH = 80;
const unsigned platanus::ConstParam::MAX_FILE_NUM = 100;
const unsigned platanus::ConstParam::MAX_FILE_LEN = 200;
const unsigned platanus::ConstParam::MAX_THREAD = 100;
const double platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR = 0.9;
const unsigned long long platanus::ConstParam::DEFAULT_CONTIG_READ_LEN = 100;
const double platanus::ConstParam::SHORT_READ_INS_SIZE_LOWER_BOUND_FACTOR = 0.25;
const double platanus::ConstParam::SHORT_READ_INS_SIZE_UPPER_BOUND_FACTOR = 1.75;
const double platanus::ConstParam::LONG_READ_INS_SIZE_LOWER_BOUND_FACTOR = 0.25;
const double platanus::ConstParam::LONG_READ_INS_SIZE_UPPER_BOUND_FACTOR = 10.0;
std::string platanus::globalTmpFileDir = ".";


void platanus::setLimitNumFileOpenMax()
{
	struct rlimit rl;

	getrlimit(RLIMIT_NOFILE, &rl);

	rl.rlim_cur = rl.rlim_max;
	setrlimit(RLIMIT_NOFILE, &rl);
}


void platanus::setFastaFileNameAndNumber(const vector<string> fastaFileName, std::vector<std::string> &scaffoldFileName, std::vector<long> &numScaffoldInFile)
{
	scaffoldFileName.assign(fastaFileName.size(), "");
	numScaffoldInFile.assign(fastaFileName.size(), 0);

	for (unsigned fileIndex = 0; fileIndex < fastaFileName.size(); ++fileIndex) {
		scaffoldFileName[fileIndex] = fastaFileName[fileIndex];
		FILE* fp = platanus::openFileAllowingCompression(fastaFileName[fileIndex], "r");

		string line;
		while (platanus::getlineFILE(line, fp) != NULL) {
			if (line[0] == '>')
				++numScaffoldInFile[fileIndex];
		}

		platanus::closeFileAllowingCompression(fp, fastaFileName[fileIndex]);
	}
}


std::string platanus::shellQuote(const std::string &argument)
{
	std::string quoted("'");
	for (std::string::const_iterator itr = argument.begin(); itr != argument.end(); ++itr) {
		if (*itr == '\'')
			quoted += "'\\''";
		else
			quoted += *itr;
	}
	quoted += "'";
	return quoted;
}


platanus::FILECOMPRESSION platanus::checkFileCompression(const std::string &filename)
{
	FILE *fp = fopen(filename.c_str(), "rb");
	if (fp == NULL)
		throw platanus::FILEError(filename);

	unsigned char magic[6] = {0};
	size_t length = fread(magic, 1, sizeof(magic), fp);
	fclose(fp);

	// gzip's magic bytes also identify BGZF, which is a gzip-compatible format.
	if (length >= 2 && magic[0] == 0x1f && magic[1] == 0x8b)
		return platanus::FILECOMPRESSION::GZIP;
	if (length >= 3 && magic[0] == 0x42 && magic[1] == 0x5a && magic[2] == 0x68)
		return platanus::FILECOMPRESSION::BZIP2;
	if (length >= 6 && magic[0] == 0xfd && magic[1] == 0x37 && magic[2] == 0x7a
		&& magic[3] == 0x58 && magic[4] == 0x5a && magic[5] == 0x00)
		return platanus::FILECOMPRESSION::XZ;
	if (length >= 4 && magic[0] == 0x28 && magic[1] == 0xb5 && magic[2] == 0x2f && magic[3] == 0xfd)
		return platanus::FILECOMPRESSION::ZSTD;
	if (length >= 4 && ((magic[0] == 0x04 && magic[1] == 0x22 && magic[2] == 0x4d && magic[3] == 0x18)
		|| (magic[0] == 0x02 && magic[1] == 0x21 && magic[2] == 0x4c && magic[3] == 0x18)))
		return platanus::FILECOMPRESSION::LZ4;

	return platanus::FILECOMPRESSION::UNCOMPRESSED;
}


std::string platanus::getDecompressionCommand(const platanus::FILECOMPRESSION format, const std::string &filename)
{
	const char *program = NULL;
	switch (format) {
	case platanus::GZIP: program = "gzip"; break;
	case platanus::BZIP2: program = "bzip2"; break;
	case platanus::XZ: program = "xz"; break;
	case platanus::ZSTD: program = "zstd"; break;
	case platanus::LZ4: program = "lz4"; break;
	default: return std::string();
	}

	// Prefix relative paths so filenames beginning with '-' are not parsed as options.
	std::string input = filename;
	if (!input.empty() && input[0] != '/')
		input = "./" + input;

	return std::string(program) + " -cd " + platanus::shellQuote(input);
}


FILE* platanus::openFileAllowingCompression(const std::string &filename, const char *mode)
{
	platanus::FILECOMPRESSION format = checkFileCompression(filename);

	FILE *fp;

	if (format == platanus::FILECOMPRESSION::UNCOMPRESSED) {
		fp = fopen(filename.c_str(), mode);
	}
	else {
		std::string command = getDecompressionCommand(format, filename);
		fp = popen(command.c_str(), mode);
	}

	if (fp == NULL)
		throw platanus::FILEError(filename);

	return fp;
}


int platanus::closeFileAllowingCompression(FILE *fp, const std::string &filename)
{
	platanus::FILECOMPRESSION format = checkFileCompression(filename);

	if (format != platanus::FILECOMPRESSION::UNCOMPRESSED) {
		return pclose(fp);
	}

	return fclose(fp);
}


FILE *platanus::getlineFILE(std::string &str, FILE *fp)
{
	char *buffer = NULL;
	size_t allocSize = 0;

	ssize_t lineLength = getline(&buffer, &allocSize, fp);
	if (lineLength <= 0)
		return NULL;

	if (buffer[lineLength - 1] == '\n')
		buffer[lineLength - 1] = '\0';

	// CRLFの場合も念の為、対処
	if (lineLength > 1) {
		if (buffer[lineLength - 2] == '\r')
			buffer[lineLength - 2] = '\0';
	}

	str = buffer;
	free(buffer);

	return fp;
}


void platanus::printNumOpenFiles()
{
	std::cerr << "DEBUG: Number of files:" << std::endl;

    std::ostringstream oss;
	oss << "ls /proc/" << getpid() << "/fd | wc -l >&2";
	system(oss.str().c_str());

	std::cerr << std::endl;
}
