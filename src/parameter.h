/**********************************************************************
 * Skewer - a fast and accurate adapter trimming tool
 *          using the bit-masked k-difference matching algorithm
 * Copyright (c) 2013-2014 by Hongshan Jiang
 * hongshan.jiang@gmail.com
 *
 * If you use this program, please cite the paper:
 * Jiang, H., Lei, R., Ding, S.W. and Zhu, S. (2014) Skewer: a fast and
 * accurate adapter trimmer for next-generation sequencing paired-end reads.
 * BMC Bioinformatics, 15, 182.
 * http://www.biomedcentral.com/1471-2105/15/182
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#ifndef _PARAMETER_H
#define _PARAMETER_H

#include <vector>
#include <string>

#include "common.h"
#include "fastq.h"

using namespace std;

typedef struct{
	const char * name;
	char chr;
}OPTION_ITEM;

typedef enum{
	COMPRESS_NONE = 0,
	COMPRESS_GZ = 1,
	COMPRESS_BZ2 = 2
}COMPRESS_FORMAT;

///////////////////////////////////////
class cParameter
{
private:
	char arr[2][MAX_PATH+1];
	int argc;
	char ** argv;
	bool bXFile;
	bool bYFile;
	bool bJFile;
	bool bAutoFormat;

public:
	const char * version;
	string x_str;
	string y_str;
	string m_str;
	string j_str;
	vector<string> adapters;
	vector<string> adapters2;
	vector< vector<bool> > bMatrix;
	vector<string> rowNames;
	vector<string> colNames;
	vector<string> juncAdapters;
	char * input[2];
	// These are the output fastq files for the
	// trimmed reads that pass filters
	vector<string> output;
	vector<string> output2;
	// These are the output fastq files for only
	// the reads that passed filters AND were trimmed.
	// The full reads appear in the file with the
	// trimmed bases in lower case rather than removed.
	vector<string> masked;
	vector<string> masked2;
	bool bWriteMasked;
	// These files contain the reads that were excluded
	// from the output because they failed filters.
	vector<string> excluded;
	vector<string> excluded2;
	bool bWriteExcluded;
	vector<string> barcodeNames;
	string barcodes;
	string mapfile;
	string untrimmed;
	string untrimmed2;
	char logfile[MAX_PATH+1+100];
	const char * pDecorate;
	TRIM_MODE trimMode;
	bool bShareAdapter;
	bool bBarcode;
	bool bClip;
	bool bQiime;
	bool bFilterNs;
	bool bFilterUndetermined;
	bool bRedistribute;
	bool bStdin;
	bool bStdout;
	bool bQuiet;
	bool bEnquireVersion;
	int nFileCnt;
	enum FASTQ_FORMAT fastqFormat;
	int baseQual;
	COMPRESS_FORMAT outputFormat;
	double epsilon, delta;
	int minLen, maxLen, minAverageQual, minEndQual, nThreads;
	int minK;
	int iCutF, iCutR;
	bool bCutTail;
	bool bFillWithNs;

private:
	char * occOfLastDot(char * str);
	bool IsDirectorySpecified (char * str);
	int ReadMatrix(const char * fileName);
	int ReadFasta(const char * fileName, vector<string> & sequences);

public:
	cParameter();
	bool IsAutoFastqFormat();
	void PrintVersion(FILE * fp);
	void PrintUsage(char * program, FILE *fp);
	void PrintSimpleUsage(char * program, FILE *fp);
	void printCommandLine(FILE *fp);
	void printRelatedFiles(FILE *fp);
	void printVersion(FILE *fp);
	void printLogo(FILE *fp, bool bLeadingRtn=false);
	void printOpt(FILE *fp, bool bLeadingRtn=false);
	int GetOpt(int argc, char *argv[], char * errMsg);
};

extern "C" int color_fprintf(int colorCode, FILE *stream, const char *format, ...);
extern "C" int color_sprintf(int colorCode, char *str, const char *format, ...);

#endif // _PARAMETER_H
