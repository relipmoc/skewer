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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <string>
#include <stdarg.h>

#include "parameter.h"
#include "fastq.h"

using namespace std;

const char * VERSION = "0.1.114";
const char * DATE = "March 18, 2014";
const char * AUTHOR = "Hongshan Jiang";

const char * ILLUMINA_ADAPTER_PREFIX = "AGATCGGAAGAGC";
const char * ILLUMINA_PAIR1_ADAPTER_PREFIX = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
const char * ILLUMINA_PAIR1_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
const char * ILLUMINA_PAIR2_ADAPTER_PREFIX = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA";
const char * ILLUMINA_PAIR2_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";
const char * ILLUMINA_SRNA_ADAPTER = "TCGTATGCCGTCTTCTGCTTGT";
const char * ILLUMINA_JUNCTION_ADAPTER = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";

// colorCode:
// -1: normal
// 0: black
// 1: red
// 2: green
// 3: yellow
// 4: blue
// 5: magenta
// 6: cyan
// 256: bold
int color_fprintf(int colorCode, FILE *stream, const char *format, ...)
{
	va_list arglist;
	bool bBold = (colorCode & 256) != 0;
	colorCode &= ~256;
	if(colorCode > 7){
		colorCode = 7;
	}
	va_start(arglist, format);
	if(colorCode >= 0)
		fprintf(stream, "\033[%d;3%dm", bBold, colorCode);
	int iRet = vfprintf(stream, format, arglist);
	if(colorCode >= 0)
		fprintf(stream, "\033[0m");
	va_end(arglist);
	return iRet;
}

int color_sprintf(int colorCode, char *str, const char *format, ...)
{
	va_list arglist;
	bool bBold = (colorCode & 256) != 0;
	colorCode &= ~256;
	if(colorCode > 7){
		colorCode = 7;
	}
	va_start(arglist, format);
	if(colorCode >= 0)
		sprintf(str, "\33[%d;3%dm", bBold, colorCode);
	int iRet = vsprintf(str, format, arglist);
	if(colorCode >= 0)
		sprintf(str, "\33[0m");
	va_end(arglist);
	return iRet;
}

void color_fprintf_sequences(int colorCode, FILE *stream, vector<string> &sequences)
{
	vector<string>::iterator it_seq;
	for(it_seq=sequences.begin(); it_seq!=sequences.end(); it_seq++){
		color_fprintf(colorCode, stream, "\t%s\n", (*it_seq).c_str());
	}
}

///////////////////////////////////////
cParameter::cParameter()
{
	version = VERSION;
	argc = 0;
	argv = NULL;

	input[0] = arr[0]; input[1] = arr[1];
	input[0][0] = input[1][0] = '\0';
	logfile[0] = '\0';

	trimMode = TRIM_DEFAULT;
	bShareAdapter = false;
	bBarcode = false;
	bFilterNs = false;
	bFilterUndetermined = false;
	bXFile = bYFile = bJFile = false;

	nFileCnt = 0;
	fastqFormat = UNKNOWN_FASTQ;
	outputFormat = COMPRESS_NONE;
	bStdin = false;
	bStdout = false;
	bQuiet = false;
	bEnquireVersion = false;
	bAutoFormat = true;
	baseQual = 33;
	epsilon = 0.1;
	delta = 0.03;
	minLen = 18;
	maxLen = 0;
	minAverageQual = 0;
	minEndQual = 0;
	minK = 5;
	nThreads = 1;
}

char * cParameter::occOfLastDot	(char * str)
{
	char * ret = NULL;
	do{
		if(*str == '.'){
			ret = str;
		}
		else if(*str == '/'){
			ret = NULL;
		}
	} while (*str++);
	return (ret != NULL) ? ret : str;
}

bool cParameter::ReadFasta(const char * fileName, vector<string> & sequences)
{
	char * line = NULL;
	size_t alloc;
	FILE * fp = fopen(fileName, "r");
	if(fp == NULL)
		return false;
	int len;
	int no = -1;
	string seq;
	while( (len = getline(&line, &alloc, fp)) > 0 ){
		if(line[0] == '>'){
			seq.clear();
			no = 0;
			continue;
		}
		line[--len] = '\0';
		if( (len > 0) && (line[len-1] == '\r') ){
			line[--len] = '\0';
		}
		if(no == 0){
			if(int(sequences.size()) == MAX_ADAPTER_CNT){
				fprintf(stderr, "\rWarning: only uses the first %d adapter sequences in \"%s\"\n", MAX_ADAPTER_CNT, fileName);
				break;
			}
			if(len > MAX_ADAPTER_LEN){
				if( (trimMode & TRIM_ANY) == TRIM_HEAD )
					seq.assign(line + len - MAX_ADAPTER_LEN);
				else
					seq.assign(line, 0, MAX_ADAPTER_LEN);
			}
			else{
				seq.assign(line);
			}
			sequences.push_back(seq);
		}
		no++;
	}
	if(line != NULL)
		free(line);
	fclose(fp);
	return true;
}

///////////////////////////////////////////
// public subroutines

bool cParameter::IsAutoFastqFormat()
{
	return bAutoFormat;
}

void cParameter::PrintVersion(FILE * fp)
{
	fprintf(fp, "\tskewer version: %s\n", VERSION);
	fprintf(fp, "\tAuthor: %s\n", AUTHOR);
	fprintf(fp, "\tLast update: %s\n", DATE);
}

void cParameter::PrintUsage(char * program, FILE * fp)
{
	fprintf(fp, "Skewer (A fast and sensitive adapter trimmer for paired-end reads)\n");
	fprintf(fp, "Version %s (updated in %s), Author: %s\n\n", VERSION, DATE, AUTHOR);
	fprintf(fp, "USAGE: %s [options] <reads.fastq> [paired-reads.fastq]\n", program);
	fprintf(fp, "    or %s [options] - (for input from STDIN)\n\n", program);
	fprintf(fp, "OPTIONS (ranges in brackets, defaults in parentheses):\n");
	fprintf(fp, " Adapter:\n");
	fprintf(fp, "          -x <str> Adapter sequence/file (%s)\n", ILLUMINA_PAIR1_ADAPTER_PREFIX);
	fprintf(fp, "          -y <str> Adapter sequence/file for pair-end reads (%s),\n", ILLUMINA_PAIR2_ADAPTER_PREFIX);
	fprintf(fp, "                   implied by -x if -x is the only one specified explicitly.\n");
	fprintf(fp, "          -j <str> Junction adapter sequence/file for Nextera Mater Pair reads (%s)\n", ILLUMINA_JUNCTION_ADAPTER);
	fprintf(fp, "          -m, --mode <str> trimming mode; 1) single-end -- head: 5' end; tail: 3' end; any: anywhere (tail)\n");
	fprintf(fp, "                           2) paired-end -- pe: paired-end; mp: mate-pair (pe)\n");
	fprintf(fp, " Tolerance:\n");
	fprintf(fp, "          -r <num> Maximum allowed error rate (normalized #errors / length of aligned region) [0, 0.5], (0.1)\n");
	fprintf(fp, "          -d <num> Maximum allowed indel error rate [0, r], (0.03)\n");
	fprintf(fp, "                   reciprocal is used for -r and -d when num > or = 2\n");
	fprintf(fp, "          -k <int> Minimum overlap length for adapter detection [1, inf);\n");
	fprintf(fp, "                   (max(1, int(4-10*r)) for single-end; (<junction length>/2) for mate-pair)\n");
	fprintf(fp, " Filtering & Post-trimming:\n");
	fprintf(fp, "          -q, --end-quality  <int> Trim 3' end until specified or higher quality reached; (0)\n");
	fprintf(fp, "          -Q, --mean-quality <int> The lowest mean quality value allowed before trimming; (0)\n");
	fprintf(fp, "          -l, --min <int> The minimum read length allowed after trimming; (18)\n");
	fprintf(fp, "          -L, --max <int> The maximum read length allowed after trimming; (no limit)\n");
	fprintf(fp, "          -n  Whether to filter out highly degerative (many Ns) reads; (no)\n");
	fprintf(fp, "          -u  Whether to filter out undetermined mate-pair reads; (no)\n");
	fprintf(fp, " Input/Output:\n");
	fprintf(fp, "          -f, --format <str> Format of FASTQ quality value: sanger|solexa|auto; (auto)\n");
	fprintf(fp, "          -b, --barcode      Use adapters to demultiplex reads to trimmed file(s) and an untrimmed file (no)\n");
	fprintf(fp, "          -o, --output <str>   Base name of output file; ('<reads>.trimmed-Q<int>L<int>')\n");
	fprintf(fp, "          -z, --compress       Compress output in GZIP format (no)\n");
	fprintf(fp, "          -1, --stdout         Redirect output to STDOUT, suppressing -b, -o, and -z options (no)\n");
	fprintf(fp, "          --quiet              No progress update (not quiet)\n");
	fprintf(fp, " Miscellaneous:\n");
	fprintf(fp, "          -t, --threads <int>    Number of concurrent threads [1, 16]; (1)\n");
	fprintf(fp, "\nEXAMPLES:\n");
	fprintf(fp, "          %s -Q 9 -t 2 -x adapters.fa sample.fastq -o trimmed\n", program);
	fprintf(fp, "          %s -x %s -q 3 sample-pair1.fq.gz sample-pair2.fq.gz\n", program, ILLUMINA_ADAPTER_PREFIX);
	fprintf(fp, "          %s -x %s -l 16 -L 30 -d 0 srna.fastq\n", program, ILLUMINA_SRNA_ADAPTER);
	fprintf(fp, "          %s -m mp lmp-pair1.fastq lmp-pair2.fastq\n", program);
}

void cParameter::PrintSimpleUsage(char * program, FILE * fp)
{
	fprintf(fp, "Usage: %s [options] <file> [file2]\n", program);
	fprintf(fp, "Try `%s --help' for more information.\n", program);
}

void cParameter::printCommandLine(FILE * fp)
{
	int i;
	if(argc <= 0) return;
	fprintf(fp, "COMMAND LINE:\t%s", argv[0]);
	for(i=1; i<argc; i++){
		fprintf(fp, " %s", argv[i]);
	}
	fprintf(fp, "\n");
}

void cParameter::printRelatedFiles(FILE * fp)
{
	int i;
	if(nFileCnt < 2){
		if(bStdin)
			fprintf(fp, "Input file:\tSTDIN\n");
		else
			fprintf(fp, "Input file:\t%s\n", input[0]);
		if(bStdout){
			fprintf(fp, "trimmed:\tSTDOUT\n");
		}
		else{
			fprintf(fp, "trimmed:\t");
			for(i=0; i<int(output.size()); i++){
				if(i > 0)
					fprintf(fp, "; ");
				fprintf(fp, "%s", output[i].c_str());
			}
			fprintf(fp, "\n");
			if(!untrimmed.empty()){
				fprintf(fp, "untrimmed:\t%s\n", untrimmed.c_str());
			}
		}
	}
	else{
		fprintf(fp, "Input file:\t%s\n", input[0]);
		fprintf(fp, "Paired file:\t%s\n", input[1]);
		fprintf(fp, "trimmed:\t");
		for(i=0; i<int(output.size()); i++){
			if(i > 0)
				fprintf(fp, "; ");
			fprintf(fp, "%s, %s", output[i].c_str(), output2[i].c_str());
		}
		fprintf(fp, "\n");
		if(!untrimmed.empty()){
			fprintf(fp, "untrimmed:\t%s, %s\n", untrimmed.c_str(), untrimmed2.c_str());
		}
	}
}

void cParameter::printOpt(FILE * fp, bool bLeadingRtn)
{
	int color = (fp == stdout) ? 3 : -1;
	if(bLeadingRtn)	fprintf(fp, "\n");
	fprintf(fp, "Parameters used:\n");

	// adapter
	const char * endInfo = ( ((trimMode & TRIM_ANY) == TRIM_HEAD) ? "5' end" : "3' end" );
	if(bXFile){
		fprintf(fp, "-- %s adapter sequences in file (-x):", endInfo);
		fprintf(fp, "\t%s\n", x_str.c_str());
		color_fprintf_sequences(color, fp, adapters);
	}
	else{
		fprintf(fp, "-- %s adapter sequence (-x):", endInfo);
		color_fprintf(color, fp, "\t%s\n", x_str.c_str());
	}
	if(nFileCnt == 2){
		if(!bShareAdapter){
			if(bYFile){
				fprintf(fp, "-- paired %s adapter sequences in file (-y):", endInfo);
				fprintf(fp, "\t%s\n", y_str.c_str());
				color_fprintf_sequences(color, fp, adapters2);
			}
			else{
				fprintf(fp, "-- paired %s adapter sequence (-y):", endInfo);
				color_fprintf(color, fp, "\t%s\n", y_str.c_str());
			}
		}
		if(trimMode == TRIM_MP){
			if(bJFile){
				fprintf(fp, "-- junction adapter sequences in file (-j):");
				fprintf(fp, "\t%s\n", j_str.c_str());
				color_fprintf_sequences(color, fp, juncAdapters);
			}
			else{
				fprintf(fp, "-- junction adapter sequence (-j):");
				color_fprintf(color, fp, "\t%s\n", j_str.c_str());
			}
		}
	}

	// penalty
	fprintf(fp, "-- maximum error ratio allowed (-r):\t%.3f\n", epsilon);
	fprintf(fp, "-- maximum indel error ratio allowed (-d):\t%.3f\n", delta);

	// filtering
	if(minAverageQual > 0){
		fprintf(fp, "-- mean quality threshold (-Q):\t\t%d\n", minAverageQual);
	}
	if(minEndQual > 0){
		fprintf(fp, "-- end quality threshold (-q):\t\t%d\n", minEndQual);
	}
	fprintf(fp, "-- minimum read length allowed after trimming (-l):\t%d\n", minLen);
	if(maxLen > 0){
		fprintf(fp, "-- maximum read length for output (-L):\t%d\n", maxLen);
	}

	// input
	fprintf(fp, "-- file format (-f):\t\t%s %s\n", FASTQ_FORMAT_NAME[fastqFormat], (bAutoFormat ? "(auto detected)" : ""));

	// misc
	if(nFileCnt < 2){
		fprintf(fp, "-- minimum overlap length for adapter detection (-k):\t%d\n", minK);
	}
	else{
		if(trimMode == TRIM_MP){
			fprintf(fp, "-- minimum overlap length for junction adapter detection (-k):\t%d\n", minK);
		}
	}
	if(nThreads > 1){
		fprintf(fp, "-- number of concurrent threads (-t):\t%d\n", nThreads);
	}
}

int cParameter::GetOpt(int argc, char *argv[], char * errMsg)
{
	const char *options = "x:y:j:m:r:d:q:l:L:nuf:bo:z1Q:k:t:*vh";
	OPTION_ITEM longOptions[] = {
		{"barcode", 'b'},
		{"mode", 'm'},
		{"end-quality", 'q'},
		{"mean-quality", 'Q'},
		{"min", 'l'},
		{"max", 'L'},
		{"output", 'o'},
		{"threads", 't'},
		{"format", 'f'},
		{"stdout", '1'},
		{"compress", 'z'},
		{"quiet", '*'},
		{"version", 'v'},
		{"help", 'h'}
	};
	char trimmed[MAX_PATH+1+100];
	this->argc = argc;
	this->argv = argv;
	bool bSetX, bSetY, bSetJ, bSetO, bSetL, bSetD, bSetK;
	bSetX = bSetY = bSetJ = bSetO = bSetL = bSetD = bSetK = false;
	int iRet = 0;
	int i, j;
	char chr;
	const char * str;
	for(i=1; i<argc; i++){
		chr = argv[i][0];
		if(chr != '-'){ // input file
			if(nFileCnt < 2){
				strcpy(input[nFileCnt], argv[i]);
				nFileCnt++;
			}
			continue;
		}
		chr = argv[i][1];
		if(chr == '-'){ // "--" for long options
			str = argv[i] + 2;
			for(j=0; j<int(sizeof(longOptions)/sizeof(OPTION_ITEM)); j++){
				if(strcmp(str, longOptions[j].name) == 0){
					chr = longOptions[j].chr;
					break;
				}
			}
			if(j >= int(sizeof(longOptions)/sizeof(OPTION_ITEM))){
				sprintf(errMsg, "No such an option --%s", str);
				iRet = -2;
				break;
			}
		}
		else if(chr == '\0'){
			bStdin = true;
			continue;
		}
		str = strchr(options, chr);
		if(str == NULL){
			sprintf(errMsg, "No such an option -%c", chr);
			iRet = -2;
			break;
		}
		if(str[1] == ':'){ // has an argument
			if(++i >= argc){
				sprintf(errMsg, "-%c needs an argument", chr);
				iRet = -2;
				break;
			}
		}
		switch(chr){
		case 'x':
			x_str.assign(argv[i]);
			str = strchr(x_str.c_str(), '.');
			bXFile = (str != NULL);
			bSetX = true;
			break;
		case 'y':
			y_str.assign(argv[i]);
			str = strchr(y_str.c_str(), '.');
			bYFile = (str != NULL);
			bSetY = true;
			break;
		case 'j':
			j_str.assign(argv[i]);
			str = strchr(j_str.c_str(), '.');
			bJFile = (str != NULL);
			bSetJ = true;
			break;
		case 'm':
			if(strcasecmp(argv[i], "head") == 0){
				trimMode = TRIM_HEAD;
			}
			else if(strcasecmp(argv[i], "any") == 0){
				trimMode = TRIM_ANY;
			}
			else if(strcasecmp(argv[i], "mp") == 0){
				trimMode = TRIM_MP;
			}
			else if(strcasecmp(argv[i], "tail") == 0){
				trimMode = TRIM_TAIL;
			}
			else if(strcasecmp(argv[i], "pe") == 0){
				trimMode = TRIM_PE;
			}
			else{
				iRet = -3;
			}
			break;
		case 'r':
			if( (argv[i][0] < '0' || argv[i][0] > '9') && argv[i][0] != '.' ){
				iRet = -3;
				break;
			}
			epsilon = atof(argv[i]);
			if(epsilon < 0) epsilon = 0;
			else if(epsilon > 2) epsilon = 1 / epsilon;
			else if(epsilon > 0.5) epsilon = 0.5;
			break;
		case 'd':
			if( (argv[i][0] < '0' || argv[i][0] > '9') && argv[i][0] != '.' ){
				iRet = -3;
				break;
			}
			delta = atof(argv[i]);
			bSetD = true;
			break;
		case 'q':
			if(argv[i][0] < '0' || argv[i][0] > '9'){
				iRet = -3;
				break;
			}
			minEndQual = atoi(argv[i]);
			break;
		case 'Q':
			if(argv[i][0] < '0' || argv[i][0] > '9'){
				iRet = -3;
				break;
			}
			minAverageQual = atoi(argv[i]);
			break;
		case 'l':
			if(argv[i][0] < '0' || argv[i][0] > '9'){
				iRet = -3;
				break;
			}
			minLen = atoi(argv[i]);
			if(minLen < 0) minLen = 0;
			break;
		case 'L':
			if(argv[i][0] < '0' || argv[i][0] > '9'){
				iRet = -3;
				break;
			}
			maxLen = atoi(argv[i]);
			bSetL = true;
			break;
		case 'n':
			bFilterNs = true;
			break;
		case 'u':
			bFilterUndetermined = true;
			break;
		case 'f':
			if(strcasecmp(argv[i], "auto") == 0){
				break;
			}
			bAutoFormat = false;
			if(strcasecmp(argv[i], "sanger") == 0){
				fastqFormat = SANGER_FASTQ;
				baseQual = 33;
			}
			else if(strcasecmp(argv[i], "solexa") == 0){
				fastqFormat = SOLEXA_FASTQ;
				baseQual = 64;
			}
			else{
				sprintf(errMsg, "unkown format '%s', please select either 'sanger' or 'solexa'", argv[i]);
				iRet = -2;
			}
			break;
		case 'o':
			gzstrncpy(trimmed, argv[i], MAX_PATH);
			bSetO = true;
			break;
		case 'z':
			outputFormat = COMPRESS_GZ;
			break;
		case '1':
			bStdout = true;
			break;
		case '*':
			bQuiet = true;
			break;
		case 'k':
			minK = atoi(argv[i]);
			if(minK < 1) minK = 1;
			bSetK = true;
			break;
		case 't':
			if(argv[i][0] < '0' || argv[i][0] > '9'){
				iRet = -3;
				break;
			}
			nThreads = atoi(argv[i]);
			if(nThreads < 1) nThreads = 1;
			else if(nThreads > 16) nThreads = 16;
			break;
		default:
			iRet = -1;
			bEnquireVersion |= (chr == 'v');
			break;
		}
		if(iRet < 0) break;
	}
	if(iRet < 0){
		if(iRet == -3){
			sprintf(errMsg, "Invalid argument of -%c: %s", chr, argv[i]);
		}
		return iRet;
	}

	// input and output
	if(nFileCnt == 0){
		if(!bStdin){
			sprintf(errMsg, "No input file specified");
			return -2;
		}
	}
	else{
		if(bStdin){
			sprintf(errMsg, "STDIN can not be specified with other input files");
			return -2;
		}
	}
	if(bStdout){
		if(nFileCnt == 2){
			sprintf(errMsg, "STDOUT can not be specified for output paired files");
			return -2;
		}
		if(bBarcode){
			sprintf(errMsg, "STDOUT can not be used for demultiplexing (-b)");
			return -2;
		}
		if(bSetO){
			sprintf(errMsg, "STDOUT can not be used with -o");
			return -2;
		}
		if(outputFormat != COMPRESS_NONE){
			sprintf(errMsg, "STDOUT can not be used for compressing (-z)");
			return -2;
		}
	}
	// trimming mode
	if(nFileCnt < 2){
		if( (trimMode & TRIM_ANY) == TRIM_DEFAULT ){
			trimMode = TRIM_TAIL;
		}
	}
	else{ // nFileCnt == 2
		if( trimMode != TRIM_MP ) {
			trimMode = TRIM_PE;
		}
	}
	// adapters
	if(bSetX){ // specified by command
		if(bXFile){
			if(!ReadFasta(x_str.c_str(), adapters)){
				sprintf(errMsg, "Can not read adapter sequences from FASTA file \"%s\"", x_str.c_str());
				return -2;
			}
		}
		else{
			if(int(x_str.length()) > MAX_ADAPTER_LEN){
				string tmpString;
				if( (trimMode & TRIM_ANY) == TRIM_HEAD )
					tmpString.assign(x_str.c_str() + x_str.length() - MAX_ADAPTER_LEN);
				else
					tmpString.assign(x_str.c_str(), 0, MAX_ADAPTER_LEN);
				x_str.assign(tmpString);
			}
			adapters.push_back(x_str);
		}
	}
	else{ // default
		x_str.assign(ILLUMINA_PAIR1_ADAPTER_PREFIX);
		adapters.push_back(x_str);
	}
	if(nFileCnt == 2){
		if(bSetY){ // specified by command
			if(bYFile){
				if(!ReadFasta(y_str.c_str(), adapters2)){
					sprintf(errMsg, "Can not read adapter sequences from FASTA file \"%s\"", y_str.c_str());
					return -2;
				}
			}
			else{
				if(int(y_str.length()) > MAX_ADAPTER_LEN){
					string tmpString;
					if( (trimMode & TRIM_ANY) == TRIM_HEAD )
						tmpString.assign(y_str.c_str() + y_str.length() - MAX_ADAPTER_LEN);
					else
						tmpString.assign(y_str.c_str(), 0, MAX_ADAPTER_LEN);
					y_str.assign(tmpString);
				}
				adapters2.push_back(y_str);
			}
		}
		else{ // default
			if(bSetX){
				bShareAdapter = true;
			}
			else{
				y_str.assign(ILLUMINA_PAIR2_ADAPTER_PREFIX);
				adapters2.push_back(y_str);
			}
		}
		if(trimMode == TRIM_MP){
			if(bSetJ){ // specified by command
				if(bXFile){
					if(!ReadFasta(j_str.c_str(), juncAdapters)){
						sprintf(errMsg, "Can not read adapter sequences from FASTA file \"%s\"", j_str.c_str());
						return -2;
					}
				}
				else{
					if(int(j_str.length()) > MAX_ADAPTER_LEN){
						string tmpString;
						if( (trimMode & TRIM_ANY) == TRIM_HEAD )
							tmpString.assign(j_str.c_str() + j_str.length() - MAX_ADAPTER_LEN);
						else
							tmpString.assign(j_str.c_str(), 0, MAX_ADAPTER_LEN);
						j_str.assign(tmpString);
					}
					juncAdapters.push_back(j_str);
				}
			}
			else{ // default
				j_str.assign(ILLUMINA_JUNCTION_ADAPTER);
				juncAdapters.push_back(j_str);
			}
		}
	}
	// penalty
	if(bSetD){
		if(delta < 0) delta = 0;
		else{
			if(delta > 2){
				delta = 1 / delta;
			}
			if(delta > epsilon){
				delta = epsilon;
			}
		}
	}
	if(!bSetK){
		minK = (trimMode == TRIM_MP) ? (j_str.length() / 2) : max(int(4 - 10 * epsilon), 1);
	}
	// filtering
	if(bSetL){
		if(maxLen < minLen){
			maxLen = minLen;
		}
	}

	// prepare output files
	if(bStdout){
		return iRet;
	}
	char * end;
	if(!bSetO){
		if(bStdin){
			if(minAverageQual > 0)
				sprintf(trimmed, "trimmed-Q%dL%d", minAverageQual, minLen);
			else
				sprintf(trimmed, "trimmed-L%d", minLen);
		}
		else{
			gzstrncpy(trimmed, input[0], MAX_PATH);
			end = occOfLastDot(trimmed);
			if(minAverageQual > 0)
				sprintf(end, ".trimmed-Q%dL%d", minAverageQual, minLen);
			else
				sprintf(end, ".trimmed-L%d", minLen);
		}
	}
	strcpy(logfile, trimmed);
	strcat(logfile, ".log");

	string fileName, fileName2;
	if(bBarcode){
		char buffer[MAX_PATH];
		if(nFileCnt >= 2){
			for(i=0; i<int(adapters.size()); i++){
				if(bShareAdapter){
					for(j=0; j<int(adapters.size()); j++){
						sprintf(buffer, "%c%02d", ('A'+i), (j+1));
						barcodes.push_back(string(buffer));
						fileName.assign(string(trimmed) + string(buffer) + string("-pair1.fastq"));
						fileName2.assign(string(trimmed) + string(buffer) + string("-pair2.fastq"));
						if(outputFormat == COMPRESS_GZ){
							fileName += string(".gz");
							fileName2 += string(".gz");
						}
						output.push_back(fileName);
						output2.push_back(fileName2);
					}
				}
				else{
					for(j=0; j<int(adapters2.size()); j++){
						sprintf(buffer, "%c%02d", ('A'+i), (j+1));
						barcodes.push_back(string(buffer));
						fileName.assign(string(trimmed) + string(buffer) + string("-pair1.fastq"));
						fileName2.assign(string(trimmed) + string(buffer) + string("-pair2.fastq"));
						if(outputFormat == COMPRESS_GZ){
							fileName += string(".gz");
							fileName2 += string(".gz");
						}
						output.push_back(fileName);
						output2.push_back(fileName2);
					}
				}
			}
			untrimmed.assign(string(trimmed) + string("untrimmed-pair1.fastq"));
			untrimmed2.assign(string(trimmed) + string("untrimmed-pair2.fastq"));
			if(outputFormat == COMPRESS_GZ){
				untrimmed += string(".gz");
				untrimmed2 += string(".gz");
			}
		}
		else{
			for(i=0; i<int(adapters.size()); i++){
				sprintf(buffer, "A%02d.fastq", (i+1));
				barcodes.push_back(string(buffer));
				fileName.assign(string(trimmed) + string(buffer));
				if(outputFormat == COMPRESS_GZ){
					fileName += string(".gz");
				}
				output.push_back(fileName);
			}
			untrimmed.assign(string(trimmed) + string("untrimmed.fastq"));
			if(outputFormat == COMPRESS_GZ){
				untrimmed += string(".gz");
			}
		}
	}
	else{
		if(nFileCnt >= 2){
			fileName.assign(string(trimmed) + string("-pair1.fastq"));
			fileName2.assign(string(trimmed) + string("-pair2.fastq"));
			if(outputFormat == COMPRESS_GZ){
				fileName += string(".gz");
				fileName2 += string(".gz");
			}
			output.push_back(fileName);
			output2.push_back(fileName2);
		}
		else{
			fileName.assign(string(trimmed) + string(".fastq"));
			if(outputFormat == COMPRESS_GZ){
				fileName += string(".gz");
			}
			output.push_back(fileName);
		}
	}

	return iRet;
}

