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
#ifndef _FASTQ_H
#define _FASTQ_H

// standard libs
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <getopt.h>
#include <assert.h>
#include <math.h>
#include <sys/stat.h>
#include <search.h>
#include <limits.h>

#include <time.h>
#ifdef __MACH__
#include <mach/mach_time.h>
#define ORWL_NANO (+1.0E-9)
#define ORWL_GIGA UINT64_C(1000000000)
#endif

#include "common.h"

typedef struct tag_LINE {
	char *s;
	int n;
	size_t a;
}LINE;

typedef enum{
	TAG_NORMAL = 0,
	TAG_BLURRY = 1,
	TAG_BADQUAL = 2,
	TAG_EMPTY = 3,
	TAG_SHORT = 4,
	TAG_CONTAMINANT = 5,
	TAG_UNDETERMINED = 6,
	TAG_LONG = 7
}REC_TAG;

typedef struct tag_REC{
	uint32 tag:7; // REC_TAG
	uint32 bExchange:1;
	uint32 nCnt:24;
	INDEX idx;
	LINE id;
	LINE seq;
	LINE com;
	LINE qual;
}RECORD;

class cFQ
{
private:
//	int size;

public:
	char tag;

	RECORD rec;
//	RECORD *pBuffer;
//	int cnt;

	FILE * in;
	int64 offset;
	int64 next_pos;
	int rno;

private:
	inline int read_line(LINE &l);

public:
	cFQ();
	~cFQ();
//	bool InitBuffer(int nBuffSize=256);
//	void DestroyBuffer();
//	void clearBuffer();
	void associateFile(FILE * fp);
	int readRecord(RECORD *pRecord=NULL);
//	int readRecord2Buffer();
//	RECORD * getLastRecord();
	inline int64 tell(){ return offset; }
};

typedef struct tag_CFILE{
	FILE * fp;
	bool bGz;
}CFILE;

enum FASTQ_FORMAT{
	SANGER_FASTQ = 0,
	SOLEXA_FASTQ = 1,
	UNKNOWN_FASTQ = 2,
	CONTRADICT_FASTQ = 3,
	FASTA = 4,
	FASTQ_FORMAT_CNT = 5
};
extern const char * FASTQ_FORMAT_NAME[FASTQ_FORMAT_CNT];

// open a file, possibly gzipped, exit on failure
extern CFILE gzopen(const char * fileName, const char *mode);
extern int gzclose(CFILE *f);
extern int64 gzsize(const char * fileName);
extern enum FASTQ_FORMAT gzformat(char * fileNames[], int nFileCnt);
extern int gzreadlen(char * fileName);
extern void gzstrncpy (char * dest, const char * src, int n);

// time function
extern void versatile_gettime(struct timespec *tp);

#endif // _FASTQ_H
