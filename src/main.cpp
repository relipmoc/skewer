/**********************************************************************
 * Skewer - a fast and accurate adapter trimming tool
 *          using the bit-masked k-difference matching algorithm
 * Copyright (c) 2013-2016 by Hongshan Jiang
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
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <pthread.h>
#include <unistd.h>
#include <assert.h>

#include "common.h"
#include "fastq.h"
#include "parameter.h"
#include "matrix.h"

using namespace std;

inline void OutputTaggedRecord(FILE * fpOut, RECORD * pRecord)
{
	// refer to "enum REC_TAG" defined in "fastq.h"
	const char * TAG_NAME[8] = { "NORMAL", "BLURRY", "BADQUAL", "EMPTY", "SHORT", "CONTAMINANT", "UNDETERMINED", "LONG" };
	if(pRecord->com.n > 0){ // fastq
		fprintf(fpOut, "@%s TAG=%s\n%s\n+\n%s\n", strtok(pRecord->id.s, "\n"), TAG_NAME[pRecord->tag], pRecord->seq.s, pRecord->qual.s);
	}
	else{ // fasta
		fprintf(fpOut, ">%s TAG=%s\n%s\n", strtok(pRecord->id.s, "\n"), TAG_NAME[pRecord->tag], pRecord->seq.s);
	}
}

inline void OutputMaskedRecord(FILE * fpOut, RECORD * pRecord, int offset, int len)
{
	int i;
	for(i=0; i<offset; i++){
		pRecord->seq.s[i] = tolower(pRecord->seq.s[i]);
	}
	for(i=offset+len; i<pRecord->seq.n; i++){
		pRecord->seq.s[i] = tolower(pRecord->seq.s[i]);
	}
	if(pRecord->com.n > 0) // fastq
		fprintf(fpOut, "@%s%s\n+\n%s\n", pRecord->id.s, pRecord->seq.s, pRecord->qual.s);
	else // fasta
		fprintf(fpOut, ">%s%s\n", pRecord->id.s, pRecord->seq.s);
}

inline void OutputEntireRecord(FILE * fpOut, RECORD * pRecord)
{
	if(pRecord->com.n > 0) // fastq
		fprintf(fpOut, "@%s%s\n+\n%s\n", pRecord->id.s, pRecord->seq.s, pRecord->qual.s);
	else // fasta
		fprintf(fpOut, ">%s%s\n", pRecord->id.s, pRecord->seq.s);
}

inline void OutputEntireRecordFilledWithNs(FILE * fpOut, RECORD * pRecord, int offset, int len)
{
	int len2 = pRecord->seq.n - offset - len;
	string s1 = string(offset, 'N');
	string s2 = string(len2, 'N');
	if(pRecord->com.n > 0){ // fastq
		string q1 = string(offset, '!');
		string q2 = string(len2, '!');
		fprintf(fpOut, "@%s%s%.*s%s\n+\n%s%.*s%s\n", pRecord->id.s, s1.c_str(), len, pRecord->seq.s + offset, s2.c_str(),
				q1.c_str(), len,  pRecord->qual.s + offset, q2.c_str());
	}
	else // fasta
		fprintf(fpOut, ">%s%s%.*s%s\n", pRecord->id.s, s1.c_str(), len, pRecord->seq.s + offset, s2.c_str());
}

inline void OutputPartialRecord(FILE * fpOut, RECORD * pRecord, int offset, int len)
{
	if(pRecord->com.n > 0) // fastq
		fprintf(fpOut, "@%s%.*s\n+\n%.*s\n", pRecord->id.s, len, pRecord->seq.s + offset, len, pRecord->qual.s + offset);
	else // fasta
		fprintf(fpOut, ">%s%.*s\n", pRecord->id.s, len, pRecord->seq.s + offset);
}

class cStats
{
	struct timespec tpstart, tpend;

private:
	bool bPaired;
	size_t minLen;
	size_t maxLen;
	size_t maxReadLen;
	size_t allocLen;
	size_t nBarcodes;
	long * pHist;
	long * pBarcode;
	vector<string> * pBarcodeNames;
	const char * pDecorate;

public:
	long nBlurry;
	long nBad;
	long nContaminant;
	long nUndetermined;
	long nEmpty;
	long nShort;
	long nLong;
	long nTrimAvail;
	long nUntrimAvail;

	CFILE *fpOuts;
	CFILE *fpOuts2;
	CFILE *fpMasked;
	CFILE *fpMasked2;
	CFILE *fpExcluded;
	CFILE *fpExcluded2;
	CFILE fpUntrim;
	CFILE fpUntrim2;
	CFILE fpBarcodes;
	CFILE fpMapfile;
	int nFiles;
	int nFiles2;
	bool bBarcode;
	bool bStdout;

	// for mutiple threads
	int64 total_file_length;
	cFQ * pfq;
	cFQ * pfq2;
	FILE * fpOut;
	FILE * fpOut2;
	FILE * fpMask;
	FILE * fpMask2;
	FILE * fpExcl;
	FILE * fpExcl2;
	FILE * fpBarcode;
	int minAverageQual;
	int minEndQual;
	bool bFivePrimeEnd;
	bool bFilterNs;
	bool bFilterUndetermined;
	bool bRedistribute;
	bool bQuiet;
	bool bMatepair;
	bool bFillWithNs;
	bool bCutTail;
	int iCutF, iCutR;

	int getMinLen(){
		return int(minLen);
	}
	int getMaxLen(){
		return int(maxLen);
	}
public:
	cStats(){
		nBlurry = nBad = nContaminant = nUndetermined = nEmpty = nShort = nLong = 0L;
		nTrimAvail = nUntrimAvail = 0;
		bPaired = false;
		pHist = NULL;
		pBarcode = NULL;
		bBarcode = false;
		bStdout = false;
		bFilterNs = false;
		bFilterUndetermined = false;
		bFillWithNs = false;
		bRedistribute = false;
		bCutTail = false;
		minLen = allocLen = 0;
		maxLen = INT_MAX;
		nBarcodes = 0;
		iCutF = iCutR = 0;
		pDecorate = "";
		
		fpOuts = fpOuts2 = NULL;
		fpMasked = fpMasked2 = NULL;
		fpMask = fpMask2 = NULL;
		fpExcl = fpExcl2 = NULL;
		fpExcluded = fpExcluded2 = NULL;
		nFiles = nFiles2 = 0;
		fpUntrim.fp = fpUntrim2.fp = NULL;
		fpBarcodes.fp = NULL;
		fpMapfile.fp = NULL;
	}
	~cStats(){
		int i;
		gzclose(&fpMapfile);
		gzclose(&fpBarcodes);
		gzclose(&fpUntrim2);
		gzclose(&fpUntrim);
		if(fpOuts2 != NULL){
			for(i=nFiles2-1; i>=0; i--){
				gzclose(&fpOuts2[i]);
			}
			free(fpOuts2);
			fpOuts2 = NULL;
			nFiles2 = 0;
		}
		if(fpOuts != NULL){
			for(i=nFiles-1; i>=0; i--){
				gzclose(&fpOuts[i]);
			}
			free(fpOuts);
			fpOuts = NULL;
			nFiles = 0;
		}
		if (fpMasked != NULL){
			gzclose(&fpMasked[0]);
			free(fpMasked);
			fpMasked = NULL;
		}
		if (fpMasked2 != NULL){
			gzclose(&fpMasked2[0]);
			free(fpMasked2);
			fpMasked2 = NULL;
		}
		if (fpExcluded != NULL){
			gzclose(&fpExcluded[0]);
			free(fpExcluded);
			fpExcluded = NULL;
		}
		if (fpExcluded2 != NULL) {
			gzclose(&fpExcluded2[0]);
			free(fpExcluded2);
			fpExcluded2 = NULL;
		}
		if(pBarcode != NULL){
			delete [] pBarcode;
			pBarcode = NULL;
			nBarcodes = 0;
		}
		if(pHist != NULL){
			delete [] pHist;
			pHist = NULL;
			allocLen = 0;
		}
	}
	bool initHist(cParameter * pParameter){
		this->minLen = pParameter->minLen;
		this->maxReadLen = 0;
		pHist = new long[50];
		if(pHist == NULL)
			return false;
		allocLen = 50;
		memset(pHist, 0, allocLen * sizeof(long));
		nBarcodes = 0;
		this->iCutF = pParameter->iCutF;
		this->iCutR = pParameter->iCutR;
		if(!pParameter->bBarcode)
			return true;
		pBarcode = new long[pParameter->output.size()];
		if(pBarcode == NULL)
			return false;
		nBarcodes = pParameter->output.size();
		memset(pBarcode, 0, nBarcodes * sizeof(long));
		pBarcodeNames = &pParameter->barcodeNames;
		bBarcode = true;
		return true;
	}
	void InitGlobalAttributes(cParameter * pParameter, int64 total_file_length, bool bPaired, cFQ * pFq, cFQ * pFq2){
		// global attributes used by threads
		this->total_file_length = total_file_length;
		this->pfq = pFq;
		this->fpOut = pParameter->bStdout ? stdout : fpOuts[0].fp;
		this->pDecorate = pParameter->pDecorate;
		if(bPaired){
			this->pfq2 = pFq2;
			this->fpOut2 = fpOuts2[0].fp;
			this->fpBarcode = fpBarcodes.fp;
		}
		this->fpMask = (fpMasked != NULL) ? fpMasked[0].fp : NULL;
		this->fpMask2 = (fpMasked2 != NULL) ? fpMasked2[0].fp : NULL;
		this->fpExcl = (fpExcluded != NULL) ? fpExcluded[0].fp : NULL;
		this->fpExcl2 = (fpExcluded2 != NULL) ? fpExcluded2[0].fp : NULL;
		this->minAverageQual = (pParameter->minAverageQual > 0) ? (pParameter->baseQual + pParameter->minAverageQual) : 0;
		this->minEndQual = (pParameter->minEndQual > 0) ? (pParameter->baseQual + pParameter->minEndQual) : 0;
		this->minLen = pParameter->minLen;
		this->maxLen = (pParameter->maxLen > 0) ? pParameter->maxLen : INT_MAX;
		this->bFivePrimeEnd = ((pParameter->trimMode & TRIM_ANY) == TRIM_HEAD);
		this->bQuiet = pParameter->bQuiet || pParameter->bStdin;
		this->bFilterNs = pParameter->bFilterNs;
		this->bFilterUndetermined = pParameter->bFilterUndetermined;
		this->bRedistribute = pParameter->bRedistribute;
		this->bCutTail = pParameter->bCutTail;
		this->bFillWithNs = pParameter->bFillWithNs;
	}
	bool openOutputFiles(cParameter * pParameter){
		bPaired = (pParameter->nFileCnt >= 2);
		bStdout = pParameter->bStdout;
		assert(!(bStdout && (bPaired | bBarcode)));
		if(bStdout){
			return true;
		}
		fpOuts = (CFILE *)calloc(pParameter->output.size(), sizeof(CFILE));
		if(fpOuts == NULL){
			fprintf(stderr, "Can not allocate memory for file handles for writing\n");
			return false;
		}
		for(nFiles=0; nFiles<int(pParameter->output.size()); nFiles++){
			fpOuts[nFiles] = gzopen(pParameter->output[nFiles].c_str(), "w");
			if(fpOuts[nFiles].fp == NULL){
				fprintf(stderr, "Can not open %s for writing\n", pParameter->output[nFiles].c_str());
				break;
			}
		}
		if(!pParameter->untrimmed.empty()){
			fpUntrim = gzopen(pParameter->untrimmed.c_str(), "w");
			if(fpUntrim.fp == NULL){
				fprintf(stderr, "Can not open %s for writing\n", pParameter->untrimmed.c_str());
				return false;
			}
		}
		else if(pParameter->bWriteMasked) {
			fpMasked = (CFILE *)calloc(1, sizeof(CFILE));
			if(fpMasked == NULL){
				fprintf(stderr, "Can not allocate memory for file handles for writing\n");
				return false;
			}
			fpMasked[0] = gzopen(pParameter->masked[0].c_str(), "w");
			if(fpMasked[0].fp == NULL){
				fprintf(stderr, "Can not open masked for writing\n");
			}
			if (bPaired) {
				fpMasked2 = (CFILE *)calloc(1, sizeof(CFILE));
				if(fpMasked2 == NULL){
					fprintf(stderr, "Can not allocate memory for file handles for writing\n");
					return false;
				}
				fpMasked2[0] = gzopen(pParameter->masked2[0].c_str(), "w");
				if(fpMasked2[0].fp == NULL){
					fprintf(stderr, "Can not open masked for writing\n");
				}
			}
		}
		if(bPaired){
			fpOuts2 = (CFILE *)calloc(pParameter->output2.size(), sizeof(CFILE));
			if(fpOuts2 == NULL){
				fprintf(stderr, "Can not allocate memory for file handles for writing\n");
				return false;
			}
			for(nFiles2=0; nFiles2<int(pParameter->output2.size()); nFiles2++){
				fpOuts2[nFiles2] = gzopen(pParameter->output2[nFiles2].c_str(), "w");
				if(fpOuts2[nFiles2].fp == NULL){
					fprintf(stderr, "Can not open %s for writing\n", pParameter->output2[nFiles2].c_str());
					break;
				}
			}
			if(!pParameter->untrimmed2.empty()){
				fpUntrim2 = gzopen(pParameter->untrimmed2.c_str(), "w");
				if(fpUntrim2.fp == NULL){
					fprintf(stderr, "Can not open %s for writing\n", pParameter->untrimmed2.c_str());
					return false;
				}
			}
			if(!pParameter->barcodes.empty()){
				fpBarcodes = gzopen(pParameter->barcodes.c_str(), "w");
				if(fpBarcodes.fp == NULL){
					fprintf(stderr, "Can not open %s for writing\n", pParameter->barcodes.c_str());
					return false;
				}
			}
			if(!pParameter->mapfile.empty()){
				fpMapfile = gzopen(pParameter->mapfile.c_str(), "w");
				if(fpMapfile.fp == NULL){
					fprintf(stderr, "Can not open %s for writing\n", pParameter->mapfile.c_str());
					return false;
				}
			}
		}
		if(pParameter->bWriteExcluded) {
			fpExcluded = (CFILE *)calloc(1, sizeof(CFILE));
			if(fpExcluded == NULL){
				fprintf(stderr, "Can not allocate memory for file handles for writing\n");
				return false;
			}
			fpExcluded[0] = gzopen(pParameter->excluded[0].c_str(), "w");
			if(fpExcluded[0].fp == NULL){
				fprintf(stderr, "Can not open excluded for writing\n");
			}
			if (bPaired) {
				fpExcluded2 = (CFILE *)calloc(1, sizeof(CFILE));
				if(fpExcluded2 == NULL){
					fprintf(stderr, "Can not allocate memory for file handles for writing\n");
					return false;
				}
				fpExcluded2[0] = gzopen(pParameter->excluded2[0].c_str(), "w");
				if(fpExcluded2[0].fp == NULL){
					fprintf(stderr, "Can not open excluded for writing\n");
				}
			}
		}
		return true;
	}
	bool incrementCount(size_t readLen){
		if(readLen + 1 > allocLen){
			size_t newAllocLen = readLen * 3 / 2 + 32;
			long * pNewHist = new long[newAllocLen];
			if(pNewHist == NULL)
				return false;
			if(pHist != NULL){
				memcpy(pNewHist, pHist, allocLen * sizeof(long));
				delete [] pHist;
			}
			memset(&pNewHist[allocLen], 0, (newAllocLen - allocLen) * sizeof(long));
			pHist = pNewHist;
			allocLen = newAllocLen;
		}
		if(readLen > maxReadLen){
			maxReadLen = readLen;
		}
		pHist[readLen]++;
		return true;
	}
	bool incrementBarcode(size_t bc){
		if(bc >= nBarcodes)
			return false;
		pBarcode[bc]++;
		return true;
	}
	void printHist(FILE * fp, bool bLeadingRtn=true){
		char buffer[100];
		long sum = nTrimAvail + nUntrimAvail;
		sprintf(buffer, "%ld", sum);
		int width = int(strlen(buffer));
		if(bLeadingRtn)
			fprintf(fp, "\n");
		int i;
		if(bBarcode && (nTrimAvail > 0) ){
			fprintf(fp, "Barcode dispatch after trimming:\n");
			fprintf(fp, "category \tcount\tpercentage:\n");
			for(i=0; i<int(nBarcodes); i++){
				fprintf(fp, "%s\t%*ld\t%6.2f%%\n", (*pBarcodeNames)[i].c_str(), width, pBarcode[i], pBarcode[i] * 100.0 / nTrimAvail);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "Length distribution of reads after trimming:\n");
		fprintf(fp, "length\tcount\tpercentage\n");
		for(i=0; (i<=int(maxReadLen)) && (pHist[i] == 0); i++);
		for( ; i<=int(maxReadLen); i++){
			fprintf(fp, "%3d\t%*ld\t%6.2f%%\n", i, width, pHist[i], pHist[i] * 100.0 / sum);
		}
	}
	void progress(double ratio, int width) {
		char bar[101];

		if(width < 25) width = 25;
		else if(width > 100) width = 100;
		int i;
		int point = int(min(ratio, 1.0) * width);
		for(i=0; i<point-1; i++){
			bar[i] = '=';
		}
		bar[i++] = '>';
		for(; i<width; i++){
			bar[i] = ' ';
		}
		bar[i] = '\0';

		fprintf(stderr, "\r|%s| (%.2f%%)", bar, ratio * 100);
		fflush(stderr);
	}
	void endProgress(){
		fprintf(stderr, "\n");
	}
	char* timeStamp(){
		static const char wday_name[][4] = {
			"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"
		};
		static const char mon_name[][4] = {
			"Jan", "Feb", "Mar", "Apr", "May", "Jun",
			"Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
		};
		static char result[26];
		time_t curtime;
		struct tm * loctime;

		curtime = time(NULL);
		loctime = localtime(&curtime);
		sprintf(result, "%.3s %.3s%3d %.2d:%.2d:%.2d %d",
			wday_name[loctime->tm_wday],
			mon_name[loctime->tm_mon],
			loctime->tm_mday, loctime->tm_hour,
			loctime->tm_min, loctime->tm_sec,
			1900 + loctime->tm_year);
		return result;
	}
	void printTime(const char * message, FILE * fp, int flag=0x1){
		if(flag & 0x02) fprintf(fp, "\n");
		fputs(timeStamp(), fp);
		int color = (fp == stdout) ? 2 : -1;
		color_fprintf(color, fp, " >> %s",  message);
		if(flag & 0x01) fprintf(fp, "\n");
	}
	void start(){
		versatile_gettime(&tpstart);
	}
	void end(){
		versatile_gettime(&tpend);
	}
	void printDiffTime(FILE * fp, bool bRnt=true){
		double timediff = (tpend.tv_sec-tpstart.tv_sec)+(tpend.tv_nsec-tpstart.tv_nsec)/1e9;
		fprintf(fp, " (%.3lfs)", timediff);
		if(bRnt) fprintf(fp, "\n");
	}
	void printSummary(FILE * fp){
		char buffer[100];
		long sum = nBlurry + nBad + nContaminant + nUndetermined + nEmpty + nShort + nLong + nTrimAvail + nUntrimAvail;
		sprintf(buffer, "%ld", sum);
		int width = int(strlen(buffer));
		const char * entity = bPaired ? "read pairs" : "reads";
		fprintf(fp, "%.*ld %s processed; of these:\n", width, sum, entity);
		if(bFilterNs)
			fprintf(fp, "%*ld (%5.2f%%) degenerative %s filtered out\n", width, nBlurry, (nBlurry * 100.0) / sum, entity);
		if(nBad > 0)
			fprintf(fp, "%*ld (%5.2f%%) %s filtered out by quality control\n", width, nBad, (nBad * 100.0) / sum, entity);
		if(nContaminant > 0)
			fprintf(fp, "%*ld (%5.2f%%) non-junction %s filtered out by contaminant control\n", width, nContaminant, (nContaminant * 100.0) / sum, entity);
		if(bFilterUndetermined)
			fprintf(fp, "%*ld (%5.2f%%) undetermined %s filtered out by contaminant control\n", width, nUndetermined, (nUndetermined * 100.0) / sum, entity);
		if(minLen > 0){
			if(minLen > 1)
				fprintf(fp, "%*ld (%5.2f%%) short %s filtered out after trimming by size control\n", width, nShort, (nShort * 100.0) / sum, entity);
			fprintf(fp, "%*ld (%5.2f%%) empty %s filtered out after trimming by size control\n", width, nEmpty, (nEmpty * 100.0) / sum, entity);
		}
		if(nLong > 0)
			fprintf(fp, "%*ld (%5.2f%%) long %s filtered out after trimming by size control\n", width, nLong, (nLong * 100.0) / sum, entity);
		long nAvailSum = nTrimAvail + nUntrimAvail;
		fprintf(fp, "%*ld (%5.2f%%) %s available", width, nAvailSum, (nAvailSum * 100.0) / sum, entity);
		if(nAvailSum > 0){
			fprintf(fp, "; of these:\n");
			if(nTrimAvail > 0)
				fprintf(fp, "%*ld (%5.2f%%) %s %s available after processing\n", width, nTrimAvail, (nTrimAvail * 100.0) / nAvailSum, pDecorate, entity);
			if(nUntrimAvail > 0)
				fprintf(fp, "%*ld (%5.2f%%) un%s %s available after processing\n", width, nUntrimAvail, (nUntrimAvail * 100.0) / nAvailSum, pDecorate, entity);
		}
		else{
			fprintf(fp, ".\n");
		}
	}
	bool writeMapFile(cParameter * pParameter){
		if(fpMapfile.fp == NULL){
			return true;
		}
		int i, bc, bc2;
		fprintf(fpMapfile.fp, "#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer\tDescription\n");
		string sampleId, barcode, fw_primer, rv_primer;
		for(i=0; i<cMatrix::iIdxCnt; i++){
			bc = cMatrix::rowBc[i];
			bc2 = cMatrix::colBc[i];
			sampleId = pParameter->rowNames[bc] + pParameter->colNames[bc2];
			barcode = cMatrix::fw_barcodes[bc] + cMatrix::rv_barcodes[bc2];
			fw_primer = cMatrix::fw_primers[bc];
			rv_primer = cMatrix::rv_primers[bc2];
			fprintf(fpMapfile.fp, "%s\t%s\t%s\t%s\tNA\n", sampleId.c_str(), barcode.c_str(), fw_primer.c_str(), rv_primer.c_str());
		}
		return true;
	}
};

typedef enum{
	TASK_READ = 0,
	TASK_WRITE = 1,
	TASK_END = 2
}TASK_TYPE;

typedef struct{
	TASK_TYPE type;
	int nItemCnt;
	int nBlockSize;
	int64 startId;
}TASK;

class cTaskManager
{
private:
	deque<TASK> queue;
	pthread_mutex_t mutex;
	pthread_mutex_t mutex_cnt;
	pthread_mutex_t mutex_item;

	int nItemCnt;
	int nBlockSize;
	int nBufferSize;
	bool bFinished;

	int64 nextId;
public:
	bool bSingleBlock;

	cTaskManager(){
		pthread_mutex_init(&mutex, NULL);
		pthread_mutex_init(&mutex_cnt, NULL);
		pthread_mutex_init(&mutex_item, NULL);
	}
	~cTaskManager(){
		pthread_mutex_destroy(&mutex_item);
		pthread_mutex_destroy(&mutex_cnt);
		pthread_mutex_destroy(&mutex);
	}
	void initialize(int nSize, int nBlockSize, int64 id = 0L){
		TASK task;

		task.type = TASK_READ;
		task.nItemCnt = 0;
		task.nBlockSize = nBlockSize;
		task.startId = id;

		nextId = id;

		nBufferSize = nSize;
		this->nBlockSize = nBlockSize;
		bSingleBlock = (nSize == nBlockSize);
		nItemCnt = 0;

		bFinished = false;
		queue.clear();
		queue.push_back(task);
	}
	void finish(){
		bFinished = true;
	}
	bool IsFinished(){
		return bFinished;
	}
	bool getTask(TASK & task){
		bool bRet;

		pthread_mutex_lock(&mutex);
		if(queue.empty()){
			bRet = false;
		}
		else{
			bRet = true;
			task = queue.front();
			queue.pop_front();
		}
		pthread_mutex_unlock(&mutex);

		return bRet;
	}
	void addTask(TASK & task){
		pthread_mutex_lock(&mutex);
		queue.push_back(task);
		pthread_mutex_unlock(&mutex);
	}
	void insertTask(TASK & task){
		pthread_mutex_lock(&mutex);
		queue.push_front(task);
		pthread_mutex_unlock(&mutex);
	}
	bool increaseCnt(){
		bool bFull;
		pthread_mutex_lock(&mutex_cnt);
		if(nItemCnt + nBlockSize <= nBufferSize){
			bFull = false;
			nItemCnt += nBlockSize;
		}
		else{
			bFull = true;
		}
		pthread_mutex_unlock(&mutex_cnt);

		return !bFull;
	}
	void decreaseCnt(){
		pthread_mutex_lock(&mutex_cnt);
		nItemCnt -= nBlockSize;
		pthread_mutex_unlock(&mutex_cnt);
	}
	int getItemCnt(int id, RECORD *pRecord){
		int nCnt;

		pthread_mutex_lock(&mutex_item);
		nCnt = pRecord->nCnt;
		if(nCnt <= 0){
			nextId = id;
		}
		pthread_mutex_unlock(&mutex_item);

		return nCnt;
	}
	bool setItemCnt(int id, RECORD *pRecord, int nCnt){
		bool bDependent;

		pthread_mutex_lock(&mutex_item);
		if(id == nextId)
			bDependent = false;
		else{
			bDependent = true;
			pRecord->nCnt = nCnt;
		}
		pthread_mutex_unlock(&mutex_item);

		return bDependent;
	}
};

class cData{
public:
	int tid;
	cStats * pStats;
	cTaskManager * pTaskMan;
	RECORD * pBuffer;
	int size;
};

typedef struct tag_mtaux_t{
	int n_threads;
	pthread_t *tid;
	cData *w;
} mtaux_t;

class cWork {
friend class cData;
	cTaskManager taskManager;
	RECORD * pBuffer;
	int size;
	bool bPaired;
	cFQ fq;
	cFQ fq2;
	//int minLen;
	//bool bFivePrimeEnd;

	mtaux_t *mt; // for multi-threading

private:
	inline bool fldEqual(char *a, char *b)
	{
		while (*a && *b){
			if(*a != *b) return false;
			if( (*a == ' ') || (*a == '/') ) return true;
			++a;
			++b;
		}
		return true;
	}

public:
	pthread_attr_t attr;

public:
	cWork(){
		mt = NULL;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
		pBuffer = NULL;
		size = 0;
	}
	~cWork(){
		if(mt != NULL){
			if(mt->w != NULL)
				delete [] mt->w;
			if(mt->tid != NULL)
				free(mt->tid);
			free(mt);
			mt = NULL;
		}
		DestroyBuffer();
		pthread_attr_destroy(&attr);
	}
	bool InitBuffer(int nBuffSize, bool bPaired=false){
		pBuffer = (RECORD *)calloc(nBuffSize * (1 + bPaired), sizeof(RECORD));
		size = nBuffSize;
		this->bPaired = bPaired;
		return (pBuffer != NULL);
	}
	void DestroyBuffer(){
		if(pBuffer == NULL) return;
		RECORD *pBuf;
		int i;
		for(pBuf=pBuffer,i=size * (1+bPaired); i>0; i--,pBuf++){
			if(pBuf->id.s != NULL)
				free(pBuf->id.s);
			if(pBuf->seq.s != NULL)
				free(pBuf->seq.s);
			if(pBuf->com.s != NULL)
				free(pBuf->com.s);
			if(pBuf->qual.s != NULL)
				free(pBuf->qual.s);
		}
		free(pBuffer);
		pBuffer = NULL;
		size = 0;
	}
	bool Init(cParameter * pParameter, cStats * pStats, int64 total_file_length, FILE * fp, FILE * fp2=NULL){
		mt = (mtaux_t *)calloc(1, sizeof(mtaux_t));
		if(mt == NULL)
			return false;
		mt->n_threads = (pParameter->nThreads <= 1) ? 1: pParameter->nThreads;
		mt->tid = (pthread_t *)calloc(mt->n_threads, sizeof(pthread_t));
		mt->w = new cData[mt->n_threads];
		if( (mt->tid == NULL) || (mt->w == NULL) )
			return false;
		bool bPaired = (fp2 != NULL);
		fq.associateFile(fp);
		if(bPaired)	fq2.associateFile(fp2);
		int nBasicSize = (total_file_length > 8 * 100 * 1024L * 1024L) ? 10 : ((total_file_length / 100 / 1024 / 1024) + 2);
		int nBlockSize = nBasicSize * mt->n_threads;
		int nSize = (mt->n_threads * 2 - 1) * nBlockSize; // 1, 3, 5, 7, ..., 31
		if(!InitBuffer(nSize, bPaired))
			return false;

		pStats->InitGlobalAttributes(pParameter, total_file_length, bPaired, &fq, &fq2);

		int i;
		for(i=0; i<mt->n_threads; i++){
			mt->w[i].tid = i;
			mt->w[i].pStats = pStats;
			mt->w[i].pTaskMan = &taskManager;
			mt->w[i].pBuffer = pBuffer;
			mt->w[i].size = size;
		}
		cMatrix::InitParameters(pParameter->trimMode, pParameter->epsilon, pParameter->delta, pParameter->baseQual, pParameter->bShareAdapter);
		cMatrix::iMinOverlap = pParameter->minK;
		vector<string> *pAdapters;
		TRIM_MODE trimMode = ((pParameter->trimMode & TRIM_ANY) == TRIM_DEFAULT) ? TRIM_TAIL : TRIM_MODE(pParameter->trimMode & TRIM_ANY);
		pAdapters = &pParameter->adapters;
		for(i=0; i<int(pAdapters->size()); i++){
			cMatrix::AddAdapter(cMatrix::firstAdapters, (char *)(*pAdapters)[i].c_str(), (*pAdapters)[i].length(), trimMode);
		}
		if(!pParameter->bShareAdapter){
			pAdapters = &pParameter->adapters2;
			for(i=0; i<int(pAdapters->size()); i++){
				cMatrix::AddAdapter(cMatrix::secondAdapters,(char *)(*pAdapters)[i].c_str(), (*pAdapters)[i].length(), trimMode);
			}
		}
		cMatrix::CalculateIndices(pParameter->bMatrix, pParameter->rowNames.size(), pParameter->colNames.size());
		if(bPaired){
			if( (pParameter->trimMode & TRIM_MP) != 0 ){
				pAdapters = &pParameter->juncAdapters;
				for(i=0; i<int(pAdapters->size()); i++){
					cMatrix::AddAdapter(cMatrix::junctionAdapters,(char *)(*pAdapters)[i].c_str(), (*pAdapters)[i].length(), TRIM_ANY);
				}
				cMatrix::CalculateJunctionLengths();
			}
			if( (pParameter->trimMode & TRIM_AP) != 0 ){
				if(pStats->fpMapfile.fp != NULL){
					cMatrix::InitBarcodes(cMatrix::firstAdapters, pParameter->iCutF, (pParameter->bShareAdapter ? cMatrix::firstAdapters : cMatrix::secondAdapters), pParameter->iCutR);
				}
			}
		}
		taskManager.initialize(nSize, nBlockSize);
		return true;
	}
	mtaux_t * getMultiThreadingPointer(){
		return mt;
	}
};

void * mt_worker(void * data)
{
	cData * pData = (cData *)data;
	cTaskManager *pTaskMan = pData->pTaskMan;
	cStats * pStats = pData->pStats;
	int64 file_length = pStats->total_file_length;
	cFQ * pfq = pStats->pfq;
	FILE *fpOut = pStats->fpOut;
	FILE *fpMask = pStats->fpMask;
	FILE *fpExcl = pStats->fpExcl;
	int minAverageQual = pStats->minAverageQual;
	int minEndQual = pStats->minEndQual;
	int minLen = pStats->getMinLen();
	int maxLen = pStats->getMaxLen();
	bool bFivePrimeEnd = pStats->bFivePrimeEnd;
	bool bBarcode = pStats->bBarcode;
	bool bCutTail = pStats->bCutTail;
	
	RECORD * pBuffer, *pRecord;
	TASK task;
	int size, rc, nItemCnt, nCnt;
	int64 startId;

	pBuffer = pData->pBuffer;
	size = pData->size;
	rc = 0;

	int64 cur_pos;
	double cur_ratio;
	int pos;

	while(true){
		while(!pTaskMan->getTask(task)){
			if(pTaskMan->IsFinished()){
				task.type = TASK_END;
				break;
			}
			usleep(1);
		}
		if(task.type == TASK_END){
			break;
		}
		startId = task.startId;
		if(task.type == TASK_READ){
			if(!pTaskMan->increaseCnt()){ // reach the buffer size
				pTaskMan->addTask(task); // perform reading later
				usleep(1);
				continue;
			}
			// read records from input file to buffer
			for(pRecord=&pBuffer[startId % size], nItemCnt=0; nItemCnt<task.nBlockSize; nItemCnt++, pRecord++){
				rc = pfq->readRecord(pRecord);
				if(rc < 0){
					break;
				}
			}
			if(!pStats->bQuiet){
				cur_pos = pfq->tell();
				if(cur_pos >= pfq->next_pos){
					cur_ratio = int64(cur_pos * 10000 / file_length) / 10000.0;
					pStats->progress(cur_ratio, 50);
					pfq->next_pos = int64(((cur_ratio * 10000 + 1) * file_length + 9999)/10000);
				}
			}
			if(rc < 0){ // error or end of file
				pTaskMan->finish();
				if(rc < -1) continue; // error
				if(nItemCnt == 0) continue; // no record read
			}
			task.startId += task.nBlockSize;
			pTaskMan->addTask(task); // save next task for parallelism

			// process the records
			for(pRecord=&pBuffer[startId % size], nCnt=0; nCnt < nItemCnt; nCnt++, pRecord++){
				if( pStats->bFilterNs && cMatrix::isBlurry(pRecord->seq.s, pRecord->seq.n)){
					pRecord->tag = TAG_BLURRY;
					continue;
				}
				if( (minAverageQual > 0) && !cMatrix::checkQualities((uchar *)pRecord->qual.s, pRecord->qual.n, minAverageQual) ){
					pRecord->tag = TAG_BADQUAL;
					continue;
				}
				pRecord->tag = TAG_NORMAL;
				pRecord->idx = cMatrix::findAdapter(pRecord->seq.s, pRecord->seq.n, (uchar *)pRecord->qual.s, pRecord->qual.n);
				if(pRecord->idx.pos < 0){
					pRecord->idx.pos = 0;
				}
				if( (minEndQual > 0) && (pRecord->idx.pos > 0) && (pRecord->qual.n > 0) ){ // not found
					pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->idx.pos, pRecord->qual.n), minEndQual);
				}
			}

			pRecord = &pBuffer[startId % size];
			if(!pTaskMan->setItemCnt(startId, pRecord, nItemCnt)){
				task.type = TASK_WRITE;
				task.startId = startId;
				task.nItemCnt = nItemCnt;
				if(pTaskMan->bSingleBlock)
					pTaskMan->insertTask(task);
				else
					pTaskMan->addTask(task);
			}
			continue;
		}
		// task.type == TASK_WRITE
		pRecord = &pBuffer[startId % size];
		nItemCnt = task.nItemCnt;
		do{
			// write to file
			pRecord->nCnt = 0; // reset
			for(nCnt=0; nCnt < nItemCnt; nCnt++, pRecord++){
				if(pRecord->tag == TAG_BLURRY){
					pStats->nBlurry++;
					if(fpExcl != NULL) {
						OutputTaggedRecord(fpExcl, pRecord);
					}
					continue;
				}
				if(pRecord->tag == TAG_BADQUAL){
					pStats->nBad++;
					if(fpExcl != NULL) {
						OutputTaggedRecord(fpExcl, pRecord);
					}
					continue;
				}
				// TAG_NORMAL
				pos = pRecord->idx.pos;
				if(pos < minLen){
					if(pos <= 0) {
						pStats->nEmpty++;
						pRecord->tag = TAG_EMPTY;
					}
					else {
						pStats->nShort++;
						pRecord->tag = TAG_SHORT;
					}
					if(fpExcl != NULL) {
						OutputTaggedRecord(fpExcl, pRecord);
					}
					continue;
				}
				if(pos > maxLen){
					if(!bCutTail){
						pStats->nLong++;
						pRecord->tag = TAG_LONG;
						if(fpExcl != NULL) {
							OutputTaggedRecord(fpExcl, pRecord);
						}
						continue;
					}
					pos = maxLen;
				}
				
				if(bBarcode){
					int bc = pRecord->idx.bc - 1;
					if(bc < 0){
						fpOut = pStats->fpUntrim.fp;
					}
					else{
						fpOut = pStats->fpOuts[bc].fp;
						pStats->incrementBarcode(bc);
					}
				}
				if(bFivePrimeEnd){
					if( (fpMask != NULL) && (pos < pRecord->seq.n) ){
						OutputMaskedRecord(fpMask, pRecord, pRecord->seq.n - pos, pos);
					}
					if(pStats->bFillWithNs){ // for equal-read-length requirement of some applications
						OutputEntireRecordFilledWithNs(fpOut, pRecord, pRecord->seq.n - pos, pos);
					}
					else{
						OutputPartialRecord(fpOut, pRecord, pRecord->seq.n - pos, pos);
					}
				}
				else{
					if( (fpMask != NULL) && (pos < pRecord->seq.n) ){
						OutputMaskedRecord(fpMask, pRecord, 0, pos);
					}
					if(pStats->bFillWithNs){ // for equal-read-length requirement of some applications
						OutputEntireRecordFilledWithNs(fpOut, pRecord, 0, pos);
					}
					else{
						OutputPartialRecord(fpOut, pRecord, 0, pos);
					}
				}
				if(bBarcode){
					if(pRecord->idx.bc == 0)
						pStats->nUntrimAvail++;
					else
						pStats->nTrimAvail++;
				}
				else{
					if(pos < pRecord->seq.n)
						pStats->nTrimAvail++;
					else
						pStats->nUntrimAvail++;
				}
				pStats->incrementCount(size_t(pos));
			}
			pTaskMan->decreaseCnt();
			startId += task.nBlockSize;
			pRecord = &pBuffer[startId % size];
			nItemCnt = pTaskMan->getItemCnt(startId, pRecord);
		}while(nItemCnt > 0);
	}
	return NULL;
}

void * mt_worker_ap(void * data)
{
	cData * pData = (cData *)data;
	cTaskManager *pTaskMan = pData->pTaskMan;
	cStats * pStats = pData->pStats;
	int64 file_length = pStats->total_file_length;
	cFQ * pfq = pStats->pfq;
	FILE *fpOut = pStats->fpOut;
	FILE *fpMask = pStats->fpMask;
	FILE *fpExcl = pStats->fpExcl;
	FILE *fpBarcode = pStats->fpBarcode;
	int minAverageQual = pStats->minAverageQual;
	int minEndQual = pStats->minEndQual;
	int minLen = pStats->getMinLen();
	int maxLen = pStats->getMaxLen();
	bool bBarcode = pStats->bBarcode;
	int iCutF = pStats->iCutF;
	int iCutR = pStats->iCutR;
	bool bCutTail = pStats->bCutTail;
	char barcodeSeq[25];
	char barcodeQua[25];
	if(fpBarcode != NULL){
		memset(barcodeSeq, 0, 25 * sizeof(char));
		memset(barcodeQua, 0, 25 * sizeof(char));
	}
	assert( (iCutF + iCutR) <= 24 );

	RECORD *pBuffer, *pRecord;
	TASK task;
	int size, rc, nItemCnt, nCnt;
	int flag;
	int64 startId;

	pBuffer = pData->pBuffer;
	size = pData->size;
	rc = 0;

	int64 cur_pos;
	double cur_ratio;
	int pos;

	while(true){
		while(!pTaskMan->getTask(task)){
			if(pTaskMan->IsFinished()){
				task.type = TASK_END;
				break;
			}
			usleep(1);
		}
		if(task.type == TASK_END){
			break;
		}
		startId = task.startId;
		if(task.type == TASK_READ){
			if(!pTaskMan->increaseCnt()){ // reach the buffer size
				pTaskMan->addTask(task); // perform reading later
				usleep(1);
				continue;
			}
			// read records from input file to buffer
			for(pRecord=&pBuffer[startId % size], nItemCnt=0; nItemCnt<task.nBlockSize; nItemCnt++, pRecord++){
				rc = pfq->readRecord(pRecord);
				if(rc < 0){
					break;
				}
			}
			if(!pStats->bQuiet){
				cur_pos = pfq->tell();
				if(cur_pos >= pfq->next_pos){
					cur_ratio = int64(cur_pos * 10000 / file_length) / 10000.0;
					pStats->progress(cur_ratio, 50);
					pfq->next_pos = int64(((cur_ratio * 10000 + 1) * file_length + 9999)/10000);
				}
			}
			if(rc < 0){ // error or end of file
				pTaskMan->finish();
				if(rc < -1) continue; // error
				if(nItemCnt == 0) continue; // no record read
			}
			task.startId += task.nBlockSize;
			pTaskMan->addTask(task); // save next task for parallelism

			// process the records
			for(pRecord=&pBuffer[startId % size], nCnt=0; nCnt < nItemCnt; nCnt++, pRecord++){
				if( pStats->bFilterNs && cMatrix::isBlurry(pRecord->seq.s, pRecord->seq.n)){
					pRecord->tag = TAG_BLURRY;
					continue;
				}
				if( (minAverageQual > 0) && !cMatrix::checkQualities((uchar *)pRecord->qual.s, pRecord->qual.n, minAverageQual) ){
					pRecord->tag = TAG_BADQUAL;
					continue;
				}
				pRecord->tag = TAG_NORMAL;
				flag = cMatrix::findAdaptersInARead(pRecord->seq.s, pRecord->seq.n, (uchar *)pRecord->qual.s, pRecord->qual.n, pRecord->idx);
				// TODO: 
				if(flag >= 0){
					pRecord->bExchange = (flag == 1);
					if(flag == 0){
						pRecord->idx.pos = pRecord->seq.n - iCutF;
					}
					else{
						pRecord->idx.pos = pRecord->seq.n - iCutR;
					}
					if(pRecord->idx.pos < 0){
						pRecord->idx.pos = 0;
					}
				}
				if( minEndQual > 0 ){ // TODO: quality trimming from 5' end
					if( (pRecord->idx.pos > 0) && (pRecord->qual.n > 0) ){
						pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->idx.pos, pRecord->qual.n), minEndQual);
					}
				}
			}

			pRecord = &pBuffer[startId % size];
			if(!pTaskMan->setItemCnt(startId, pRecord, nItemCnt)){
				task.type = TASK_WRITE;
				task.startId = startId;
				task.nItemCnt = nItemCnt;
				if(pTaskMan->bSingleBlock)
					pTaskMan->insertTask(task);
				else
					pTaskMan->addTask(task);
			}
			continue;
		}
		// task.type == TASK_WRITE
		pRecord = &pBuffer[startId % size];
		nItemCnt = task.nItemCnt;
		do{
			// write to file
			pRecord->nCnt = 0; // reset
			for(nCnt=0; nCnt<nItemCnt; nCnt++, pRecord++){
				if(pRecord->tag == TAG_BLURRY){
					pStats->nBlurry++;
					if(fpExcl != NULL) {
						OutputTaggedRecord(fpExcl, pRecord);
					}
					continue;
				}
				if(pRecord->tag == TAG_BADQUAL){
					pStats->nBad++;
					if(fpExcl != NULL) {
						OutputTaggedRecord(fpExcl, pRecord);
					}
					continue;
				}
				// TAG_NORMAL
				pos = pRecord->idx.pos;
				if(pos < minLen){
					if(pos <= 0) {
						pStats->nEmpty++;
						pRecord->tag = TAG_EMPTY;
					}
					else{
						pStats->nShort++;
						pRecord->tag = TAG_SHORT;
					}
					if(fpExcl != NULL) {
						OutputTaggedRecord(fpExcl, pRecord);
					}
					continue;
				}
				if(pos > maxLen){
					if(!bCutTail){
						pStats->nLong++;
						pRecord->tag = TAG_LONG;
						if(fpExcl != NULL) {
							OutputTaggedRecord(fpExcl, pRecord);
						}
						continue;
					}
					pos = maxLen;
				}
				if(bBarcode){
					if(pRecord->idx.bc < 0){ // to be checked
						fpOut = pStats->fpUntrim.fp;
					}
					else{
						fpOut = pStats->fpOuts[pRecord->idx.bc].fp;
						pStats->incrementBarcode(pRecord->idx.bc);
					}
					OutputEntireRecord(fpOut, pRecord);
				}
				else{
					int iCut = (pRecord->bExchange ? iCutR : iCutF);
					if(fpMask != NULL){
						if(pos < pRecord->seq.n - iCut)
							OutputMaskedRecord(fpMask, pRecord, iCut, pos);
					}
					if(pStats->bFillWithNs){ // for equal-read-length requirement of some applications
						OutputEntireRecordFilledWithNs(fpOut, pRecord, iCut, pos);
					}
					else{
						OutputPartialRecord(fpOut, pRecord, iCut, pos);
					}
				}
				
				if(pRecord->idx.bc < 0){ // assigned
					pStats->nUntrimAvail++;
				}
				else{
					pStats->nTrimAvail++;
				}
				pStats->incrementCount(size_t(pos));
			}
			pTaskMan->decreaseCnt();
			startId += task.nBlockSize;
			pRecord = &pBuffer[startId % size];
			nItemCnt = pTaskMan->getItemCnt(startId, pRecord);
		}while(nItemCnt > 0);
	}
	return NULL;
}

void * mt_worker2(void * data)
{
	cData * pData = (cData *)data;
	cTaskManager *pTaskMan = pData->pTaskMan;
	cStats * pStats = pData->pStats;
	int64 file_length = pStats->total_file_length;
	cFQ * pfq = pStats->pfq;
	cFQ * pfq2 = pStats->pfq2;
	FILE *fpOut = pStats->fpOut;
	FILE *fpOut2 = pStats->fpOut2;
	FILE *fpMask = pStats->fpMask;
	FILE *fpMask2 = pStats->fpMask2;
	FILE *fpExcl = pStats->fpExcl;
	FILE *fpExcl2 = pStats->fpExcl2;
	int minAverageQual = pStats->minAverageQual;
	int minEndQual = pStats->minEndQual;
	int minLen = pStats->getMinLen();
	int maxLen = pStats->getMaxLen();
	bool bBarcode = pStats->bBarcode;
	bool bCutTail = pStats->bCutTail;

	RECORD *pBuffer, *pRecord, *pRecord2;
	TASK task;
	int size2, rc, rc2, nItemCnt, nCnt;
	int64 startId;

	pBuffer = pData->pBuffer;
	size2 = pData->size * 2;
	rc = rc2 = 0;

	int64 cur_pos;
	double cur_ratio;
	int pos, pos2, mLen;

	int rLen, qLen;
	int rLen2, qLen2;

	while(true){
		while(!pTaskMan->getTask(task)){
			if(pTaskMan->IsFinished()){
				task.type = TASK_END;
				break;
			}
			usleep(1);
		}
		if(task.type == TASK_END){
			break;
		}
		startId = task.startId;
		if(task.type == TASK_READ){
			if(!pTaskMan->increaseCnt()){ // reach the buffer size
				pTaskMan->addTask(task); // perform reading later
				usleep(1);
				continue;
			}
			// read records from input file to buffer
			for(pRecord=&pBuffer[(startId << 1) % size2], nItemCnt=0; nItemCnt<task.nBlockSize; nItemCnt++, pRecord+=2){
				pRecord2 = pRecord+1;
				rc = pfq->readRecord(pRecord);
				rc2 = pfq2->readRecord(pRecord2);
				if( (rc < 0) || (rc2 < 0) ){
					break;
				}
			}
			if(!pStats->bQuiet){
				cur_pos = pfq->tell() + pfq2->tell();
				if(cur_pos >= pfq->next_pos){
					cur_ratio = int64(cur_pos * 10000 / file_length) / 10000.0;
					pStats->progress(cur_ratio, 50);
					pfq->next_pos = int64(((cur_ratio * 10000 + 1) * file_length + 9999)/10000);
				}
			}
			if( (rc < 0) || (rc2 < 0) ){ // error or end of file
				pTaskMan->finish();
				if( (rc < -1) || (rc2 < -1) ) continue; // error
				if(nItemCnt == 0) continue; // no record read
			}
			task.startId += task.nBlockSize;
			pTaskMan->addTask(task); // save next task for parallelism

			// process the records
			for(pRecord=&pBuffer[(startId << 1) % size2], nCnt=0; nCnt < nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if(minAverageQual > 0){
					if( !cMatrix::checkQualities((uchar *)pRecord->qual.s, pRecord->qual.n, minAverageQual) &&
						!cMatrix::checkQualities((uchar *)pRecord2->qual.s, pRecord2->qual.n, minAverageQual) ){
						pRecord->tag = pRecord2->tag = TAG_BADQUAL;
						continue;
					}
				}
				pRecord->tag = TAG_NORMAL;
				rLen = pRecord->seq.n;
				qLen = pRecord->qual.n;
				rLen2 = pRecord2->seq.n;
				qLen2 = pRecord2->qual.n;
				if(cMatrix::findAdapterWithPE(pRecord->seq.s, pRecord2->seq.s, rLen, rLen2,
						 (uchar *)pRecord->qual.s, (uchar *)pRecord2->qual.s, qLen, qLen2,
						 pRecord->idx, pRecord2->idx)){ // trimmed
					pos = pRecord->idx.pos;
					pos2 = pRecord2->idx.pos;
					if( (pos >= minLen) && (pos2 >= minLen) ){
						cMatrix::combinePairSeqs(pRecord->seq.s, pRecord2->seq.s, pos, pos2,
							(uchar *)pRecord->qual.s, (uchar *)pRecord2->qual.s, (int)qLen, (int)qLen2);
						if(pStats->bFilterNs){
							if(cMatrix::isBlurry(pRecord->seq.s, pos) && cMatrix::isBlurry(pRecord2->seq.s, pos2)){
								pRecord->tag = pRecord2->tag = TAG_BLURRY;
							}
						}
					}
				}
				else{
					pos = rLen;
					pos2 = rLen2;
					if(pStats->bFilterNs){
						if( cMatrix::isBlurry(pRecord->seq.s, rLen) && cMatrix::isBlurry(pRecord2->seq.s, rLen2) ){
							pRecord->tag = pRecord2->tag = TAG_BLURRY;
						}
					}
				}
				if( (minEndQual > 0) && (pRecord->tag == TAG_NORMAL) ){ // trimmed by quality
					if( (pRecord->qual.n > 0) && (pos > 0) )
						pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pos, pRecord->qual.n), minEndQual);
					if( (pRecord2->qual.n > 0) && (pos2 > 0) )
						pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, min(pos2, pRecord2->qual.n), minEndQual);
				}
			}

			pRecord = &pBuffer[(startId << 1) % size2];
			if(!pTaskMan->setItemCnt(startId, pRecord, nItemCnt)){
				task.type = TASK_WRITE;
				task.startId = startId;
				task.nItemCnt = nItemCnt;
				if(pTaskMan->bSingleBlock)
					pTaskMan->insertTask(task);
				else
					pTaskMan->addTask(task);
			}
			continue;
		}
		// task.type == TASK_WRITE
		pRecord = &pBuffer[(startId << 1) % size2];
		nItemCnt = task.nItemCnt;
		do{
			// write to file
			pRecord->nCnt = 0; // reset
			for(nCnt=0; nCnt<nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if(pRecord->tag == TAG_BLURRY){
					pStats->nBlurry++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pRecord->tag == TAG_BADQUAL){
					pStats->nBad++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				// TAG_NORMAL
				if(bBarcode){
					pos = pRecord->seq.n;
					pos2 = pRecord2->seq.n;
				}
				else{
					pos = pRecord->idx.pos;
					pos2 = pRecord2->idx.pos;
				}
				if( (pos < minLen) || (pos2 < minLen) ){
					if( (pos <= 0) || (pos2 <= 0) ) {
						pStats->nEmpty++;
						pRecord->tag = pRecord2->tag = TAG_EMPTY;
					}
					else {
						pStats->nShort++;
						pRecord->tag = pRecord2->tag = TAG_SHORT;
					}
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if( (pos > maxLen) || (pos2 > maxLen) ){
					if(!bCutTail){
						pStats->nLong++;
						pRecord->tag = pRecord2->tag = TAG_LONG;
						if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
							OutputTaggedRecord(fpExcl, pRecord);
							OutputTaggedRecord(fpExcl2, pRecord2);
						}
						continue;
					}
					if(pos > maxLen) pos = maxLen;
					if(pos2 > maxLen) pos2 = maxLen;
				}

				if(bBarcode){
					if(pRecord->idx.bc < 0){
						fpOut = pStats->fpUntrim.fp;
						fpOut2 = pStats->fpUntrim2.fp;
					}
					else{
						fpOut = pStats->fpOuts[pRecord->idx.bc].fp;
						fpOut2 = pStats->fpOuts2[pRecord->idx.bc].fp;
						pStats->incrementBarcode(pRecord->idx.bc);
					}
				}
				if( (fpMask != NULL) && (fpMask2 != NULL) ){
					if( (pos < pRecord->seq.n) || (pos2 < pRecord2->seq.n) ){
						OutputMaskedRecord(fpMask, pRecord, 0, pos);
						OutputMaskedRecord(fpMask2, pRecord2, 0, pos2);
					}
				}
				if(pStats->bFillWithNs){ // for equal-read-length requirement of some applications
					OutputEntireRecordFilledWithNs(fpOut, pRecord, 0, pos);
					OutputEntireRecordFilledWithNs(fpOut2, pRecord2, 0, pos2);
				}
				else{
					OutputPartialRecord(fpOut, pRecord, 0, pos);
					OutputPartialRecord(fpOut2, pRecord2, 0, pos2);
				}
				if(bBarcode){
					if(pRecord->idx.bc < 0){ // assigned
						pStats->nUntrimAvail++;
					}
					else{
						pStats->nTrimAvail++;
					}
				}
				else{
					if(pRecord->idx.pos < pRecord->seq.n || pRecord2->idx.pos < pRecord2->seq.n) // trimmed
						pStats->nTrimAvail++;
					else
						pStats->nUntrimAvail++;
				}
				mLen = (pos + pos2) / 2;
				pStats->incrementCount(size_t(mLen));
			}
			pTaskMan->decreaseCnt();
			startId += task.nBlockSize;
			pRecord = &pBuffer[(startId << 1) % size2];
			nItemCnt = pTaskMan->getItemCnt(startId, pRecord);
		}while(nItemCnt > 0);
	}
	return NULL;
}

void * mt_worker2_sep(void * data)
{
	cData * pData = (cData *)data;
	cTaskManager *pTaskMan = pData->pTaskMan;
	cStats * pStats = pData->pStats;
	int64 file_length = pStats->total_file_length;
	cFQ * pfq = pStats->pfq;
	cFQ * pfq2 = pStats->pfq2;
	FILE *fpOut = pStats->fpOut;
	FILE *fpOut2 = pStats->fpOut2;
	FILE *fpMask = pStats->fpMask;
	FILE *fpMask2 = pStats->fpMask2;
	FILE *fpExcl = pStats->fpExcl;
	FILE *fpExcl2 = pStats->fpExcl2;
	int minAverageQual = pStats->minAverageQual;
	int minEndQual = pStats->minEndQual;
	int minLen = pStats->getMinLen();
	int maxLen = pStats->getMaxLen();
	bool bFivePrimeEnd = pStats->bFivePrimeEnd;
	bool bBarcode = pStats->bBarcode;
	bool bCutTail = pStats->bCutTail;

	RECORD *pBuffer, *pRecord, *pRecord2;
	TASK task;
	int size2, rc, rc2, nItemCnt, nCnt;
	int64 startId;

	pBuffer = pData->pBuffer;
	size2 = pData->size * 2;
	rc = rc2 = 0;

	int64 cur_pos;
	double cur_ratio;
	int pos, pos2, mLen;

	while(true){
		while(!pTaskMan->getTask(task)){
			if(pTaskMan->IsFinished()){
				task.type = TASK_END;
				break;
			}
			usleep(1);
		}
		if(task.type == TASK_END){
			break;
		}
		startId = task.startId;
		if(task.type == TASK_READ){
			if(!pTaskMan->increaseCnt()){ // reach the buffer size
				pTaskMan->addTask(task); // perform reading later
				usleep(1);
				continue;
			}
			// read records from input file to buffer
			for(pRecord=&pBuffer[(startId << 1) % size2], nItemCnt=0; nItemCnt<task.nBlockSize; nItemCnt++, pRecord+=2){
				pRecord2 = pRecord+1;
				rc = pfq->readRecord(pRecord);
				rc2 = pfq2->readRecord(pRecord2);
				if( (rc < 0) || (rc2 < 0) ){
					break;
				}
			}
			if(!pStats->bQuiet){
				cur_pos = pfq->tell() + pfq2->tell();
				if(cur_pos >= pfq->next_pos){
					cur_ratio = int64(cur_pos * 10000 / file_length) / 10000.0;
					pStats->progress(cur_ratio, 50);
					pfq->next_pos = int64(((cur_ratio * 10000 + 1) * file_length + 9999)/10000);
				}
			}
			if( (rc < 0) || (rc2 < 0) ){ // error or end of file
				pTaskMan->finish();
				if( (rc < -1) || (rc2 < -1) ) continue; // error
				if(nItemCnt == 0) continue; // no record read
			}
			task.startId += task.nBlockSize;
			pTaskMan->addTask(task); // save next task for parallelism

			// process the records
			for(pRecord=&pBuffer[(startId << 1) % size2], nCnt=0; nCnt < nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if( pStats->bFilterNs &&
					(cMatrix::isBlurry(pRecord->seq.s, pRecord->seq.n) && 
					cMatrix::isBlurry(pRecord2->seq.s, pRecord2->seq.n)) ){
					pRecord->tag = pRecord2->tag = TAG_BLURRY;
					continue;
				}
				if(minAverageQual > 0){
					if( !cMatrix::checkQualities((uchar *)pRecord->qual.s, pRecord->qual.n, minAverageQual) &&
						!cMatrix::checkQualities((uchar *)pRecord2->qual.s, pRecord2->qual.n, minAverageQual) ){
						pRecord->tag = pRecord2->tag = TAG_BADQUAL;
						continue;
					}
				}
				pRecord->tag = TAG_NORMAL;
				pRecord->idx = cMatrix::findAdapter(pRecord->seq.s, pRecord->seq.n, (uchar *)pRecord->qual.s, pRecord->qual.n);
				pRecord2->idx = cMatrix::findAdapter2(pRecord2->seq.s, pRecord2->seq.n, (uchar *)pRecord2->qual.s, pRecord2->qual.n);
				if(pRecord->idx.pos < 0){
					pRecord->idx.pos = 0;
				}
				if(pRecord2->idx.pos < 0){
					pRecord2->idx.pos = 0;
				}
				if( minEndQual > 0 ){
					if( (pRecord->idx.pos > 0) && (pRecord->qual.n > 0) ){
						pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->idx.pos, pRecord->qual.n), minEndQual);
					}
					if( (pRecord2->idx.pos > 0) && (pRecord2->qual.n > 0) ){
						pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, min(pRecord2->idx.pos, pRecord2->qual.n), minEndQual);
					}
				}
			}

			pRecord = &pBuffer[(startId << 1) % size2];
			if(!pTaskMan->setItemCnt(startId, pRecord, nItemCnt)){
				task.type = TASK_WRITE;
				task.startId = startId;
				task.nItemCnt = nItemCnt;
				if(pTaskMan->bSingleBlock)
					pTaskMan->insertTask(task);
				else
					pTaskMan->addTask(task);
			}
			continue;
		}
		// task.type == TASK_WRITE
		pRecord = &pBuffer[(startId << 1) % size2];
		nItemCnt = task.nItemCnt;
		do{
			// write to file
			pRecord->nCnt = 0; // reset
			for(nCnt=0; nCnt<nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if(pRecord->tag == TAG_BLURRY){
					pStats->nBlurry++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pRecord->tag == TAG_BADQUAL){
					pStats->nBad++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				// TAG_NORMAL
				if(bBarcode){
					pos = pRecord->seq.n;
					pos2 = pRecord2->seq.n;
				}
				else{
					pos = pRecord->idx.pos;
					pos2 = pRecord2->idx.pos;
				}
				if( (pos < minLen) || (pos2 < minLen) ){
					if( (pos <= 0) || (pos2 <= 0) ) {
						pStats->nEmpty++;
						pRecord->tag = pRecord2->tag = TAG_EMPTY;
					}
					else {
						pStats->nShort++;
						pRecord->tag = pRecord2->tag = TAG_SHORT;
					}
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if( (pos > maxLen) || (pos2 > maxLen) ){
					if(!bCutTail){
						pStats->nLong++;
						pRecord->tag = pRecord2->tag = TAG_LONG;
						if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
							OutputTaggedRecord(fpExcl, pRecord);
							OutputTaggedRecord(fpExcl2, pRecord2);
						}
						continue;
					}
					if(pos > maxLen) pos = maxLen;
					if(pos2 > maxLen) pos2 = maxLen;
				}
				if(bBarcode){
					int bc = cMatrix::indices[pRecord->idx.bc][pRecord2->idx.bc];
					if( bc < 0){
						fpOut = pStats->fpUntrim.fp;
						fpOut2 = pStats->fpUntrim2.fp;
					}
					else{
						fpOut = pStats->fpOuts[bc].fp;
						fpOut2 = pStats->fpOuts2[bc].fp;
						pStats->incrementBarcode(bc);
					}
				}
				if(bFivePrimeEnd){
					if( (fpMask != NULL) && (fpMask2 != NULL) ){
						if( (pos < pRecord->seq.n) || (pos2 < pRecord2->seq.n) ){
							OutputMaskedRecord(fpMask, pRecord, pRecord->seq.n - pos, pos);
							OutputMaskedRecord(fpMask2, pRecord2, pRecord2->seq.n - pos2, pos2);
						}
					}
					if(pStats->bFillWithNs){ // for equal-read-length requirement of some applications
						OutputEntireRecordFilledWithNs(fpOut, pRecord, pRecord->seq.n - pos, pos);
						OutputEntireRecordFilledWithNs(fpOut2, pRecord2, pRecord2->seq.n - pos2, pos2);
					}
					else{
						OutputPartialRecord(fpOut, pRecord, pRecord->seq.n - pos, pos);
						OutputPartialRecord(fpOut2, pRecord2, pRecord2->seq.n - pos2, pos2);
					}
				}
				else{
					if( (fpMask != NULL) && (fpMask2 != NULL) ){
						if( (pos < pRecord->seq.n) || (pos2 < pRecord2->seq.n) ){
							OutputMaskedRecord(fpMask, pRecord, 0, pos);
							OutputMaskedRecord(fpMask2, pRecord2, 0, pos2);
						}
					}
					if(pStats->bFillWithNs){ // for equal-read-length requirement of some applications
						OutputEntireRecordFilledWithNs(fpOut, pRecord, 0, pos);
						OutputEntireRecordFilledWithNs(fpOut2, pRecord2, 0, pos2);
					}
					else{
						OutputPartialRecord(fpOut, pRecord, 0, pos);
						OutputPartialRecord(fpOut2, pRecord2, 0, pos2);
					}
				}
				if(bBarcode){
					if(cMatrix::indices[pRecord->idx.bc][pRecord2->idx.bc] < 0){ // assigned
						pStats->nUntrimAvail++;
					}
					else{
						pStats->nTrimAvail++;
					}
				}
				else{
					if(pRecord->idx.pos < pRecord->seq.n || pRecord2->idx.pos < pRecord2->seq.n) // trimmed
						pStats->nTrimAvail++;
					else
						pStats->nUntrimAvail++;
				}
				mLen = (pos + pos2) / 2;
				pStats->incrementCount(size_t(mLen));
			}
			pTaskMan->decreaseCnt();
			startId += task.nBlockSize;
			pRecord = &pBuffer[(startId << 1) % size2];
			nItemCnt = pTaskMan->getItemCnt(startId, pRecord);
		}while(nItemCnt > 0);
	}
	return NULL;
}

void * mt_worker2_ap(void * data)
{
	cData * pData = (cData *)data;
	cTaskManager *pTaskMan = pData->pTaskMan;
	cStats * pStats = pData->pStats;
	int64 file_length = pStats->total_file_length;
	cFQ * pfq = pStats->pfq;
	cFQ * pfq2 = pStats->pfq2;
	FILE *fpOut = pStats->fpOut;
	FILE *fpOut2 = pStats->fpOut2;
	FILE *fpMask = pStats->fpMask;
	FILE *fpMask2 = pStats->fpMask2;
	FILE *fpExcl = pStats->fpExcl;
	FILE *fpExcl2 = pStats->fpExcl2;
	FILE *fpBarcode = pStats->fpBarcode;
	int minAverageQual = pStats->minAverageQual;
	int minEndQual = pStats->minEndQual;
	int minLen = pStats->getMinLen();
	int maxLen = pStats->getMaxLen();
	bool bBarcode = pStats->bBarcode;
	int iCutF = pStats->iCutF;
	int iCutR = pStats->iCutR;
	char barcodeSeq[25];
	char barcodeQua[25];
	bool bCutTail = pStats->bCutTail;
	if(fpBarcode != NULL){
		memset(barcodeSeq, 0, 25 * sizeof(char));
		memset(barcodeQua, 0, 25 * sizeof(char));
	}
	assert( (iCutF + iCutR) <= 24 );

	RECORD *pBuffer, *pRecord, *pRecord2;
	TASK task;
	int size2, rc, rc2, nItemCnt, nCnt;
	int flag;
	int64 startId;

	pBuffer = pData->pBuffer;
	size2 = pData->size * 2;
	rc = rc2 = 0;

	int64 cur_pos;
	double cur_ratio;
	int pos, pos2, mLen;

	while(true){
		while(!pTaskMan->getTask(task)){
			if(pTaskMan->IsFinished()){
				task.type = TASK_END;
				break;
			}
			usleep(1);
		}
		if(task.type == TASK_END){
			break;
		}
		startId = task.startId;
		if(task.type == TASK_READ){
			if(!pTaskMan->increaseCnt()){ // reach the buffer size
				pTaskMan->addTask(task); // perform reading later
				usleep(1);
				continue;
			}
			// read records from input file to buffer
			for(pRecord=&pBuffer[(startId << 1) % size2], nItemCnt=0; nItemCnt<task.nBlockSize; nItemCnt++, pRecord+=2){
				pRecord2 = pRecord+1;
				rc = pfq->readRecord(pRecord);
				rc2 = pfq2->readRecord(pRecord2);
				if( (rc < 0) || (rc2 < 0) ){
					break;
				}
			}
			if(!pStats->bQuiet){
				cur_pos = pfq->tell() + pfq2->tell();
				if(cur_pos >= pfq->next_pos){
					cur_ratio = int64(cur_pos * 10000 / file_length) / 10000.0;
					pStats->progress(cur_ratio, 50);
					pfq->next_pos = int64(((cur_ratio * 10000 + 1) * file_length + 9999)/10000);
				}
			}
			if( (rc < 0) || (rc2 < 0) ){ // error or end of file
				pTaskMan->finish();
				if( (rc < -1) || (rc2 < -1) ) continue; // error
				if(nItemCnt == 0) continue; // no record read
			}
			task.startId += task.nBlockSize;
			pTaskMan->addTask(task); // save next task for parallelism

			// process the records
			for(pRecord=&pBuffer[(startId << 1) % size2], nCnt=0; nCnt < nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if( pStats->bFilterNs &&
					(cMatrix::isBlurry(pRecord->seq.s, pRecord->seq.n) && 
					cMatrix::isBlurry(pRecord2->seq.s, pRecord2->seq.n)) ){
					pRecord->tag = pRecord2->tag = TAG_BLURRY;
					continue;
				}
				if(minAverageQual > 0){
					if( !cMatrix::checkQualities((uchar *)pRecord->qual.s, pRecord->qual.n, minAverageQual) &&
						!cMatrix::checkQualities((uchar *)pRecord2->qual.s, pRecord2->qual.n, minAverageQual) ){
						pRecord->tag = pRecord2->tag = TAG_BADQUAL;
						continue;
					}
				}
				pRecord->tag = TAG_NORMAL;
				flag = cMatrix::findAdaptersBidirectionally(pRecord->seq.s, pRecord->seq.n, (uchar *)pRecord->qual.s, pRecord->qual.n,
									pRecord2->seq.s, pRecord2->seq.n, (uchar *)pRecord2->qual.s, pRecord2->qual.n, pRecord->idx, pRecord2->idx);
				if(flag >= 0){
					pRecord->bExchange = (flag == 1);
					if(flag == 0){
						pRecord->idx.pos = pRecord->seq.n - iCutF;
						pRecord2->idx.pos = pRecord2->seq.n - iCutR;
					}
					else{
						pRecord->idx.pos = pRecord->seq.n - iCutR;
						pRecord2->idx.pos = pRecord2->seq.n - iCutF;
					}
					if(pRecord->idx.pos < 0){
						pRecord->idx.pos = 0;
					}
					if(pRecord2->idx.pos < 0){
						pRecord2->idx.pos = 0;
					}
				}
				if( minEndQual > 0 ){ // TODO: quality trimming from 5' end
					if( (pRecord->idx.pos > 0) && (pRecord->qual.n > 0) ){
						pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->idx.pos, pRecord->qual.n), minEndQual);
					}
					if( (pRecord2->idx.pos > 0) && (pRecord2->qual.n > 0) ){
						pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, min(pRecord2->idx.pos, pRecord2->qual.n), minEndQual);
					}
				}
			}

			pRecord = &pBuffer[(startId << 1) % size2];
			if(!pTaskMan->setItemCnt(startId, pRecord, nItemCnt)){
				task.type = TASK_WRITE;
				task.startId = startId;
				task.nItemCnt = nItemCnt;
				if(pTaskMan->bSingleBlock)
					pTaskMan->insertTask(task);
				else
					pTaskMan->addTask(task);
			}
			continue;
		}
		// task.type == TASK_WRITE
		pRecord = &pBuffer[(startId << 1) % size2];
		nItemCnt = task.nItemCnt;
		do{
			// write to file
			pRecord->nCnt = 0; // reset
			for(nCnt=0; nCnt<nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if(pRecord->tag == TAG_BLURRY){
					pStats->nBlurry++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pRecord->tag == TAG_BADQUAL){
					pStats->nBad++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				// TAG_NORMAL
				pos = pRecord->idx.pos;
				pos2 = pRecord2->idx.pos;
				if( (pos < minLen) || (pos2 < minLen) ){
					if( (pos <= 0) || (pos2 <= 0) ) {
						pStats->nEmpty++;
						pRecord->tag = pRecord2->tag = TAG_EMPTY;
					}
					else {
						pStats->nShort++;
						pRecord->tag = pRecord2->tag = TAG_SHORT;
					}
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if( (pos > maxLen) || (pos2 > maxLen) ){
					if(!bCutTail){
						pStats->nLong++;
						pRecord->tag = pRecord2->tag = TAG_LONG;
						if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
							OutputTaggedRecord(fpExcl, pRecord);
							OutputTaggedRecord(fpExcl2, pRecord2);
						}
						continue;
					}
					if(pos > maxLen) pos = maxLen;
					if(pos2 > maxLen) pos2 = maxLen;
				}
				
				if(bBarcode){
					if(pRecord->idx.bc < 0){
						fpOut = pStats->fpUntrim.fp;
						fpOut2 = pStats->fpUntrim2.fp;
					}
					else{
						if(pRecord->bExchange){
							fpOut2 = pStats->fpOuts[pRecord->idx.bc].fp;
							fpOut = pStats->fpOuts2[pRecord->idx.bc].fp;
						}
						else{
							fpOut = pStats->fpOuts[pRecord->idx.bc].fp;
							fpOut2 = pStats->fpOuts2[pRecord->idx.bc].fp;
						}
						pStats->incrementBarcode(pRecord->idx.bc);
					}
					OutputEntireRecord(fpOut, pRecord);
					OutputEntireRecord(fpOut2, pRecord2);
				}
				else{
					int iCut = (pRecord->bExchange ? iCutR : iCutF);
					int iCut2 = (pRecord->bExchange ? iCutF : iCutR);
					if( (fpMask != NULL) && (fpMask2 != NULL) ){
						if( (pos < pRecord->seq.n - iCut) || (pos2 < pRecord2->seq.n - iCut2) ){
							OutputMaskedRecord(fpMask, pRecord, iCut, pos);
							OutputMaskedRecord(fpMask2, pRecord2, iCut2, pos2);
						}
					}
					if(pStats->bFillWithNs){
						OutputEntireRecordFilledWithNs(fpOut, pRecord, iCut, pos);
						OutputEntireRecordFilledWithNs(fpOut2, pRecord2, iCut2, pos2);
					}
					else{
						OutputPartialRecord(fpOut, pRecord, iCut, pos);
						OutputPartialRecord(fpOut2, pRecord2, iCut2, pos2);
					}
				}
				if( (fpBarcode != NULL) && (pRecord->idx.bc >= 0) ){
					if( (pRecord->com.n > 0) && (pRecord2->com.n > 0) ){ // fastq
						if(cMatrix::PrepareBarcode(barcodeSeq, pRecord->idx.bc, pRecord->seq.s, iCutF, pRecord2->seq.s, iCutR, barcodeQua, pRecord->qual.s, pRecord2->qual.s)){
							fprintf(fpBarcode, "@%s%s\n+\n%s\n", pRecord->id.s, barcodeSeq, barcodeQua);
						}
						else{
							fprintf(fpBarcode, "@%s%.*s%.*s\n+\n%.*s%.*s\n", pRecord->id.s, iCutF, pRecord->seq.s, iCutR, pRecord2->seq.s,
										iCutF, pRecord->qual.s, iCutR, pRecord2->qual.s);
						}
					}
					else{ // fasta
						if(cMatrix::PrepareBarcode(barcodeSeq, pRecord->idx.bc, pRecord->seq.s, iCutF, pRecord2->seq.s, iCutR)){
							fprintf(fpBarcode, ">%s%s\n", pRecord->id.s, barcodeSeq);
						}
						else{
							fprintf(fpBarcode, ">%s%.*s%.*s\n", pRecord->id.s, iCutF, pRecord->seq.s, iCutR, pRecord2->seq.s);
						}
					}
				}
				if(pRecord->idx.bc < 0){ // assigned
					pStats->nUntrimAvail++;
				}
				else{
					pStats->nTrimAvail++;
				}
				mLen = (pos + pos2) / 2;
				pStats->incrementCount(size_t(mLen));
			}
			pTaskMan->decreaseCnt();
			startId += task.nBlockSize;
			pRecord = &pBuffer[(startId << 1) % size2];
			nItemCnt = pTaskMan->getItemCnt(startId, pRecord);
		}while(nItemCnt > 0);
	}
	return NULL;
}

void * mt_worker2_mp(void * data)
{
	cData * pData = (cData *)data;
	cTaskManager *pTaskMan = pData->pTaskMan;
	cStats * pStats = pData->pStats;
	int64 file_length = pStats->total_file_length;
	cFQ * pfq = pStats->pfq;
	cFQ * pfq2 = pStats->pfq2;
	FILE *fpOut = pStats->fpOut;
	FILE *fpOut2 = pStats->fpOut2;
	FILE *fpMask = pStats->fpMask;
	FILE *fpMask2 = pStats->fpMask2;
	FILE *fpExcl = pStats->fpExcl;
	FILE *fpExcl2 = pStats->fpExcl2;
	int minAverageQual = pStats->minAverageQual;
	int minEndQual = pStats->minEndQual;
	int minLen = pStats->getMinLen();
	int maxLen = pStats->getMaxLen();
	int maxLen2 = (maxLen == INT_MAX) ? INT_MAX : (maxLen * 2);
	bool bBarcode = pStats->bBarcode;
	bool bCutTail = pStats->bCutTail;
	bool bRedistribute = pStats->bRedistribute;

	RECORD *pBuffer, *pRecord, *pRecord2;
	TASK task;
	int size2, rc, rc2, nItemCnt, nCnt;
	int64 startId;

	pBuffer = pData->pBuffer;
	size2 = pData->size * 2;
	rc = rc2 = 0;

	int64 cur_pos;
	double cur_ratio;
	int pos, pos2;

	int rLen, qLen;
	int rLen2, qLen2;

	while(true){
		while(!pTaskMan->getTask(task)){
			if(pTaskMan->IsFinished()){
				task.type = TASK_END;
				break;
			}
			usleep(1);
		}
		if(task.type == TASK_END){
			break;
		}
		startId = task.startId;
		if(task.type == TASK_READ){
			if(!pTaskMan->increaseCnt()){ // reach the buffer size
				pTaskMan->addTask(task); // perform reading later
				usleep(1);
				continue;
			}
			// read records from input file to buffer
			for(pRecord=&pBuffer[(startId << 1) % size2], nItemCnt=0; nItemCnt<task.nBlockSize; nItemCnt++, pRecord+=2){
				pRecord2 = pRecord+1;
				rc = pfq->readRecord(pRecord);
				rc2 = pfq2->readRecord(pRecord2);
				if( (rc < 0) || (rc2 < 0) ){
					break;
				}
			}
			if(!pStats->bQuiet){
				cur_pos = pfq->tell() + pfq2->tell();
				if(cur_pos >= pfq->next_pos){
					cur_ratio = int64(cur_pos * 10000 / file_length) / 10000.0;
					pStats->progress(cur_ratio, 50);
					pfq->next_pos = int64(((cur_ratio * 10000 + 1) * file_length + 9999)/10000);
				}
			}
			if( (rc < 0) || (rc2 < 0) ){ // error or end of file
				pTaskMan->finish();
				if( (rc < -1) || (rc2 < -1) ) continue; // error
				if(nItemCnt == 0) continue; // no record read
			}
			task.startId += task.nBlockSize;
			pTaskMan->addTask(task); // save next task for parallelism

			// process the records
			for(pRecord=&pBuffer[(startId << 1) % size2], nCnt=0; nCnt < nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if(minAverageQual > 0){
					if( !cMatrix::checkQualities((uchar *)pRecord->qual.s, pRecord->qual.n, minAverageQual) &&
						!cMatrix::checkQualities((uchar *)pRecord2->qual.s, pRecord2->qual.n, minAverageQual) ){
						pRecord->tag = pRecord2->tag = TAG_BADQUAL;
						continue;
					}
				}
				pRecord->tag = TAG_NORMAL;
				rLen = pRecord->seq.n;
				qLen = pRecord->qual.n;
				rLen2 = pRecord2->seq.n;
				qLen2 = pRecord2->qual.n;
				if(cMatrix::findAdapterWithPE(pRecord->seq.s, pRecord2->seq.s, rLen, rLen2,
							(uchar *)pRecord->qual.s, (uchar *)pRecord2->qual.s, qLen, qLen2,
						 	pRecord->idx, pRecord2->idx)){ // trimmed
					pos = pRecord->idx.pos;
					pos2 = pRecord2->idx.pos;
					if( (pos >= minLen) && (pos2 >= minLen) ){
						cMatrix::combinePairSeqs(pRecord->seq.s, pRecord2->seq.s, pos, pos2,
							(uchar *)pRecord->qual.s, (uchar *)pRecord2->qual.s, (int)qLen, (int)qLen2);
						if(pStats->bFilterNs){
							if( cMatrix::isBlurry(pRecord->seq.s, pos) && cMatrix::isBlurry(pRecord2->seq.s, pos2) ){
								pRecord->tag = pRecord2->tag = TAG_BLURRY;
							}
						}
					}
				}
				else{
					pos = rLen;
					pos2 = rLen2;
					if(pStats->bFilterNs){
						if( cMatrix::isBlurry(pRecord->seq.s, rLen) && cMatrix::isBlurry(pRecord2->seq.s, rLen) ){
							pRecord->tag = pRecord2->tag = TAG_BLURRY;
						}
					}
				}
				if( (pRecord->tag == TAG_NORMAL) && (pos >= minLen) && (pos2 >= minLen) ) {
					pRecord->idx = cMatrix::findJuncAdapter(pRecord->seq.s, pos, (uchar *)pRecord->qual.s, qLen);
					if(pRecord->idx.pos < 0){
						pRecord->idx.pos = 0;
					}
					pRecord2->idx = cMatrix::findJuncAdapter(pRecord2->seq.s, pos2, (uchar *)pRecord2->qual.s, qLen2);
					if(pRecord2->idx.pos < 0){
						pRecord2->idx.pos = 0;
					}
					if( (pos < rLen) || (pos2 < rLen2) ) { // trimmed
						if( (pRecord->idx.bc == 0) || (pRecord2->idx.bc == 0) || (pRecord->idx.bc != pRecord2->idx.bc) ){
							pRecord->tag = pRecord2->tag = TAG_CONTAMINANT;
						}
						else{
							if(pRecord->idx.pos + pRecord2->idx.pos + cMatrix::junctionLengths[pRecord->idx.bc] != max(pos, pos2)){
								pRecord->tag = pRecord2->tag = TAG_CONTAMINANT;
							}
						}
						if(minEndQual > 0){
							if(pRecord->qual.n > 0)
								pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->idx.pos, pRecord->qual.n), minEndQual);
							if(pRecord2->qual.n > 0)
								pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, min(pRecord2->idx.pos, pRecord2->qual.n), minEndQual);
						}
					}
					else{
						if( (pRecord->idx.bc == 0) && (pRecord2->idx.bc == 0) ){ // case D
							if(pStats->bFilterUndetermined)
								pRecord->tag = pRecord2->tag = TAG_UNDETERMINED;
							else{
								if(minEndQual > 0){
									if(pRecord->qual.n > 0)
										pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->qual.n, pRecord->idx.pos), minEndQual);
									if(pRecord2->qual.n > 0)
										pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, min(pRecord2->qual.n, pRecord2->idx.pos), minEndQual);
								}
							}
						}
						else{
							if(pRecord->idx.bc == 0){ // case B
								if( (pRecord2->idx.pos >= minLen) && (pRecord2->seq.n >= rLen) && (pRecord2->qual.n >= qLen) ){
									if(bRedistribute){
										pRecord->idx = cMatrix::mergePE(pRecord->seq.s, pRecord2->seq.s, rLen, (uchar *)pRecord->qual.s, (uchar *)pRecord2->qual.s, qLen, pRecord2->idx.pos, cMatrix::junctionLengths[pRecord2->idx.bc]);
									}
									if(minEndQual > 0){
										int pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s + pRecord2->idx.pos, pRecord->idx.pos - pRecord->seq.n, minEndQual);
										if(pos != pRecord->idx.pos - pRecord->seq.n){
											if(pos == 0){
												pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, pRecord->qual.n, minEndQual);
											}
											else{
												pRecord->idx.pos = pRecord->seq.n + pos;
											}
										}
									}
								}
							}
							else if(pRecord2->idx.bc == 0){ // case C
								if( (pRecord->idx.pos >= minLen) && (pRecord->seq.n >= rLen) && (pRecord->qual.n >= qLen) ){
									if(bRedistribute){
										pRecord2->idx = cMatrix::mergePE(pRecord2->seq.s, pRecord->seq.s, rLen, (uchar *)pRecord2->qual.s, (uchar *)pRecord->qual.s, qLen, pRecord->idx.pos, cMatrix::junctionLengths[pRecord->idx.bc]);
									}
									if(minEndQual > 0){
										int pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s + pRecord->idx.pos, pRecord2->idx.pos - pRecord2->seq.n, minEndQual);
										if(pos != pRecord2->idx.pos - pRecord2->seq.n){
											if(pos == 0){
												pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, pRecord2->qual.n, minEndQual);
											}
											else{
												pRecord2->idx.pos = pRecord2->seq.n + pos;
											}
										}
									}
								}
							} // case C
							else{
								if(minEndQual > 0){
									if(pRecord->qual.n > 0)
										pRecord->idx.pos = cMatrix::trimByQuality((uchar *)pRecord->qual.s, min(pRecord->qual.n, pRecord->idx.pos), minEndQual);
									if(pRecord2->qual.n > 0)
										pRecord2->idx.pos = cMatrix::trimByQuality((uchar *)pRecord2->qual.s, min(pRecord2->qual.n, pRecord2->idx.pos), minEndQual);
								}
							}
						} // not case D
					} // pos >= rLen
				}
			}

			pRecord = &pBuffer[(startId << 1) % size2];
			if(!pTaskMan->setItemCnt(startId, pRecord, nItemCnt)){
				task.type = TASK_WRITE;
				task.startId = startId;
				task.nItemCnt = nItemCnt;
				if(pTaskMan->bSingleBlock)
					pTaskMan->insertTask(task);
				else
					pTaskMan->addTask(task);
			}
			continue;
		}
		// task.type == TASK_WRITE
		pRecord = &pBuffer[(startId << 1) % size2];
		nItemCnt = task.nItemCnt;
		do{
			// write to file
			pRecord->nCnt = 0; // reset
			for(nCnt=0; nCnt<nItemCnt; nCnt++, pRecord+=2){
				pRecord2 = pRecord + 1;
				if(pRecord->tag == TAG_BLURRY){
					pStats->nBlurry++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pRecord->tag == TAG_BADQUAL){
					pStats->nBad++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pRecord->tag == TAG_CONTAMINANT){
					pStats->nContaminant++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pRecord->tag == TAG_UNDETERMINED){
					pStats->nUndetermined++;
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				// TAG_NORMAL
				pos = pRecord->idx.pos;
				pos2 = pRecord2->idx.pos;
				if( (pos < minLen) || (pos2 < minLen) ){
					if( (pos <= 0) || (pos2 <= 0) ) {
						pStats->nEmpty++;
						pRecord->tag = pRecord2->tag = TAG_EMPTY;
					}
					else {
						pStats->nShort++;
						pRecord->tag = pRecord2->tag = TAG_SHORT;
					}
					if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
						OutputTaggedRecord(fpExcl, pRecord);
						OutputTaggedRecord(fpExcl2, pRecord2);
					}
					continue;
				}
				if(pos + pos2 > maxLen2){
					if(!bCutTail){
						pStats->nLong++;
						pRecord->tag = pRecord2->tag = TAG_LONG;
						if( (fpExcl != NULL) && (fpExcl2 != NULL) ){
							OutputTaggedRecord(fpExcl, pRecord);
							OutputTaggedRecord(fpExcl2, pRecord2);
						}
						continue;
					}
					if(pos > maxLen) pos = maxLen;
					if(pos2 > maxLen) pos2 = maxLen;
				}
				
				if(bBarcode){
					if(pRecord->idx.bc < 0){
						fpOut = pStats->fpUntrim.fp;
						fpOut2 = pStats->fpUntrim2.fp;
					}
					else{
						fpOut = pStats->fpOuts[pRecord->idx.bc].fp;
						fpOut2 = pStats->fpOuts2[pRecord->idx.bc].fp;
						pStats->incrementBarcode(pRecord->idx.bc);
					}
				}
				rLen = pRecord->seq.n;
				if( (fpMask != NULL) && (fpMask2 != NULL) ){
					OutputMaskedRecord(fpMask, pRecord, 0, pos);
					OutputMaskedRecord(fpMask2, pRecord2, 0, pos2);
				}
				if(pRecord->com.n > 0){ // fastq
					if(pos <= rLen){
						fprintf(fpOut, "@%s%.*s\n+\n%.*s\n", pRecord->id.s, pos, pRecord->seq.s, pos, pRecord->qual.s);
					}
					else{
						fprintf(fpOut, "@%s%.*s%.*s\n", pRecord->id.s, rLen, pRecord->seq.s, pos - rLen, pRecord2->seq.s + pos2);
						fprintf(fpOut, "+\n%.*s%.*s\n", rLen, pRecord->qual.s, pos - rLen, pRecord2->qual.s + pos2);
					}
					if(pos2 <= rLen){
						fprintf(fpOut2, "@%s%.*s\n+\n%.*s\n", pRecord2->id.s, pos2, pRecord2->seq.s, pos2, pRecord2->qual.s);
					}
					else{
						fprintf(fpOut2, "@%s%.*s%.*s\n", pRecord2->id.s, rLen, pRecord2->seq.s, pos2 - rLen, pRecord->seq.s + pos);
						fprintf(fpOut2, "+\n%.*s%.*s\n", rLen, pRecord2->qual.s, pos2 - rLen, pRecord->qual.s + pos);
					}
				}
				else{ // fasta
					if(pos <= rLen){
						fprintf(fpOut, ">%s%.*s\n", pRecord->id.s, pos, pRecord->seq.s);
					}
					else{
						fprintf(fpOut, ">%s%.*s%.*s\n", pRecord->id.s, rLen, pRecord->seq.s, pos - rLen, pRecord2->seq.s + pos2);
					}
					if(pos2 <= rLen){
						fprintf(fpOut2, ">%s%.*s\n", pRecord2->id.s, pos2, pRecord2->seq.s);
					}
					else{
						fprintf(fpOut2, ">%s%.*s%.*s\n", pRecord2->id.s, rLen, pRecord2->seq.s, pos2 - rLen, pRecord->seq.s + pos);
					}
				}
				if(bBarcode){
					if(pRecord->idx.bc < 0){ // assigned
						pStats->nUntrimAvail++;
					}
					else{
						pStats->nTrimAvail++;
					}
				}
				else{
					if(pos + pos2 < rLen + rLen) // trimmed
						pStats->nTrimAvail++;
					else
						pStats->nUntrimAvail++;
				}
				pStats->incrementCount(size_t((pos + pos2) / 2));
			}
			pTaskMan->decreaseCnt();
			startId += task.nBlockSize;
			pRecord = &pBuffer[(startId << 1) % size2];
			nItemCnt = pTaskMan->getItemCnt(startId, pRecord);
		}while(nItemCnt > 0);
	}
	return NULL;
}

int processFile(cParameter * pParameter, cStats * pStats)
{
	CFILE cf;
	int i;

	int64 file_length;
	if(pParameter->bStdin){
		cf.fp = stdin;
		file_length = -1;
	}
	else{
		char * inFile = pParameter->input[0];
		file_length = gzsize(inFile);
		cf = gzopen(inFile, "r");
		if(cf.fp == NULL){
			fprintf(stderr, "Can not open %s for reading\n", inFile);
			return 1;
		}
	}
	cWork wk;
	if(!wk.Init(pParameter, pStats, file_length, cf.fp)){
		fprintf(stderr, "Can not allocate memory for workset\n");
		gzclose(&cf);
		return 1;
	}
	mtaux_t *mt = wk.getMultiThreadingPointer();

	int rc;
	void *status;
	if( (pParameter->trimMode & TRIM_AP) != 0 ){
		for(i=1; i<mt->n_threads; i++){ // worker 0 is effectively launched by the master thread
			rc = pthread_create(&mt->tid[i], &wk.attr, mt_worker_ap, &mt->w[i]);
			if(rc != 0){
				fprintf(stderr, "Can not create thread %d\n", i);
				break;
			}
		}
		mt_worker_ap(&mt->w[0]);
	}
	else{
		for(i=1; i<mt->n_threads; i++){ // worker 0 is effectively launched by the master thread
			rc = pthread_create(&mt->tid[i], &wk.attr, mt_worker, &mt->w[i]);
			if(rc != 0){
				fprintf(stderr, "Can not create thread %d\n", i);
				break;
			}
		}
		mt_worker(&mt->w[0]);
	}
	for(i=1; i<mt->n_threads; ++i){ // waits for termination of other threads
		rc = pthread_join(mt->tid[i], &status);
	}
	if(!pParameter->bStdin){
		gzclose(&cf);
	}
	return 0;
}

int processPairedFiles(cParameter * pParameter, cStats * pStats)
{
	char * inFile = pParameter->input[0];
	char * inFile2 = pParameter->input[1];
	CFILE cf, cf2;
	int i;

	int64 file_length = gzsize(inFile) + gzsize(inFile2);
	cf = gzopen(inFile, "r");
	cf2 = gzopen(inFile2, "r");
	if( (cf.fp == NULL) || (cf2.fp == NULL) ){
		if(cf.fp == NULL)
			fprintf(stderr, "Can not open %s for reading\n", inFile);
		else
			gzclose(&cf);
		if(cf2.fp == NULL)
			fprintf(stderr, "Can not open %s for reading\n", inFile2);
		else
			gzclose(&cf2);
		return 1;
	}

	cWork wk;
	if(!wk.Init(pParameter, pStats, file_length, cf.fp, cf2.fp)){
		fprintf(stderr, "Can not allocate memory for workset\n");
		gzclose(&cf2);
		gzclose(&cf);
		return 1;
	}
	mtaux_t *mt = wk.getMultiThreadingPointer();
	int rc;
	void *status;
	if( (pParameter->trimMode & TRIM_PE) != 0 ){
		if( (pParameter->trimMode & TRIM_ANY) == TRIM_DEFAULT ){
			for(i=1; i<mt->n_threads; i++){ // worker 0 is effectively launched by the master thread
				rc = pthread_create(&mt->tid[i], &wk.attr, mt_worker2, &mt->w[i]);
				if(rc != 0){
					fprintf(stderr, "Can not create thread %d\n", i);
					break;
				}
			}
			mt_worker2(&mt->w[0]);
		}
		else{
			for(i=1; i<mt->n_threads; i++){ // worker 0 is effectively launched by the master thread
				rc = pthread_create(&mt->tid[i], &wk.attr, mt_worker2_sep, &mt->w[i]);
				if(rc != 0){
					fprintf(stderr, "Can not create thread %d\n", i);
					break;
				}
			}
			mt_worker2_sep(&mt->w[0]);
		}
	}
	else if( (pParameter->trimMode & TRIM_AP) != 0 ){
		for(i=1; i<mt->n_threads; i++){ // worker 0 is effectively launched by the master thread
			rc = pthread_create(&mt->tid[i], &wk.attr, mt_worker2_ap, &mt->w[i]);
			if(rc != 0){
				fprintf(stderr, "Can not create thread %d\n", i);
				break;
			}
		}
		mt_worker2_ap(&mt->w[0]);
	}
	else{ // TRIM_MP
		for(i=1; i<mt->n_threads; i++){ // worker 0 is effectively launched by the master thread
			rc = pthread_create(&mt->tid[i], &wk.attr, mt_worker2_mp, &mt->w[i]);
			if(rc != 0){
				fprintf(stderr, "Can not create thread %d\n", i);
				break;
			}
		}
		mt_worker2_mp(&mt->w[0]);
	}
	for(i=1; i<mt->n_threads; ++i){ // waits for termination of other threads
		rc = pthread_join(mt->tid[i], &status);
	}
	gzclose(&cf2);
	gzclose(&cf);
	return 0;
}

int main(int argc, char * argv[])
{
	cParameter para;
	cStats stats;
	char errMsg[256];
	// process the input parameters
	int iRet = para.GetOpt(argc, argv, errMsg);
	if(iRet < 0){
		char * program = strrchr(argv[0], '/');
		program = (program == NULL) ? argv[0] : (program + 1);
		if(iRet == -1){
			if(para.bEnquireVersion){
				para.PrintVersion(stdout);
				return 0;
			}
			para.PrintUsage(program, stdout);
		}
		else{
			fprintf(stderr, "%s (%s): %s\n\n", program, para.version, errMsg);
			para.PrintSimpleUsage(program, stderr);
		}
		return 1;
	}
	if(para.IsAutoFastqFormat()){
		para.fastqFormat = gzformat(para.input, para.nFileCnt);
		if(para.fastqFormat == CONTRADICT_FASTQ){
			fprintf(stderr, "Error: the FASTQ quality formats of input files are different\n");
			return 1;
		}
		para.baseQual = (para.fastqFormat == SOLEXA_FASTQ) ? 64 : 33;
	}
	if(!stats.initHist(&para)){
		fprintf(stderr, "Error: can not allocate memory for audit\n");
		return 1;
	}
	if(!stats.openOutputFiles(&para)){
		return 1;
	}
	FILE * hLog;
	if(para.bStdout){
		hLog = stderr;
	}
	else{
		hLog = fopen(para.logfile, "w");
		if(hLog == NULL){
			fprintf(stderr, "Error: can not open %s for writing\n", para.logfile);
			return 1;
		}
	}
	para.printVersion(hLog);
	para.printCommandLine(hLog);
	para.printRelatedFiles(hLog);

	para.printOpt(hLog, true);
	if(!stats.bStdout){
		para.printLogo(stdout);
		para.printVersion(stdout);
		para.printOpt(stdout);
	}

	stats.printTime("started", hLog);
	if(!stats.bStdout) stats.printTime("started", stdout);
	stats.start();

	////////////// process the input file(s)
	if(para.nFileCnt <= 1){
		iRet = processFile(&para, &stats);
	}
	else{
		iRet = processPairedFiles(&para, &stats);
	}
	if(iRet != 0){
		if(!stats.bStdout) fclose(hLog);
		return iRet;
	}

	stats.end();
	stats.printTime("done", hLog, 0x02);
	if(!stats.bStdout) stats.printTime("done", stdout, 0x02);
	stats.printDiffTime(hLog);
	if(!stats.bStdout) stats.printDiffTime(stdout);
	stats.printSummary(hLog);
	if(!stats.bStdout) stats.printSummary(stdout);

	stats.printHist(hLog);
	if(!stats.bStdout){
		fclose(hLog);
		fprintf(stdout, "log has been saved to \"%s\".\n", para.logfile);
	}
	if(!stats.writeMapFile(&para)){
		fprintf(stderr, "Can not write Mapping file\n");
		return 1;
	}

	return 0;
}
