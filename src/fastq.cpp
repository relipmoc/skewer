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
#include "fastq.h"

const char * FASTQ_FORMAT_NAME[FASTQ_FORMAT_CNT] = {
	"Sanger/Illumina 1.8+ FASTQ",
	"Solexa/Illumina 1.3+/Illumina 1.5+ FASTQ",
	"Unknown format",
	"Contradict formats",
	"FASTA"
};

////////////////////////
cFQ::cFQ()
{
	memset(&rec, 0, sizeof(RECORD));

	in = NULL;
	offset = 0L;
	next_pos = 0L;
	rno = 0;
}

cFQ::~cFQ()
{
	if(rec.id.s != NULL)
		free(rec.id.s);
	if(rec.seq.s != NULL)
		free(rec.seq.s);
	if(rec.com.s != NULL)
		free(rec.com.s);
	if(rec.qual.s != NULL)
		free(rec.qual.s);
}

void cFQ::associateFile(FILE * fp)
{
	in = fp;
	offset = 0L;
	rno = 0;
}

inline int cFQ::read_line(LINE &l)
{
	return (l.n = getline(&l.s, &l.a, in));
}

// readRecord 
// return value:
// -1: EOF
// -2: ERROR
//  1: sucess
int cFQ::readRecord(RECORD *pRecord)
{
	if(pRecord == NULL){
		pRecord = &rec;
	}
	tag = fgetc(in);
	if(tag == EOF)
		return -1;
	read_line(pRecord->id);
	if(tag == '>') {
		pRecord->qual.n = 0;
		// read fasta instead
		char c = fgetc(in);
		pRecord->seq.n = 0;
		pRecord->com.n = 0;
		while (c != '>' && c != EOF) {
			if (pRecord->seq.a <= size_t(pRecord->seq.n+1)) {
				pRecord->seq.s=(char *)realloc(pRecord->seq.s, pRecord->seq.a=(pRecord->seq.a * 3 / 2 + 64));
			}
			if (!isspace(c)) 
				pRecord->seq.s[pRecord->seq.n++]=c;
			c = fgetc(in);
		}
		if (c != EOF) {
			ungetc(c, in);
		}
	}
	else{
		read_line(pRecord->seq);
		read_line(pRecord->com);
		read_line(pRecord->qual);
	}
	rno++;
	offset += 1 + pRecord->id.n + pRecord->seq.n + pRecord->com.n + pRecord->qual.n;
	if(tag == '>'){
		if(pRecord->seq.n <= 0){
			fprintf(stderr, "Malformed fasta record %d: empty sequence\n", rno);
			return -2;
		}
		return 1;
	}

	if (tag != '@' || pRecord->com.s[0] != '+' || pRecord->seq.n != pRecord->qual.n) {
		const char *errtyp = (pRecord->seq.n != pRecord->qual.n) ?  "length mismatch" :
			 pRecord->id.s[0] != '@' ? "no '@' for id" : "no '+' for comment";
		fprintf(stderr, "Malformed fastq record %d: %s\n", rno, errtyp);
		return -2;
	}
	// win32-safe chomp
	pRecord->seq.s[--pRecord->seq.n] = '\0';
	if (pRecord->seq.s[pRecord->seq.n-1] == '\r') {
		pRecord->seq.s[--pRecord->seq.n] = '\0';
	}
	pRecord->qual.s[--pRecord->qual.n] = '\0';
	if (pRecord->qual.s[pRecord->qual.n-1] == '\r') {
		pRecord->qual.s[--pRecord->qual.n] = '\0';
	}
	return 1;
}

///////////////////
// subroutines
const char *fext(const char *f) {
	const char *x = strrchr(f,'.');
	return (x != NULL) ? (x + 1) : "";
}

///////////////////////
// external functions
CFILE gzopen(const char * fileName, const char * mode)
{
	// maybe use zlib some day?
	CFILE cf;
	const char * ext = fext(fileName);
	if(strcmp(mode, "r") == 0){
		cf.fp = fopen(fileName, "r");
		if(cf.fp == NULL){
			return cf;
		}
		fclose(cf.fp);
	}
	if (strcmp(ext,"gz") == 0) {
		char *tmp=(char *)malloc(strlen(fileName)+100);
		if (strchr(mode, 'w')) {
			strcpy(tmp, "gzip --rsyncable > '");
			strcat(tmp, fileName);
			strcat(tmp, "'");
		} else {
			strcpy(tmp, "gunzip -c '");
			strcat(tmp, fileName);
			strcat(tmp, "'");
		}
		cf.fp = popen(tmp, mode);
		cf.bGz = true;
		free(tmp);
	} else if (strcmp(ext,"zip") == 0) {
		char *tmp=(char *)malloc(strlen(fileName)+100);
		if (strchr(mode,'w')) {
			strcpy(tmp, "zip -q '");
			strcat(tmp, fileName);
			strcat(tmp, "' -");
		} else {
			strcpy(tmp, "unzip -p '");
			strcat(tmp, fileName);
			strcat(tmp, "'");
		}
		cf.fp = popen(tmp, mode);
		cf.bGz = true;
		free(tmp);
	} else {
		cf.fp = fopen(fileName, mode);
		cf.bGz = false;
	}
	return cf;
}

int gzclose(CFILE *f)
{
	if( (f == NULL) || (f->fp == NULL) )
		return -1;
	int iRet = (f->bGz ? pclose(f->fp) : fclose(f->fp));
	f->fp = NULL;
	return iRet;
}

int64 gzsize(const char * fileName)
{
	FILE * fp = fopen(fileName, "rb");
	if(fp == NULL)
		return 0L;
	size_t len = strlen(fileName);
	bool bCompressed = (len > 3) && (strcmp(fileName + len - 3, ".gz") == 0);
	int64 file_length = 0L;
	if(bCompressed){
		if(fseek(fp, -4L, SEEK_END) == 0){
			int64 compress_length;
			unsigned char buffer[4];
			uint32 x;
			fread(buffer, 1, 4, fp);
			x = buffer[3];
			x <<= 8;
			x |= buffer[2];
			x <<= 8;
			x |= buffer[1];
			x <<= 8;
			x |= buffer[0];
			compress_length = ftell(fp);
			file_length = x;
			if(file_length < 2 * compress_length){
				file_length += (((2 * compress_length) - file_length + (1L << 32) - 1) / (1L << 32)) * (1L << 32);
			}
		}
	}
	else{
		if(fseek(fp, 0, SEEK_END) == 0){
			file_length = ftell(fp);
		}
	}
	fclose(fp);
	return file_length;
}

int gzreadlen(char * fileName)
{
	int nReadLen = 0;
	CFILE cf;
	cFQ fq;

	cf = gzopen(fileName, "r");
	if(cf.fp == NULL){
		fprintf(stderr, "Can not open %s for reading\n", fileName);
		return -1;
	}
	fq.associateFile(cf.fp);
	int iRet;
	for(int i=0; i<100; i++){
		iRet = fq.readRecord();
		if(iRet <= 0){
			if(iRet < -1){ // error
				nReadLen = -1;
			}
			break;
		}
		if(fq.rec.seq.n > nReadLen){
			nReadLen = fq.rec.seq.n;
		}
	}
	gzclose(&cf);
	return nReadLen;
}

enum FASTQ_FORMAT gzformat(char * fileNames[], int nFileCnt)
{
	CFILE cf;
	cFQ fq;
	int i, j;
	char *str;
	char chr;

	FASTQ_FORMAT format = UNKNOWN_FASTQ, format_new;
	for(i=0; i<nFileCnt; i++){
		cf = gzopen(fileNames[i], "r");
		if(cf.fp == NULL) break;
		fq.associateFile(cf.fp);
		format_new = UNKNOWN_FASTQ;
		while(fq.readRecord() > 0){
			if(fq.rec.com.n == 0){
				format_new = FASTA;
				break;
			}
			for(j=0, str=fq.rec.qual.s; j<fq.rec.qual.n; j++){
				chr = *(str++);
				if(chr < 59){
					format_new = SANGER_FASTQ;
					break;
				}
				if(chr > 74){
					format_new = SOLEXA_FASTQ;
					break;
				}
			}
			if(j < fq.rec.qual.n) break;
		}
		gzclose(&cf);
		if(format == UNKNOWN_FASTQ){
			format = format_new;
		}
		else{
			if(format_new != format && format_new != UNKNOWN_FASTQ){
				format = CONTRADICT_FASTQ;
				break;
			}
		}
	}

	return format;
}

void gzstrncpy (char * dest, const char * src, int n)
{
	while( (*dest++ = *src++) && --n );
	*dest = '\0';
}

void versatile_gettime(struct timespec *tp)
{
#ifdef __MACH__
	static double orwl_timebase = 0.0;
	static uint64_t orwl_timestart = 0; 
	if(!orwl_timestart){
		mach_timebase_info_data_t tb = { 0 };
		mach_timebase_info(&tb);
		orwl_timebase = tb.numer;
		orwl_timebase /= tb.denom;
		orwl_timestart = mach_absolute_time();
	}
	struct timespec t;
	double diff = (mach_absolute_time() - orwl_timestart) * orwl_timebase;
	tp->tv_sec = diff * ORWL_NANO;
	tp->tv_nsec = diff - (t.tv_sec * ORWL_GIGA);
#else
	clock_gettime(CLOCK_MONOTONIC, tp);
#endif
}
