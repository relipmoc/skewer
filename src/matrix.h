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
#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdlib.h>
#include <deque>
#include <vector>
#include <set>
#include <string>
#include "common.h"

using namespace std;

typedef struct{
	double score;
	int nIndel;
	INDEX idx;
}ELEMENT;

class ElementComparator
{
public:
	bool operator()(const ELEMENT &elem1, const ELEMENT &elem2){
		return elem1.idx.pos < elem2.idx.pos;
	}
};

typedef set<ELEMENT, ElementComparator> ELEMENT_SET;

class cElementSet : public ELEMENT_SET
{
public:
	bool insert(const ELEMENT& val);
};

typedef enum{
	CD_NONE = 0,
	CD_A = 1, // Adenosine
	CD_C = 2, // Cytidine
	CD_G = 3, // Guanosine
	CD_T = 4, // Thymidine (or Uridine)
	CD_R = 5, // puRine, A or G
	CD_Y = 6, // pYrimidine, T or C
	CD_S = 7, // Strong, G or C
	CD_W = 8, // Weak, A or T
	CD_K = 9, // Keto, G or T
	CD_M = 10, // aMino, A or C
	CD_B = 11, // not A
	CD_D = 12, // not C
	CD_H = 13, // not G
	CD_V = 14, // not T
	CD_N = 15, // any base
	CD_CNT = CD_N+1,
	CD_BASIC_CNT = 5
}CODE;

class cAdapter
{
	char sequence[MAX_ADAPTER_LEN+1]; // for debug only
	char barcode[MAX_ADAPTER_LEN+1];
	char primer[MAX_ADAPTER_LEN+1];
	bool masked[MAX_ADAPTER_LEN+1];
	inline void UPDATE_COLUMN(deque<ELEMENT> & queue, uint64 &d0bits, uint64 &lbits, uint64 &unbits, uint64 &dnbits, double &penal, double &dMaxPenalty, int &iMaxIndel);

public:
	size_t len;
	TRIM_MODE trimMode;
	bool bBestAlign;
	uint64 matchBits[CD_CNT];

public:
	cAdapter();
	~cAdapter();
	void Init(char * seq, size_t sLen, TRIM_MODE trimMode);
	void Init2(char * seq, size_t sLen);
	bool align(char * read, size_t rLen, uchar * qual, size_t qLen, cElementSet &result, int bc, bool bBestAlign=true);

public:
	void initBarcode(int iCut);
	char * getBarcode() { return barcode; }
	char * getPrimer() { return primer; }
	bool * getMasked() { return masked; }
};

///////////////////////////////////////
class cMatrix
{
	friend class cAdapter;

	static bool bShareAdapter;

	static double dEpsilon, dEpsilonIndel;
	static double dPenaltyPerErr;
	static double dDelta, dMu;
	static double penalty[256];
	static bool bSensitive;

public:
	static vector<bool *> fw_masked;
	static vector<bool *> rv_masked;
	static vector<string> fw_barcodes;
	static vector<string> rv_barcodes;
	static vector<string> fw_primers;
	static vector<string> rv_primers;
	static vector<int> rowBc;
	static vector<int> colBc;

public:
	static deque<cAdapter> firstAdapters;
	static deque<cAdapter> secondAdapters;
	static deque<cAdapter> junctionAdapters;

	static vector<int> junctionLengths;
	static vector< vector<int> > indices;
	static int iIdxCnt;
	static int iMinOverlap;

public:
	cMatrix();
	~cMatrix();

private:
	static bool CalcRevCompScore(char * seq, char * seq2, int len, uchar * qual, uchar * qual2, size_t qLen, double &score);
	static string GetRevComp(char * seq, int len);

public:
	static void InitParameters(enum TRIM_MODE trimMode, double dEpsilon, double dEpsilonIndel, int baseQual, bool bShareAdapter);
	static void AddAdapter(deque<cAdapter> & adapters, char * vector, size_t len, TRIM_MODE trimMode);
	static void CalculateJunctionLengths();
	static void CalculateIndices(vector< vector<bool> > &bMatrix, int nRow, int nCol);
	static void InitBarcodes(deque<cAdapter> & fw_primers, int iCutF, deque<cAdapter> & rv_primers, int iCutR);

	static bool isBlurry(char * seq, size_t len);
	static bool checkQualities(uchar * quals, size_t len, int minQual);
	static int trimByQuality(uchar * quals, size_t len, int minQual);

	static INDEX findAdapter(char * read, size_t rLen, uchar * qual, size_t qLen);
	static INDEX findAdapter2(char * read, size_t rLen, uchar * qual, size_t qLen);
	static INDEX findJuncAdapter(char * read, size_t rLen, uchar * qual, size_t qLen);

	static bool findAdapterWithPE(char * read, char * read2, size_t rLen, size_t rLen2, uchar * qual, uchar * qual2, size_t qLen, size_t qLen2, INDEX &index, INDEX & index2);
	static int findAdaptersBidirectionally(char * read, size_t rLen, uchar * qual, size_t qLen,
			char * read2, size_t rLen2, uchar * qual2, size_t qLen2, INDEX &index, INDEX &index2);
	static int findAdaptersInARead(char * read, size_t rLen, uchar * qual, size_t qLen, INDEX &index);
	static bool PrepareBarcode(char * barcodeSeq, int bcIdx, char * seq, int len, char * seq2, int len2, char * barcodeQual, char * qual, char * qual2);
	static bool PrepareBarcode(char * barcodeSeq, int bcIdx, char * seq, int len, char * seq2, int len2);
	static INDEX mergePE(char * read, char * read2, size_t rLen, uchar * qual, uchar * qual2, size_t qLen, size_t startPos, size_t jLen);
	static bool combinePairSeqs(char * read, char * read2, int len, int len2, uchar * qual, uchar * qual2, int qLen, int qLen2);
};

#endif // _MATRIX_H
