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
#include <stdio.h>
#include <float.h>
#include <algorithm>
#include <string.h>
#include "matrix.h"
#include "fastq.h"

CODE codeMap[256] = {
	//  0        1        2        3        4        5        6        7        8        9        A        B        C        D        E        F
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 0
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 1
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 2
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 3
	// 0x41~0x5A, A~Z
	CD_NONE,    CD_A,    CD_B,    CD_C,    CD_D, CD_NONE, CD_NONE,    CD_G,    CD_H, CD_NONE, CD_NONE,    CD_K, CD_NONE,    CD_M,    CD_N, CD_NONE, // 4
	CD_NONE, CD_NONE,    CD_R,    CD_S,    CD_T,    CD_T,    CD_V,    CD_W, CD_NONE,    CD_Y, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 5
	// 0x61~0x7A, a~z
	CD_NONE,    CD_A,    CD_B,    CD_C,    CD_D, CD_NONE, CD_NONE,    CD_G,    CD_H, CD_NONE, CD_NONE,    CD_K, CD_NONE,    CD_M,    CD_N, CD_NONE, // 6
	CD_NONE, CD_NONE,    CD_R,    CD_S,    CD_T,    CD_T,    CD_V,    CD_W, CD_NONE,    CD_Y, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 7

	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 8
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // 9
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // A
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // B
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // C
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // D
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, // E
	CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE, CD_NONE  // F
};

bool blurry[256] = {
	// 0     1     2     3     4     5     6     7     8     9     A     B     C     D     E     F
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 0
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 1
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 2
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 3

	true,false, true,false, true, true, true,false, true, true, true, true, true, true, true, true, // 4
	true, true, true, true,false,false, true, true, true, true, true, true, true, true, true, true, // 5

	true,false, true,false, true, true, true,false, true, true, true, true, true, true, true, true, // 6
	true, true, true, true,false,false, true, true, true, true, true, true, true, true, true, true, // 7

	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 8
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // 9
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // A
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // B
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // C
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // D
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, // E
	true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true  // F
};

CODE complement[CD_CNT] = {
	CD_NONE, CD_T, CD_G, CD_C, CD_A, CD_Y, CD_R, CD_W, CD_S, CD_M, CD_K, CD_V, CD_H, CD_D, CD_B, CD_N
};

char character[CD_CNT] = {
	'N', 'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N'
};

double scoring[CD_CNT][CD_CNT] = {
	//   -    A    C    G    T  | R    Y    S    W    K    M |  B    D    H    V |  N
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1, 0.05}, //-
	{    1,   0,   1,   1,   1,   0,   1,   1,   0,   1,   0,   1,   0,   0,   0, 0.05}, //A
	{    1,   1,   0,   1,   1,   1,   0,   0,   1,   1,   0,   0,   1,   0,   0, 0.05}, //C
	{    1,   1,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   0,   1,   0, 0.05}, //G
	{    1,   1,   1,   1,   0,   1,   0,   1,   0,   0,   1,   0,   0,   0,   1, 0.05}, //T

	{    1,   0,   1,   0,   1,   0,   1,0.75,0.75,0.75,0.75, 0.5,   0, 0.5,   0, 0.05}, //R
	{    1,   1,   0,   1,   0,   1,   0,0.75,0.75,0.75,0.75,   0, 0.5,   0, 0.5, 0.05}, //Y
	{    1,   1,   0,   0,   1,0.75,0.75,   0,   1,0.75,0.75,   0, 0.5, 0.5,   0, 0.05}, //S
	{    1,   0,   1,   1,   0,0.75,0.75,   1,   0,0.75,0.75, 0.5,   0,   0, 0.5, 0.05}, //W
	{    1,   1,   1,   0,   0,0.75,0.75,0.75,0.75,   0,   1,   0,   0, 0.5, 0.5, 0.05}, //K
	{    1,   0,   0,   1,   1,0.75,0.75,0.75,0.75,   1,   0, 0.5, 0.5,   0,   0, 0.05}, //M

	{    1,   1,   0,   0,   0, 0.5,   0,   0, 0.5,   0, 0.5,   0, 0.4, 0.4, 0.4, 0.05}, //B
	{    1,   0,   1,   0,   0,   0, 0.5, 0.5,   0,   0, 0.5, 0.4,   0, 0.4, 0.4, 0.05}, //D
	{    1,   0,   0,   1,   0, 0.5,   0, 0.5,   0, 0.5,   0, 0.4, 0.4,   0, 0.4, 0.05}, //H
	{    1,   0,   0,   0,   1,   0, 0.5,   0, 0.5, 0.5,   0, 0.4, 0.4, 0.4,   0, 0.05}, //V

	{ 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,    0}  //N
};

uint64 chrVadp[CD_CNT][CD_CNT] = {
// adp   -    A    C    G    T  | R    Y    S    W    K    M |  B    D    H    V |  N   //chr
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0}, //-
	{    1,   0,   1,   1,   1,   0,   1,   1,   0,   1,   0,   1,   0,   0,   0,   0}, //A
	{    1,   1,   0,   1,   1,   1,   0,   0,   1,   1,   0,   0,   1,   0,   0,   0}, //C
	{    1,   1,   1,   0,   1,   0,   1,   0,   1,   0,   1,   0,   0,   1,   0,   0}, //G
	{    1,   1,   1,   1,   0,   1,   0,   1,   0,   0,   1,   0,   0,   0,   1,   0}, //T

	{    1,   1,   1,   1,   1,   0,   1,   1,   1,   1,   1,   1,   0,   1,   0,   0}, //R
	{    1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   1,   0,   1,   0,   1,   0}, //Y
	{    1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0,   1,   1,   0,   0}, //S
	{    1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0,   0,   1,   0}, //W
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0,   0,   1,   1,   0}, //K
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   0,   0,   0}, //M

	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   1,   0}, //B
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   1,   0}, //D
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0}, //H
	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0}, //V

	{    1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0}  //N
};

const double MIN_PENALTY = 0.477121255;
const double MEAN_PENALTY = 2.477121255;
const double MAX_PENALTY = 4.477121255;
const double EPSILON = (MIN_PENALTY / 10);
///////////////////////////////////////

bool cElementSet::insert (const ELEMENT& val)
{
	pair<ELEMENT_SET::iterator,bool> ret;
	ELEMENT_SET::iterator it = this->find(val);
	if(it == this->end()){
		ret = ELEMENT_SET::insert(val);
		return ret.second;
	}
	if(val.score < it->score){
		return true;
	}
	this->erase(it++);
	ELEMENT_SET::insert(it, val);
	return true;
}

///////////////////////////////////////
cAdapter::cAdapter()
{
	len = 0;
}

cAdapter::~cAdapter()
{
}

void cAdapter::Init(char * seq, size_t sLen, TRIM_MODE trimMode)
{
	int i;
	// construct sequence
	this->len = (int(sLen) > MAX_ADAPTER_LEN) ? MAX_ADAPTER_LEN : sLen;
	gzstrncpy(sequence, seq, len);
	this->trimMode = trimMode;

	// construct mismatch bits
	int code, code2;
	uint64 bits;
	for(code=0; code<CD_CNT; code++){
		bits = 0;
		for(i=int(len)-1; i>=0; i--){
			code2 = codeMap[uchar(sequence[i])];
			bits = (bits << 1) | chrVadp[code][code2];
		}
		matchBits[code] = ~bits;
	}
}

void cAdapter::Init2(char * seq, size_t sLen)
{
	int i;
	// construct sequence
	if(int(sLen) > MAX_ADAPTER_LEN){
		// seq += (sLen - MAX_ADAPTER_LEN);
		// the head instead of the tail should be used
		this->len = MAX_ADAPTER_LEN;
	}
	else{
		this->len = sLen;
	}
	for(i=0; i<int(len); i++)
		sequence[i] = character[complement[codeMap[uchar(seq[len-1-i])]]];
	sequence[len] = '\0';
	this->trimMode = TRIM_TAIL;

	// construct mismatch bits
	int code, code2;
	uint64 bits;
	for(code=CD_BASIC_CNT-1; code>=0; code--){
		bits = 0;
		for(i=int(len)-1; i>=0; i--){
			code2 = codeMap[uchar(sequence[i])];
			bits = (bits << 1) | chrVadp[code][code2];
		}
		matchBits[code] = ~bits;
	}
	for(bits=~bits,code=CD_BASIC_CNT; code<CD_CNT; code++){
		matchBits[code] = bits;
	}
}

inline void cAdapter::UPDATE_COLUMN(deque<ELEMENT> & queue, uint64 &d0bits, uint64 &lbits, uint64 &unbits, uint64 &dnbits, double &penal, double &dMaxPenalty, int &iMaxIndel)
{
	int i;
	double score;
	uint64 bits = ~lbits | d0bits;
	for(bits>>=1,i=1; i<int(queue.size())-1; i++,bits>>=1){
		if((bits & 0x01) == 0){
			if(cMatrix::bSensitive){
				score = queue[i].score + (penal - cMatrix::dDelta);
				if( (queue[i-1].score < score) && (queue[i-1].nIndel < iMaxIndel) ){
					if( (queue[i+1].score < score) && (queue[i+1].nIndel < iMaxIndel) ){
						if(queue[i-1].score < queue[i+1].score){
							queue[i] = queue[i-1];
							dnbits |= (1L << (i-1));
						}
						else{
							queue[i] = queue[i+1];
							unbits |= (1L << (i+1));
						}
					}
					else{
						queue[i] = queue[i-1];
						dnbits |= (1L << (i-1));
					}
					queue[i].nIndel++;
				}
				else{
					if( (queue[i+1].score < score) && (queue[i+1].nIndel < iMaxIndel) ){
						queue[i] = queue[i+1];
						unbits |= (1L << (i+1));
						queue[i].nIndel++;
					}
					else{
						queue[i].score = score;
					}
				}
				queue[i].score += cMatrix::dDelta;
			}
			else{ // !cMatrix::bSensitive
				queue[i].score += penal;
			}
			if(queue[i].score >= dMaxPenalty){
				lbits &= ~(1L << i);
			}
		}
	}
	if(queue.size() > 1){
		if((bits & 0x01) == 0){
			if(cMatrix::bSensitive){
				if( (queue[i-1].nIndel < iMaxIndel) && (queue[i-1].score + cMatrix::dDelta < queue[i].score + penal) ){
					queue[i] = queue[i-1];
					dnbits |= (1L << (i-1));
					queue[i].score += cMatrix::dDelta;
					queue[i].nIndel++;
				}
				else{
					queue[i].score += penal;
				}
			}
			else{ // !cMatrix::bSensitive
				queue[i].score += penal;
			}
			if(queue[i].score >= dMaxPenalty){
				lbits &= ~(1L << i);
			}
		}
		for(; i>0; i--){
			if(queue.back().score < dMaxPenalty) break;
			queue.pop_back();
		}
	}
}

bool cAdapter::align(char * read, size_t rLen, uchar * qual, size_t qLen, cElementSet &result, int bc, bool bBestAlign)
{
	bool bDetermined = false;
	ELEMENT elem;
	double dMaxPenalty = cMatrix::dPenaltyPerErr * len + 0.001;
	int iMaxIndel = ceil(cMatrix::dEpsilonIndel * len);
	int minK = bBestAlign ? ((cMatrix::iMinOverlap >= (int)(len - iMaxIndel + 1)) ? (int)(len - iMaxIndel + 1) : cMatrix::iMinOverlap) : 1;
	double dMu = (bc >= 0) ? cMatrix::dMu : MIN_PENALTY;

	deque<ELEMENT> queue;
	ELEMENT element;
	double score;
	uint64 legalBits = 0;
	int i, j, jj;
	element.idx.bc = bc + 1;
	if(trimMode & TRIM_HEAD){
		for(i=1; i<=int(len)-minK; i++){
			element.idx.pos = -i;
			element.score = cMatrix::dPenaltyPerErr * i;
			element.nIndel = 0;
			queue.push_back(element);
			legalBits = (legalBits << 1) | 1;
		}
	}
	else{
		for(i=1,score=cMatrix::dDelta; i<int(len); i++,score+=cMatrix::dDelta){
			if(i > iMaxIndel) break;
			element.idx.pos = -i;
			element.score = score;
			element.nIndel = i;
			queue.push_back(element);
			legalBits = (legalBits << 1) | 1;
		}
	}
	element.nIndel = 0;
	uint64 mbits, xbits, unbits, dnbits, d0bits;
	unbits = dnbits = 0L;
	double penal;
	for(j=0; j<int(rLen); j++){
		jj = j;
		mbits = matchBits[codeMap[uchar(read[jj])]];
		penal = ((qLen > 0) ? cMatrix::penalty[qual[jj]] : dMu);

		element.idx.pos = j;
		element.score = ((mbits & 0x01) == 0) ? penal : 0;
		queue.push_front(element);

		xbits = mbits | unbits;
		dnbits <<= 1;
		unbits <<= 1;
		d0bits = ((dnbits + (xbits & dnbits)) ^ dnbits) | xbits;
		legalBits = (legalBits << 1) | 1;

		UPDATE_COLUMN(queue, d0bits, legalBits, unbits, dnbits, penal, dMaxPenalty, iMaxIndel);

		dnbits &= d0bits;
		unbits &= d0bits;

		if(queue.size() == len){
			if(bBestAlign){
				if(trimMode == TRIM_HEAD){
					i = (queue.back().idx.pos < 0) ? (len + queue.back().idx.pos) : len;
					if( !bDetermined || (i * cMatrix::dMu - queue.back().score) > elem.score * (i+1) ){
						elem = queue.back();
						dMaxPenalty = elem.score;
						elem.score = (i * cMatrix::dMu - elem.score) / (i+1); // normalization
						elem.idx.pos = rLen - 1 - j;
					}
				}
				else{
					elem = queue.back();
					dMaxPenalty = elem.score + ((trimMode == TRIM_TAIL) ? EPSILON : 0);
					elem.score = (len * cMatrix::dMu - elem.score) / (len + 1); // normalization
				}
				bDetermined = true;
				if(dMaxPenalty == 0) break;
			}
			else{
				elem = queue.back();
				elem.score = len * cMatrix::dMu - elem.score; // normalization
				result.insert(elem);
			}
			queue.pop_back();
		}
	}
	if(dMaxPenalty > 0){ // not the case of "perfect match for single-end reads trimming"
		if(bBestAlign){
			if(trimMode & TRIM_TAIL){
				dMaxPenalty = (cMatrix::dPenaltyPerErr * queue.size() + 0.001);
				for(i=queue.size(); i>=minK; i--, dMaxPenalty-=cMatrix::dPenaltyPerErr){
					if(dMaxPenalty <= 0) break;
					if(queue.back().score < dMaxPenalty){
						if(!bDetermined || ((i * cMatrix::dMu - queue.back().score) > elem.score * (i+1)) ){
							elem = queue.back();
							dMaxPenalty = elem.score;
							elem.score = (i * cMatrix::dMu - elem.score) / (i+1); // normalization
							bDetermined = true;
						}
					}
					queue.pop_back();
				}
			}
			else{
				dMaxPenalty -= (len - queue.size()) * cMatrix::dDelta;
				iMaxIndel -= (len - queue.size());
				for(i=queue.size(); i>=minK; i--, dMaxPenalty-=cMatrix::dDelta, iMaxIndel--){
					if( (dMaxPenalty <= 0) || (iMaxIndel < 0) ) break;
					if( (queue.back().score < dMaxPenalty) && (queue.back().nIndel <= iMaxIndel) ){
						if(!bDetermined || ((i * cMatrix::dMu - queue.back().score) > elem.score * (i+1)) ){
							elem = queue.back();
							dMaxPenalty = elem.score;
							elem.score = (i * cMatrix::dMu - elem.score) / (i+1); // normalization
							elem.idx.pos = -(len - i);
							bDetermined = true;
						}
					}
					queue.pop_back();
				}
			}
		}
		else{
			dMaxPenalty = cMatrix::dPenaltyPerErr * queue.size() + 0.001;
			for(i=queue.size(); i>=minK; i--, dMaxPenalty-=cMatrix::dPenaltyPerErr){
				if(queue.back().score < dMaxPenalty){
					elem = queue.back();
					elem.score = i * cMatrix::dMu - elem.score; // normalization
					result.insert(elem);
				}
				queue.pop_back();
			}
		}
	}
	if(bDetermined){
		result.clear();
		result.insert(elem);
	}

	return bDetermined;
}

void cAdapter::initBarcode(int iCut)
{
	int i, k;
	if(iCut > (int)len) iCut = (int)len;
	for(k=0,i=0; i<iCut; i++){
		masked[i] = (sequence[i] == 'N');
		if(masked[i]) continue;
		barcode[k++] = sequence[i];
	}
	barcode[k] = '\0';
	for(k=0; i<(int)len; i++){
		if(sequence[i] == '\0') break;
		primer[k++] = sequence[i];
	}
	primer[k] = '\0';
}

deque<cAdapter> cMatrix::firstAdapters;
deque<cAdapter> cMatrix::secondAdapters;
deque<cAdapter> cMatrix::junctionAdapters;
vector<int> cMatrix::junctionLengths;
vector< vector<int> > cMatrix::indices;
int cMatrix::iIdxCnt = 0;
vector<bool *> cMatrix::fw_masked;
vector<bool *> cMatrix::rv_masked;
vector<string> cMatrix::fw_barcodes;
vector<string> cMatrix::rv_barcodes;
vector<string> cMatrix::fw_primers;
vector<string> cMatrix::rv_primers;
vector<int> cMatrix::rowBc;
vector<int> cMatrix::colBc;
bool cMatrix::bShareAdapter = false;
double cMatrix::dEpsilon = 0.15;
double cMatrix::dEpsilonIndel = 0.03;
double cMatrix::dPenaltyPerErr = cMatrix::dEpsilon * MEAN_PENALTY;
double cMatrix::dDelta = MAX_PENALTY;
double cMatrix::dMu = MEAN_PENALTY;
double cMatrix::penalty[256];
bool cMatrix::bSensitive = false;
int cMatrix::iMinOverlap = 3;

///////////////////////////////////////
cMatrix::cMatrix()
{
}

cMatrix::~cMatrix()
{
}

bool cMatrix::CalcRevCompScore(char * seq, char * seq2, int len, uchar * qual, uchar * qual2, size_t qLen, double &score)
{
	double dMaxPenalty = dPenaltyPerErr * len;
	double penal;
	CODE code, code2;
	if(len <= 0){
		score = (qLen > 0) ? (dMu * qLen / 2) : 0.0;
		// prefer to detect empty reads even if an error ratio of 0.5 is specified
		return true;
	}
	score = 0.0;
	for(int i=0; i<len; i++){
		code = codeMap[uchar(seq[i])];
		code2 = complement[codeMap[uchar(seq2[len-1-i])]];
		penal = scoring[code][code2];
		if(penal > 0.0){
			if(qLen > 0){
				if(cMatrix::penalty[qual[i]] <= cMatrix::penalty[qual2[len-1-i]]){
					penal *= cMatrix::penalty[qual[i]];
				}
				else{
					penal *= cMatrix::penalty[qual2[len-1-i]];
				}
			}
			else{
				penal *= dMu;
			}
			score += penal;
			if(score > dMaxPenalty){
				return false;
			}
		}
	}
	score = len * dMu - score; // normalization
	return true;
}

string cMatrix::GetRevComp(char * seq, int len)
{
	char sequence[MAX_ADAPTER_LEN+1];
	if(len > MAX_ADAPTER_LEN){
		seq += (len - MAX_ADAPTER_LEN);
		len = MAX_ADAPTER_LEN;
	}
    for(int i=0; i<int(len); i++) 
        sequence[i] = character[complement[codeMap[uchar(seq[len-1-i])]]];
    sequence[len] = '\0';
	return string(sequence);
}

//// public functions
void cMatrix::InitParameters(enum TRIM_MODE trimMode, double dEpsilon, double dEpsilonIndel, int baseQual, bool bShareAdapter)
{
	cMatrix::dDelta = (trimMode & TRIM_AP) ? MEAN_PENALTY : MAX_PENALTY;
	cMatrix::dEpsilon = dEpsilon;
	cMatrix::dEpsilonIndel = dEpsilonIndel;
	cMatrix::dPenaltyPerErr = dEpsilon * MEAN_PENALTY;
	cMatrix::bSensitive = (dEpsilonIndel > 0);
	// pre-calcualte the penalties corresponding to quality values
	int chr;
	for(chr=0; chr<=baseQual; chr++){
		cMatrix::penalty[chr] = MIN_PENALTY;
	}
	int i;
	for(i=1; i<40; i++,chr++){
		cMatrix::penalty[chr] = MIN_PENALTY + i / 10.0;
	}
	for(; chr<256; chr++){
		cMatrix::penalty[chr] = MAX_PENALTY;
	}
	cMatrix::bShareAdapter = bShareAdapter;
	cMatrix::firstAdapters.clear();
	cMatrix::secondAdapters.clear();
	cMatrix::junctionAdapters.clear();
}

void cMatrix::AddAdapter(deque<cAdapter> & adapters, char * vector, size_t len, TRIM_MODE trimMode)
{
	cAdapter adapter;
	adapter.Init(vector, len, trimMode);
	adapters.push_back(adapter);
}

void cMatrix::CalculateJunctionLengths()
{
	deque<cAdapter>::iterator it_adapter;
	junctionLengths.push_back(0);
	for(it_adapter=junctionAdapters.begin(); it_adapter!=junctionAdapters.end(); it_adapter++){
		junctionLengths.push_back((*it_adapter).len);
	}
}

void cMatrix::CalculateIndices(vector< vector<bool> > &bMatrix, int nRow, int nCol)
{
	int i, j;
	rowBc.clear();
	colBc.clear();
	indices.resize(nRow, vector<int>(nCol, -1));
	iIdxCnt = 0;
	for(i=0; i<nRow; i++){
		for(j=0; j<nCol; j++){
			if(bMatrix[i][j]){
				rowBc.push_back(i);
				colBc.push_back(j);
				indices[i][j] = iIdxCnt++;
			}
		}
	}
}

void cMatrix::InitBarcodes(deque<cAdapter> & fw_adapters, int iCutF, deque<cAdapter> & rv_adapters, int iCutR)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	fw_masked.clear(); fw_masked.push_back(NULL);
	fw_barcodes.clear(); fw_barcodes.push_back(string(""));
	fw_primers.clear(); fw_primers.push_back(string(""));
	for(it_adapter=fw_adapters.begin(); it_adapter!=fw_adapters.end(); it_adapter++){
		pAdapter = &(*it_adapter);
		pAdapter->initBarcode(iCutF);
		fw_masked.push_back(pAdapter->getMasked());
		fw_barcodes.push_back(string(pAdapter->getBarcode()));
		fw_primers.push_back(string(pAdapter->getPrimer()));
	}
	rv_masked.clear(); rv_masked.push_back(NULL);
	rv_barcodes.clear(); rv_barcodes.push_back(string(""));
	rv_primers.clear(); rv_primers.push_back(string(""));
	for(it_adapter=rv_adapters.begin(); it_adapter!=rv_adapters.end(); it_adapter++){
		pAdapter = &(*it_adapter);
		pAdapter->initBarcode(iCutR);
		rv_masked.push_back(pAdapter->getMasked());
		rv_barcodes.push_back(string(pAdapter->getBarcode()));
		rv_primers.push_back(string(pAdapter->getPrimer()));
	}
}

bool cMatrix::isBlurry(char * seq, size_t len)
{
	size_t u;
	int iMaxBlurry = ceil(cMatrix::dEpsilon * len);
	int iBlurry = 0;
	for(u=0; u<len; u++){
		if(blurry[int(seq[u])]){
			if(++iBlurry > iMaxBlurry){
				return true;
			}
		}
	}
	return false;
}

bool cMatrix::checkQualities(uchar * quals, size_t len, int minQual)
{
	size_t u;
	if(len == 0) return true;
	int total = 0;
	for(u=0; u<len; u++){
		total += quals[u];
	}   
	return (double(total) / len) >= minQual;
}

int cMatrix::trimByQuality(uchar * quals, size_t len, int minQual)
{
	int i;
	for(i=(int)len-1; i>=0; i--){
		if(quals[i] >= minQual)
			break;
	}
	return (i+1);
}

INDEX cMatrix::findAdapter(char * read, size_t rLen, uchar * qual, size_t qLen)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	double maxScore = -1;
	INDEX index;
	index.pos = int(rLen);
	index.bc = 0;
	int i;
	for(i=0,it_adapter=firstAdapters.begin(); it_adapter!=firstAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		if(pAdapter->align(read, rLen, qual, qLen, result, i)){
			if(result.begin()->score > maxScore){
				index = result.begin()->idx;
				maxScore = result.begin()->score;
			}
		}
	}
	return index;
}

INDEX cMatrix::findAdapter2(char * read, size_t rLen, uchar * qual, size_t qLen)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	double maxScore = -1;
	INDEX index;
	index.pos = int(rLen);
	index.bc = 0;
	int i;
	deque<cAdapter> *pAdapters = (bShareAdapter ? &firstAdapters : &secondAdapters);
	for(i=0,it_adapter=pAdapters->begin(); it_adapter!=pAdapters->end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		if(pAdapter->align(read, rLen, qual, qLen, result, i)){
			if(result.begin()->score > maxScore){
				index = result.begin()->idx;
				maxScore = result.begin()->score;
			}
		}
	}
	return index;
}

INDEX cMatrix::findJuncAdapter(char * read, size_t rLen, uchar * qual, size_t qLen)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	double maxScore = -1;
	INDEX index;
	index.pos = int(rLen);
	index.bc = 0;
	int i;
	for(i=0,it_adapter=junctionAdapters.begin(); it_adapter!=junctionAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		if(pAdapter->align(read, rLen, qual, qLen, result, i)){
			if(result.begin()->score > maxScore){
				index = result.begin()->idx;
				maxScore = result.begin()->score;
			}
		}
	}
	return index;
}

bool cMatrix::findAdapterWithPE(char * read, char * read2, size_t rLen, size_t rLen2, uchar * qual, uchar * qual2, size_t qLen, size_t qLen2, INDEX &index, INDEX &index2)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result, result2;
	index.pos = int(rLen);
	index.bc = -1;
	index2.pos = int(rLen2);
	index2.bc = -1;
	int i;
	for(i=0,it_adapter=firstAdapters.begin(); it_adapter!=firstAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		pAdapter->align(read, rLen, qual, qLen, result, i, false);
	}
	deque<cAdapter> *pAdapters = (bShareAdapter ? &firstAdapters : &secondAdapters);
	for(i=0,it_adapter=pAdapters->begin(); it_adapter!=pAdapters->end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		pAdapter->align(read2, rLen2, qual2, qLen2, result2, i, false);
	}
	if(result.empty() && result2.empty()){
		return false;
	}
	size_t minQLen = (qLen <= qLen2) ? qLen : qLen2;
	double maxScore = -1;
	double score;
	cElementSet::iterator it_element, it_element2, it_ele, it_ele2;
	bool bRevComplement;
	int pos, pos2, cpos, iStart, iStart2;
	int apos, apos2;
	int bc, bc2;
	bc = bc2 = 0;
	apos = int(rLen);
	apos2 = int(rLen2);
	it_element = result.begin();
	it_element2 = result2.begin();
	while( true ){
		pos = (it_element == result.end()) ? INT_MAX : it_element->idx.pos;
		pos2 = (it_element2 == result2.end()) ? INT_MAX : it_element2->idx.pos;
		if( (pos == INT_MAX) && (pos2 == INT_MAX) )
			break;
		if(pos <= pos2){ // partial ordering: pos < rLen && pos2 < rLen2 iff. pos2 != INT_MAX
			cpos = pos;
			iStart = (pos <= int(rLen2)) ? 0 : (pos - int(rLen2));
			iStart2 = 0;
		}
		else{ // pos > pos2
			cpos = pos2;
			iStart = 0;
			iStart2 = (pos2 <= int(rLen)) ? 0 : (pos2 - int(rLen));
		}
		bRevComplement = CalcRevCompScore(read + iStart, read2 + iStart2, cpos, qual + iStart, qual2 + iStart2, minQLen, score);
		if(pos < pos2){
			if(bRevComplement){
				do{
					if( (indices[it_element->idx.bc][0] >= 0) && (score + it_element->score > maxScore) ){
						maxScore = score + it_element->score;
						bc = it_element->idx.bc;
						bc2 = 0;
						apos = cpos;
						apos2 = (cpos <= int(rLen2)) ? cpos : int(rLen2);
					}
					it_element++;
				}while( (it_element != result.end()) && (it_element->idx.pos == cpos) );
			}
			else{
				do{
					it_element++;
				}while( (it_element != result.end()) && (it_element->idx.pos == cpos) );
			}
		}
		else if(pos > pos2){
			if(bRevComplement){
				do{
					if( (indices[0][it_element2->idx.bc] >= 0) && (score + it_element2->score > maxScore) ){
						maxScore = score + it_element2->score;
						bc = 0;
						bc2 = it_element2->idx.bc;
						apos = (cpos <= int(rLen)) ? cpos : int(rLen);
						apos2 = cpos;
					}
					it_element2++;
				}while( (it_element2 != result2.end()) && (it_element2->idx.pos == cpos) );
			}
			else{
				do{
					it_element2++;
				}while( (it_element2 != result2.end()) && (it_element2->idx.pos == cpos) );
			}
		}
		else{ // ==
			if(bRevComplement){
				for(it_ele = it_element; (it_ele != result.end()) && (it_ele->idx.pos == cpos); it_ele++){
					for(it_ele2 = it_element2; (it_ele2 != result2.end()) && (it_ele2->idx.pos == cpos); it_ele2++){
						if(indices[it_ele->idx.bc][it_ele2->idx.bc] < 0)
							continue;
						if(score + it_ele->score + it_ele2->score <= maxScore)
							continue;
						maxScore = score + it_ele->score + it_ele2->score;
						bc = it_ele->idx.bc;
						bc2 = it_ele2->idx.bc;
						apos = cpos;
						apos2 = cpos;
					}
				}
			}
			do{
				it_element++;
			}while( (it_element != result.end()) && (it_element->idx.pos == cpos) );
			do{
				it_element2++;
			}while( (it_element2 != result2.end()) && (it_element2->idx.pos == cpos) );
		}
	}
	index.bc = index2.bc = indices[bc][bc2];
	if(index.bc < 0){
		return false;
	}
	if( (apos <= 0) || (apos2 <= 0) ){
		index.pos = index2.pos = 0;
	}
	else{
		index.pos = apos;
		index2.pos = apos2;
	}
	return true;
}

// return value
//  0: forward-reverse
//  1: reverse-forward
// -1: no match
int cMatrix::findAdaptersInARead(char * read, size_t rLen, uchar * qual, size_t qLen, INDEX &index)
{
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	index.pos = int(rLen);
	index.bc = -1;
	int i;
	size_t nLen;
	double maxScore = -1;
	int flag = 0;
	deque<ELEMENT>::iterator it_element, it_element2;
	for(i=0,it_adapter=firstAdapters.begin(); it_adapter!=firstAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		nLen = (pAdapter->len < rLen ? pAdapter->len : rLen);
		if(pAdapter->align(read, nLen, qual, qLen, result, i)){
			if(result.begin()->score > maxScore){
				index = result.begin()->idx;
				index.bc = indices[index.bc][0];
				maxScore = result.begin()->score;
			}
		}
	}
	if(!bShareAdapter){
		for(i=0,it_adapter=secondAdapters.begin(); it_adapter!=secondAdapters.end(); it_adapter++,i++){
			pAdapter = &(*it_adapter);
			nLen = (pAdapter->len < rLen ? pAdapter->len : rLen);
			if(pAdapter->align(read, nLen, qual, qLen, result, i)){
				if(result.begin()->score > maxScore){
					index = result.begin()->idx;
					index.bc = indices[0][index.bc];
					maxScore = result.begin()->score;
					flag = 1;
				}
			}
		}
	}
	return (index.bc < 0) ? -1 : flag;
}

// return value
//  0: forward-reverse
//  1: reverse-forward
// -1: no match
int cMatrix::findAdaptersBidirectionally(char * read, size_t rLen, uchar * qual, size_t qLen,
char * read2, size_t rLen2, uchar * qual2, size_t qLen2, INDEX &index, INDEX &index2)
{
	int bc = -1;
	deque<cAdapter>::iterator it_adapter;
	cAdapter * pAdapter;
	cElementSet result;
	deque<ELEMENT> result1, result2, result3, result4;
	index.pos = index2.pos = int(rLen);
	index.bc = index2.bc = 0;
	int i;
	size_t nLen;
	deque<ELEMENT>::iterator it_element, it_element2;
	for(i=0,it_adapter=firstAdapters.begin(); it_adapter!=firstAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		nLen = (pAdapter->len < rLen ? pAdapter->len : rLen);
		if(pAdapter->align(read, nLen, qual, qLen, result, i)){
			result1.push_back(*result.begin());
		}
		nLen = (pAdapter->len < rLen2 ? pAdapter->len : rLen2);
		if(pAdapter->align(read2, nLen, qual2, qLen2, result, i)){
			result3.push_back(*result.begin());
		}
	}
	ELEMENT eElement;
	memset((void *)&eElement, 0, sizeof(ELEMENT));
	double maxScore = -1;
	if(bShareAdapter){
		if(result1.empty() ^ result3.empty()){
			if(result1.empty()){
				result1.push_back(eElement);
			}
			else{ // result3.empty()
				result3.push_back(eElement);
			}
		}
		for(it_element=result1.begin(); it_element!=result1.end(); it_element++){
			for(it_element2=result3.begin(); it_element2!=result3.end(); it_element2++){
				if(indices[it_element->idx.bc][it_element2->idx.bc] < 0)
					continue;
				if(it_element->score + it_element2->score > maxScore){
					maxScore = it_element->score + it_element2->score;
					index = it_element->idx;
					index2 = it_element2->idx;
					bc = indices[it_element->idx.bc][it_element2->idx.bc];
				}
			}
		}
		index.bc = bc;
		return (bc < 0) ? -1 : 0;
	}
	for(i=0,it_adapter=secondAdapters.begin(); it_adapter!=secondAdapters.end(); it_adapter++,i++){
		pAdapter = &(*it_adapter);
		nLen = (pAdapter->len < rLen2 ? pAdapter->len : rLen2);
		if(pAdapter->align(read2, nLen, qual2, qLen2, result, i)){
			result2.push_back(*result.begin());
		}
		nLen = (pAdapter->len < rLen ? pAdapter->len : rLen);
		if(pAdapter->align(read, nLen, qual, qLen, result, i)){
			result4.push_back(*result.begin());
		}
	}
	if(result1.empty() ^ result2.empty()){
		if(result1.empty()){
			result1.push_back(eElement);
		}
		else{ // result2.empty()
			result2.push_back(eElement);
		}
	}
	for(it_element=result1.begin(); it_element!=result1.end(); it_element++){
		for(it_element2=result2.begin(); it_element2!=result2.end(); it_element2++){
			if(indices[it_element->idx.bc][it_element2->idx.bc] < 0)
				continue;
			if(it_element->score + it_element2->score > maxScore){
				maxScore = it_element->score + it_element2->score;
				index = it_element->idx;
				index2 = it_element2->idx;
				bc = indices[it_element->idx.bc][it_element2->idx.bc];
			}
		}
	}
	bool bReverse = false;
	if(result3.empty() ^ result4.empty()){
		if(result3.empty()){
			result3.push_back(eElement);
		}
		else{ // result4.empty()
			result4.push_back(eElement);
		}
	}
	for(it_element=result3.begin(); it_element!=result3.end(); it_element++){
		for(it_element2=result4.begin(); it_element2!=result4.end(); it_element2++){
			if(indices[it_element->idx.bc][it_element2->idx.bc] < 0)
				continue;
			if(it_element->score + it_element2->score > maxScore){
				maxScore = it_element->score + it_element2->score;
				index2 = it_element->idx;
				index = it_element2->idx;
				bc = indices[it_element->idx.bc][it_element2->idx.bc];
				bReverse = true;
			}
		}
	}
	index.bc = bc;
	return (bc < 0) ? -1 : bReverse;
}

bool cMatrix::PrepareBarcode(char * barcodeSeq, int bcIdx, char * seq, int len, char * seq2, int len2, char * barcodeQual, char * qual, char * qual2)
{
	assert(bcIdx >= 0);
	bool *mask = fw_masked[rowBc[bcIdx]];
	bool *mask2 = rv_masked[colBc[bcIdx]];
	if( (mask == NULL) || (mask2 == NULL) ){
		return false;
	}
	int n = 0;
	int i;
	for(i=0; i<len; i++){
		if(!mask[i]){
			barcodeSeq[n] = seq[i];
			barcodeQual[n] = qual[i];
			n++;
		}
	}
	for(i=0; i<len2; i++){
		if(!mask2[i]){
			barcodeSeq[n] = seq2[i];
			barcodeQual[n] = qual2[i];
			n++;
		}
	}
	return true;
}

bool cMatrix::PrepareBarcode(char * barcodeSeq, int bcIdx, char * seq, int len, char * seq2, int len2)
{
	assert(bcIdx >= 0);
	bool *mask = fw_masked[rowBc[bcIdx]];
	bool *mask2 = rv_masked[colBc[bcIdx]];
	if( (mask == NULL) || (mask2 == NULL) ){
		return false;
	}
	int n = 0;
	int i;
	for(i=0; i<len; i++){
		if(!mask[i]){
			barcodeSeq[n++] = seq[i];
		}
	}
	for(i=0; i<len2; i++){
		if(!mask2[i]){
			barcodeSeq[n++] = seq2[i];
		}
	}
	return true;
}

INDEX cMatrix::mergePE(char * read, char * read2, size_t rLen, uchar * qual, uchar * qual2, size_t qLen, size_t startPos, size_t jLen)
{
	INDEX index;
	cElementSet result;
	cAdapter adapter;
	double score;
	int pos, clen;
	bool bRevComplement = false;
	size_t endPos, eLen;
	char chr;
	uchar uchr;
	int i;
	index.pos = int(rLen);
	index.bc = 0;
	adapter.Init2(read2, rLen);
	if(adapter.align(read, rLen, NULL, 0, result, -1)){
		pos = result.begin()->idx.pos;
		if(pos >= 0){
			clen = rLen - pos;
			bRevComplement = CalcRevCompScore(read + pos, read2 + pos, clen, qual + pos, qual2 + pos, qLen, score);
		}
	}
	if(bRevComplement){ // overlap detected
		if(qLen > 0)
			combinePairSeqs(read+pos, read2+pos, clen, clen, qual+pos, qual2+pos, qLen, qLen);
		endPos = startPos + jLen;
		if(endPos + clen >= rLen){ // junction adapter locates in overlapping region
			index.pos = rLen - (endPos + clen - rLen);
		}
		else{
			eLen = rLen - (endPos + clen);
			index.pos += eLen;
			for(i=eLen/2; i>=0; i--){ //reverse
				chr = read2[endPos + i];
				read2[endPos + i] = read2[endPos + eLen - 1 - i];
				read2[endPos + eLen - 1 - i] = chr;
			}
			for(i=0; i<(int)eLen; i++){
				read2[startPos + i] = character[complement[codeMap[uchar(read2[endPos + i])]]];
			}
			if(qLen > 0){
				for(i=eLen/2; i>=0; i--){ //reverse
					uchr = qual2[endPos + i];
					qual2[endPos + i] = qual2[endPos + eLen - 1 -i];
					qual2[endPos + eLen - 1 - i] = uchr;
				}
				for(i=0; i<(int)eLen; i++){
					qual2[startPos + i] = qual2[endPos + i];
				}
			}
		}
	}
	else{
		index.pos -= cMatrix::iMinOverlap;
	}
	return index;
}

bool cMatrix::combinePairSeqs(char * read, char * read2, int len, int len2, uchar * qual, uchar * qual2, int qLen, int qLen2)
{
	CODE code, code2;
	if(len != len2){
		int offset;
		if(len > len2){
			offset = len - len2;
			read += offset;
			len -= offset;
			qual += offset;
			qLen -= offset;
		}
		else{
			offset = len2 - len;
			read2 += offset;
			qual2 += offset;
			qLen2 -= offset;
		}
	}
	int minQLen = (qLen < qLen2) ? qLen : qLen2;
	if(minQLen < len){
		return false;
	}
	for(int i=0; i<len; i++){
		code = codeMap[uchar(read[i])];
		code2 = complement[codeMap[uchar(read2[len-1-i])]];
		if(qual2[len-1-i] > qual[i]){
			qual[i] = qual2[len-1-i];
			if(code != code2){
				read[i] = character[code2];
			}
		}
	}
	return true;
}
