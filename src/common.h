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
#ifndef _COMMON_H
#define _COMMON_H

typedef unsigned long long uint64;
typedef unsigned long uint32;
typedef unsigned short uint16;
typedef long long int64;
typedef long int32;
typedef short int16;

#ifndef uint
typedef unsigned int uint;
#endif

typedef unsigned char uchar;

const int MAX_PATH = 255;
const int MAX_ADAPTER_LEN = 64;
const int MAX_ADAPTER_CNT = 96;

typedef struct tag_INDEX{
	int pos;
	int pos2;
	int bc;
}INDEX;

enum TRIM_MODE{
	TRIM_DEFAULT = 0,
	TRIM_HEAD = 1,
	TRIM_TAIL = 2,
	TRIM_ANY = 3,
	TRIM_PE = 4,
	TRIM_PE_HEAD = 5,
	TRIM_PE_TAIL = 6,
	TRIM_PE_ANY = 7,
	TRIM_MP = 8,
	TRIM_MP_HEAD = 9,
	TRIM_MP_TAIL = 10,
	TRIM_MP_ANY = 11,
	TRIM_AP = 16,
	TRIM_AP_HEAD = 17,
	TRIM_MODE_CNT = 14
};

#endif // _COMMON_H
