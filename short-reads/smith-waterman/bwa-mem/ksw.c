/* The MIT License

   Copyright (c) 2011 by Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <emmintrin.h>
#include "ksw.h"

#define PROFILE 1

int64_t numCellsComputed = 0;

unsigned char seq_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// #define MAX_SEQ 1000
#define MAX_SEQ 20000000

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif

/********************
 *** SW extension ***
 ********************/

typedef struct {
	int32_t h, e;
} eh_t;

int ksw_extend2(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w, int end_bonus, int zdrop, int h0, int *_qle, int *_tle, int *_gtle, int *_gscore, int *_max_off)
{
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	assert(h0 > 0);
	// allocate memory
	qp = malloc(qlen * m);
	eh = calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
	for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
		eh[j].h = eh[j-1].h - e_ins;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del; // TODO: is this necessary?
	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
#ifdef PROFILE
	numCellsComputed = 0;
#endif
	for (i = 0; i < tlen; ++i) {
		int t, f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		// compute the first column
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		for (j = beg; j < end; ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
#ifdef PROFILE
			numCellsComputed++;
#endif
			eh_t *p = &eh[j];
			int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
			h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			mj = m > h? mj : j; // record the position where max score is achieved
			m = m > h? m : h;   // m is stored at eh[mj+1]
			t = M - oe_del;
			t = t > 0? t : 0;
			e -= e_del;
			e = e > t? e : t;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			t = M - oe_ins;
			t = t > 0? t : 0;
			f -= e_ins;
			f = f > t? f : t;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m == 0) break;
		if (m > max) {
			max = m, max_i = i, max_j = mj;
			max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
		} else if (zdrop > 0) {
			if (i - max_i > mj - max_j) {
				if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
			}
		}
		// update beg and end for the next round
		for (j = beg; j < end && eh[j].h == 0 && eh[j].e == 0; ++j);
		beg = j;
		for (j = end; j >= beg && eh[j].h == 0 && eh[j].e == 0; --j);
		end = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	free(eh); free(qp);
#ifdef PROFILE
	printf("Qlen:%d,Rlen:%d,ScalarNumCellsComputed:%lld\n", qlen, tlen, numCellsComputed);
#endif
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}

#define MAX_LINE_LEN 2048

void __parsec_roi_begin() {
	printf("ROI Begin.");
}

void __parsec_roi_end() {
	printf("ROI End.");
}

int main(int argc, char *argv[])
{
	int c, sa = 1, sb = 4, i, j, k;
	int8_t mat[25];
	int gapo = 6, gape = 1;
	uint8_t *read_ar = 0, *ref_ar = 0;

	read_ar = (uint8_t*) malloc(MAX_SEQ * sizeof(uint8_t));
	ref_ar = (uint8_t*) malloc(MAX_SEQ * sizeof(uint8_t));

    FILE* fp_in, *fp_out;

	if (argc < 3) {
		fprintf(stderr, "Usage: ksw <input.txt> <output.txt>\n");
		return 1;
	}

	// initialize scoring matrix
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? sa : -sb;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;

	// open file
	fp_in = fopen(argv[1], "r");
	fp_out = fopen(argv[2], "w");

	if (fp_in == NULL) {
		fprintf(stderr, "Unable to open %s for reading\n", argv[1]);
		return 1;
	}

	if (fp_out == NULL) {
		fprintf(stderr, "Unable to open %s for writing\n", argv[2]);
		return 1;
	}

	char line[MAX_LINE_LEN], *linePtr;
	char* delim = ",";
	int numLinesRead = 0;

	//
	// Read input data
	//
	while (fgets(line, sizeof(line), fp_in) != NULL && numLinesRead < MAX_SEQ) {
		char* query = 0, *target = 0;
		int h0, qlen, tlen, w;
		linePtr = line;
		query = strsep(&linePtr, delim);
		target = strsep(&linePtr, delim);
		if (query == NULL || target == NULL) continue;
		h0 = atoi(strsep(&linePtr, delim));
		w = atoi(strsep(&linePtr, delim));
		qlen = strlen(query);
		tlen = strlen(target);
		read_ar = (uint8_t*) malloc(qlen * sizeof(uint8_t));
		ref_ar = (uint8_t*) malloc(tlen * sizeof(uint8_t));
		for (i = 0; i < qlen; ++i) {
			read_ar[i] = seq_nt4_table[(int)query[i]]; 
		}
		for (i = 0; i < tlen; ++i) {
			ref_ar[i] = seq_nt4_table[(int)target[i]]; 
		}
		int score, qle, tle, gtle, gscore, max_off;
		// __parsec_roi_begin();
		score = ksw_extend2(qlen, read_ar, tlen, ref_ar, 5, mat, gapo, gape, gapo, gape, w, 5, 100, h0, &qle, &tle, &gtle, &gscore, &max_off); 
		// __parsec_roi_end();
		fprintf(fp_out, "%d,%d,%d,%d,%d,%d\n", score, qle, tle, gtle, gscore, max_off);
		numLinesRead += 1;
		free(read_ar);
		free(ref_ar);
	}

	fclose(fp_in);
	fclose(fp_out);
	return 0;
}
