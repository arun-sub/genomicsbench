/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
*****************************************************************************************/

#include "omp.h" 
#include "bandedSWA.h"
#ifdef VTUNE_ANALYSIS
#include <ittnotify.h> 
#endif

#if defined(__clang__) || defined(__GNUC__)
#define __mmask8 uint8_t
#define __mmask16 uint16_t
#define __mmask32 uint32_t
#endif

// ------------------------------------------------------------------------------------
// MACROs for vector code
extern uint64_t prof[10][112];
#define AMBIG 4
#define DUMMY1 99
#define DUMMY2 100

//-----------------------------------------------------------------------------------
// constructor
BandedPairWiseSW::BandedPairWiseSW(const int o_del, const int e_del, const int o_ins,
                                   const int e_ins, const int zdrop,
                                   const int end_bonus, const int8_t *mat_,
                                   const int8_t w_match, const int8_t w_mismatch, int numThreads)
{
    mat = mat_;
    this->m = 5;
    this->end_bonus = end_bonus;
    this->zdrop = zdrop;
    this->o_del = o_del;
    this->o_ins = o_ins;
    this->e_del = e_del;
    this->e_ins = e_ins;
    
    this->w_match    = w_match;
    this->w_mismatch = w_mismatch*-1;
    this->w_open     = o_del;  // redundant, used in vector code.
    this->w_extend   = e_del;  // redundant, used in vector code.
    this->w_ambig    = DEFAULT_AMBIG;
    this->swTicks = 0;
    this->SW_cells = 0;
    setupTicks = 0;
    sort1Ticks = 0;
    swTicks = 0;
    swComputationTicks = 0;
    swBandAdjustmentTicks = 0;
    sort2Ticks = 0;
#ifdef PROFILE
    numCellsComputed = 0;
#endif
    this->F8_ = this->H8_  = this->H8__ = NULL;
    this->F16_ = this->H16_  = this->H16__ = NULL;
    
    F8_ = H8_ = H8__ = NULL;
    F8_ = (int8_t *)_mm_malloc(MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(int8_t), 64);
    H8_ = (int8_t *)_mm_malloc(MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(int8_t), 64);
    H8__ = (int8_t *)_mm_malloc(MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(int8_t), 64);

    F16_ = H16_ = H16__ = NULL;
    F16_ = (int16_t *)_mm_malloc(MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    H16_ = (int16_t *)_mm_malloc(MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);
    H16__ = (int16_t *)_mm_malloc(MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(int16_t), 64);

    if (F8_ == NULL || H8_ == NULL || H8__ == NULL) {
        printf("BSW8 Memory not alloacted!!!\n"); exit(EXIT_FAILURE);
    }       
    if (F16_ == NULL || H16_ == NULL || H16__ == NULL) {
        printf("BSW16 Memory not alloacted!!!\n"); exit(EXIT_FAILURE);
    }       
}

// destructor 
BandedPairWiseSW::~BandedPairWiseSW() {
    _mm_free(F8_); _mm_free(H8_); _mm_free(H8__);
    _mm_free(F16_);_mm_free(H16_); _mm_free(H16__);
}

int64_t BandedPairWiseSW::getTicks()
{
    //printf("oneCount = %ld, totalCount = %ld\n", oneCount, totalCount);
    int64_t totalTicks = sort1Ticks + setupTicks + swTicks + sort2Ticks;
    #ifdef PROFILE
        printf("cost breakup: Sort1:%ld, Setup+SOA:%ld, SWComputation:%ld, SWBandAdjustment:%ld, Sort2:%ld, Total:%ld\n",
            sort1Ticks, setupTicks + (swTicks - swComputationTicks), prof[DP1][0], prof[DP2][0] + prof[DP3][0], sort2Ticks,
            totalTicks);
    #else 
        // printf("cost breakup: Sort1:%ld, Setup:%ld, SWComputation:%ld, SWComputationOnly:%ld, Sort2:%ld, Total:%ld\n",
        //        sort1Ticks, setupTicks, swTicks, swComputationTicks, sort2Ticks,
        //        totalTicks);
    #endif
    return totalTicks;
}

// ------------------------------------------------------------------------------------
// Banded SWA - scalar code
// ------------------------------------------------------------------------------------

int BandedPairWiseSW::scalarBandedSWA(int qlen, const uint8_t *query,
                                      int tlen, const uint8_t *target,
                                      int32_t w, int h0, int *_qle, int *_tle,
                                      int *_gtle, int *_gscore,
                                      int *_max_off) {
    
    // uint64_t sw_cells = 0;
    eh_t *eh; // score array
    int8_t *qp; // query profile
    int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
    
    // assert(h0 > 0); //check !!!
    
    // allocate memory
    qp = (int8_t *) malloc(qlen * m);
    assert(qp != NULL);
    eh = (eh_t *) calloc(qlen + 1, 8);
    assert(eh != NULL);

    // generate the query profile
    for (k = i = 0; k < m; ++k) {
        const int8_t *p = &mat[k * m];
        //for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]-48];  //sub 48
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
    for (i = 0; (i < tlen); ++i) {
        int t, f = 0, h1, m = 0, mj = -1;
        //int8_t *q = &qp[(target[i]-48) * qlen];   // sub 48
        int8_t *q = &qp[(target[i]) * qlen];
        // apply the band and the constraint (if provided)
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > qlen) end = qlen;
        // compute the first column
        if (beg == 0) {
            h1 = h0 - (o_del + e_del * (i + 1));
            if (h1 < 0) h1 = 0;
        } else h1 = 0;
        for (j = beg; (j < end); ++j) {
            // At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
            // Similar to SSE2-SW, cells are computed in the following order:
            //   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
            //   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
            //   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
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
            // SW_cells++;
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
        for (j = beg; (j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
        beg = j;
        for (j = end; (j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
        end = j + 2 < qlen? j + 2 : qlen;
        //beg = 0; end = qlen; // uncomment this line for debugging
    }
    free(eh); free(qp);
    if (_qle) *_qle = max_j + 1;
    if (_tle) *_tle = max_i + 1;
    if (_gtle) *_gtle = max_ie + 1;
    if (_gscore) *_gscore = gscore;
    if (_max_off) *_max_off = max_off;
    
#if MAXI
    fprintf(stderr, "%d (%d %d) %d %d %d\n", max, max_i+1, max_j+1, gscore, max_off, max_ie+1);
#endif

    // return sw_cells;
    return max;
}

// -------------------------------------------------------------
// Banded SWA, wrapper function
//-------------------------------------------------------------
void BandedPairWiseSW::scalarBandedSWAWrapper(SeqPair *seqPairArray,
                                              uint8_t *seqBufRef,
                                              uint8_t *seqBufQer,
                                              int numPairs,
                                              int nthreads,
                                              int32_t w) {

    for (int i=0; i<numPairs; i++)
    {
        SeqPair *p = seqPairArray + i;
        uint8_t *seq1 = seqBufRef + p->idr;
        uint8_t *seq2 = seqBufQer + p->idq;
        
        p->score = scalarBandedSWA(p->len2, seq2, p->len1,
                                   seq1, w, p->h0, &p->qle, &p->tle,
                                   &p->gtle, &p->gscore, &p->max_off);      
    }

}
