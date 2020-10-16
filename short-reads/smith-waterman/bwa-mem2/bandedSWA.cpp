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
    sort2Ticks = 0;
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
    printf("cost breakup: %ld, %ld, %ld, %ld, %ld\n",
            sort1Ticks, setupTicks, swTicks, sort2Ticks,
            totalTicks);

    return totalTicks;
}

#if ((!__AVX512BW__) & (__AVX2__))

//------------------------------------------------------------------------------
// MACROs
// ------------------------ vec-8 ---------------------------------------------
#define ZSCORE8(i4_256, y4_256)                                         \
    {                                                                   \
        __m256i tmpi = _mm256_sub_epi8(i4_256, x256);                   \
        __m256i tmpj = _mm256_sub_epi8(y4_256, y256);                   \
        cmp = _mm256_cmpgt_epi8(tmpi, tmpj);                            \
        score256 = _mm256_sub_epi8(maxScore256, maxRS1);                \
        __m256i insdel = _mm256_blendv_epi8(e_ins256, e_del256, cmp);   \
        __m256i sub_a256 = _mm256_sub_epi8(tmpi, tmpj);                 \
        __m256i sub_b256 = _mm256_sub_epi8(tmpj, tmpi);                 \
        tmp = _mm256_blendv_epi8(sub_b256, sub_a256, cmp);              \
        tmp = _mm256_sub_epi8(score256, tmp);                           \
        cmp = _mm256_cmpgt_epi8(tmp, zdrop256);                         \
        exit0 = _mm256_blendv_epi8(exit0, zero256, cmp);                \
    }


#define MAIN_CODE8(s1, s2, h00, h11, e11, f11, f21, zero256,  maxScore256, e_ins256, oe_ins256, e_del256, oe_del256, y256, maxRS) \
    {                                                                   \
        __m256i cmp11 = _mm256_cmpeq_epi8(s1, s2);                      \
        __m256i sbt11 = _mm256_blendv_epi8(mismatch256, match256, cmp11); \
        __m256i tmp256 = _mm256_max_epu8(s1, s2);                       \
        /*tmp256 = _mm256_cmpeq_epi8(tmp256, val102);*/                 \
        sbt11 = _mm256_blendv_epi8(sbt11, w_ambig_256, tmp256);         \
        __m256i m11 = _mm256_add_epi8(h00, sbt11);                      \
        cmp11 = _mm256_cmpeq_epi8(h00, zero256);                        \
        m11 = _mm256_blendv_epi8(m11, zero256, cmp11);                  \
        h11 = _mm256_max_epi8(m11, e11);                                \
        h11 = _mm256_max_epi8(h11, f11);                                \
        __m256i temp256 = _mm256_sub_epi8(m11, oe_ins256);              \
        __m256i val256  = _mm256_max_epi8(temp256, zero256);            \
        e11 = _mm256_sub_epi8(e11, e_ins256);                           \
        e11 = _mm256_max_epi8(val256, e11);                             \
        temp256 = _mm256_sub_epi8(m11, oe_del256);                      \
        val256  = _mm256_max_epi8(temp256, zero256);                    \
        f21 = _mm256_sub_epi8(f11, e_del256);                           \
        f21 = _mm256_max_epi8(val256, f21);                             \
    }

// ------------------------ vec 16 --------------------------------------------------
#define _mm256_blendv_epi16(a,b,c)              \
        _mm256_blendv_epi8(a, b, c);            


#define ZSCORE16(i4_256, y4_256)                                            \
    {                                                                   \
        __m256i tmpi = _mm256_sub_epi16(i4_256, x256);                  \
        __m256i tmpj = _mm256_sub_epi16(y4_256, y256);                  \
        cmp = _mm256_cmpgt_epi16(tmpi, tmpj);                           \
        score256 = _mm256_sub_epi16(maxScore256, maxRS1);               \
        __m256i insdel = _mm256_blendv_epi16(e_ins256, e_del256, cmp);  \
        __m256i sub_a256 = _mm256_sub_epi16(tmpi, tmpj);                    \
        __m256i sub_b256 = _mm256_sub_epi16(tmpj, tmpi);                    \
        tmp = _mm256_blendv_epi16(sub_b256, sub_a256, cmp);             \
        tmp = _mm256_sub_epi16(score256, tmp);                          \
        cmp = _mm256_cmpgt_epi16(tmp, zdrop256);                            \
        exit0 = _mm256_blendv_epi16(exit0, zero256, cmp);               \
    }


#define MAIN_CODE16(s1, s2, h00, h11, e11, f11, f21, zero256,  maxScore256, e_ins256, oe_ins256, e_del256, oe_del256, y256, maxRS) \
    {                                                                   \
        __m256i cmp11 = _mm256_cmpeq_epi16(s1, s2);                     \
        __m256i sbt11 = _mm256_blendv_epi16(mismatch256, match256, cmp11); \
        __m256i tmp256 = _mm256_max_epu16(s1, s2);                      \
        sbt11 = _mm256_blendv_epi16(sbt11, w_ambig_256, tmp256);        \
        __m256i m11 = _mm256_add_epi16(h00, sbt11);                     \
        cmp11 = _mm256_cmpeq_epi16(h00, zero256);                       \
        m11 = _mm256_blendv_epi16(m11, zero256, cmp11);                 \
        h11 = _mm256_max_epi16(m11, e11);                               \
        h11 = _mm256_max_epi16(h11, f11);                               \
        __m256i temp256 = _mm256_sub_epi16(m11, oe_ins256);             \
        __m256i val256  = _mm256_max_epi16(temp256, zero256);           \
        e11 = _mm256_sub_epi16(e11, e_ins256);                          \
        e11 = _mm256_max_epi16(val256, e11);                            \
        temp256 = _mm256_sub_epi16(m11, oe_del256);                     \
        val256  = _mm256_max_epi16(temp256, zero256);                   \
        f21 = _mm256_sub_epi16(f11, e_del256);                          \
        f21 = _mm256_max_epi16(val256, f21);                            \
    }

// MACROs section ends
// ------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------
// B-SWA - Vector code
// ------------------------- AVX2 - 8 bit SIMD_LANES ---------------------------
inline void sortPairsLen(SeqPair *pairArray, int32_t count, SeqPair *tempArray,
                         int16_t *hist)
{

    int32_t i;
    __m256i zero256 = _mm256_setzero_si256();
    for(i = 0; i <= MAX_SEQ_LEN16; i+=16)
        _mm256_store_si256((__m256i *)(hist + i), zero256);
    
    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        hist[sp.len1]++;
    }

    int32_t cumulSum = 0;
    for(i = 0; i <= MAX_SEQ_LEN16; i++)
    {
        int32_t cur = hist[i];
        hist[i] = cumulSum;
        // histb[i] = cumulSum;
        cumulSum += cur;
    }

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = hist[sp.len1];
        tempArray[pos] = sp;
        hist[sp.len1]++;
    }

    for(i = 0; i < count; i++) {
        pairArray[i] = tempArray[i];
    }
}

inline void sortPairsId(SeqPair *pairArray, int32_t first, int32_t count,
                        SeqPair *tempArray)
{

    int32_t i;

    for(i = 0; i < count; i++)
    {
        SeqPair sp = pairArray[i];
        int32_t pos = sp.id - first;
        tempArray[pos] = sp;
    }

    for(i = 0; i < count; i++)
        pairArray[i] = tempArray[i];    
}

/******************* Vector code, version 2.0 *************************/
#define PFD 2
void BandedPairWiseSW::getScores8(SeqPair *pairArray,
                                  uint8_t *seqBufRef,
                                  uint8_t *seqBufQer,
                                  int32_t numPairs,
                                  uint16_t numThreads,
                                  int32_t w)
{
    int64_t startTick, endTick;
    
    smithWatermanBatchWrapper8(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);

#if MAXI
    printf("AVX2 Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].tle, pairArray[l].qle,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].gtle);
    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif
    
}

void BandedPairWiseSW::smithWatermanBatchWrapper8(SeqPair *pairArray,
                                                  uint8_t *seqBufRef,
                                                  uint8_t *seqBufQer,
                                                  int32_t numPairs,
                                                  uint16_t numThreads,
                                                  int32_t w)
{
    int64_t st1, st2, st3, st4, st5;
#if RDT
    st1 = __rdtsc();
#endif
    
    uint8_t *seq1SoA = NULL;
    seq1SoA = (uint8_t *)_mm_malloc(MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    
    uint8_t *seq2SoA = NULL;
    seq2SoA = (uint8_t *)_mm_malloc(MAX_SEQ_LEN8 * SIMD_WIDTH8 * numThreads * sizeof(uint8_t), 64);
    
    if (seq1SoA == NULL || seq2SoA == NULL) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }
    
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH8 - 1)/SIMD_WIDTH8 ) * SIMD_WIDTH8;
    // assert(roundNumPairs < BATCH_SIZE * SEEDS_PER_READ);
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

#if RDT
    st2 = __rdtsc();
#endif
    
#if SORT_PAIRS     // disbaled in bwa-mem2 (only used in separate benchmark bsw code)

    // Sort the sequences according to decreasing order of lengths
    SeqPair *tempArray = (SeqPair *)_mm_malloc(SORT_BLOCK_SIZE * numThreads *
                                               sizeof(SeqPair), 64);
    int16_t *hist = (int16_t *)_mm_malloc((MAX_SEQ_LEN8 + 32) * numThreads *
                                          sizeof(int16_t), 64);

#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
        int16_t *myHist = hist + tid * (MAX_SEQ_LEN8 + 32);

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsLen(pairArray + first, last - first, myTempArray, myHist);
        }
    }
    _mm_free(hist);
#endif

#if RDT
    st3 = __rdtsc();
#endif
    
    int eb = end_bonus;
#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid = omp_get_thread_num();
        // uint16_t tid = 0;
        uint8_t *mySeq1SoA = NULL;
        mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;

        uint8_t *mySeq2SoA = NULL;
        mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN8 * SIMD_WIDTH8;
        assert(mySeq1SoA != NULL && mySeq2SoA != NULL);
        
        uint8_t *seq1;
        uint8_t *seq2;
        uint8_t h0[SIMD_WIDTH8]   __attribute__((aligned(64)));
        uint8_t band[SIMD_WIDTH8];      
        uint8_t qlen[SIMD_WIDTH8] __attribute__((aligned(64)));
        int32_t bsize = 0;
        
        int8_t *H1 = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
        int8_t *H2 = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
        
        __m256i zero256   = _mm256_setzero_si256();
        __m256i e_ins256  = _mm256_set1_epi8(e_ins);
        __m256i oe_ins256 = _mm256_set1_epi8(o_ins + e_ins);
        __m256i o_del256  = _mm256_set1_epi8(o_del);
        __m256i e_del256  = _mm256_set1_epi8(e_del);
        __m256i eb_ins256 = _mm256_set1_epi8(eb - o_ins);
        __m256i eb_del256 = _mm256_set1_epi8(eb - o_del);
        
        int8_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;
        
        int nstart = 0, nend = numPairs;

        
#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH8)
        {
            int32_t j, k;
            uint8_t maxLen1 = 0;
            uint8_t maxLen2 = 0;
            bsize = w;
            
            uint64_t tim;
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                h0[j] = sp.h0;
                seq1 = seqBufRef + (int64_t)sp.idr;

                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = (seq1[k] == AMBIG?0xFF:seq1[k]);
                    H2[k * SIMD_WIDTH8 + j] = 0;
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
            
            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH8 + j] = DUMMY1;
                    H2[k * SIMD_WIDTH8 + j] = DUMMY1;
                }
            }

//--------------------
            __m256i h0_256 = _mm256_load_si256((__m256i*) h0);
            _mm256_store_si256((__m256i *) H2, h0_256);
            __m256i tmp256 = _mm256_sub_epi8(h0_256, o_del256);

            for(k = 1; k < maxLen1; k++) {
                tmp256 = _mm256_sub_epi8(tmp256, e_del256);
                __m256i tmp256_ = _mm256_max_epi8(tmp256, zero256);
                _mm256_store_si256((__m256i *)(H2 + k* SIMD_WIDTH8), tmp256_);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH8; j++)
            {               
                SeqPair sp = pairArray[i + j];
                seq2 = seqBufQer + (int64_t)sp.idq;
                
                if (sp.len2 > MAX_SEQ_LEN8) fprintf(stderr, "Error !! : %d %d\n", sp.id, sp.len2);
                assert(sp.len2 < MAX_SEQ_LEN8);
                
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = (seq2[k]==AMBIG?0xFF:seq2[k]);
                    H1[k * SIMD_WIDTH8 + j] = 0;                    
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }

            for(j = 0; j < SIMD_WIDTH8; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH8 + j] = DUMMY2;
                    H1[k * SIMD_WIDTH8 + j] = 0;
                }
            }

//------------------------
            _mm256_store_si256((__m256i *) H1, h0_256);
            __m256i cmp256 = _mm256_cmpgt_epi8(h0_256, oe_ins256);
            tmp256 = _mm256_sub_epi8(h0_256, oe_ins256);
            // _mm256_store_si256((__m256i *) (H1 + SIMD_WIDTH8), tmp256);

            tmp256 = _mm256_blendv_epi8(zero256, tmp256, cmp256);
            _mm256_store_si256((__m256i *) (H1 + SIMD_WIDTH8), tmp256);
            for(k = 2; k < maxLen2; k++)
            {
                // __m256i h1_256 = _mm256_load_si256((__m256i *) (H1 + (k-1) * SIMD_WIDTH8));
                __m256i h1_256 = tmp256;
                tmp256 = _mm256_sub_epi8(h1_256, e_ins256);
                tmp256 = _mm256_max_epi8(tmp256, zero256);
                _mm256_store_si256((__m256i *)(H1 + k*SIMD_WIDTH8), tmp256);
            }           
//------------------------
            /* Banding calculation in pre-processing */
            uint8_t myband[SIMD_WIDTH8] __attribute__((aligned(64)));
            uint8_t temp[SIMD_WIDTH8] __attribute__((aligned(64)));
            {
                __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
                __m256i sum256 = _mm256_add_epi8(qlen256, eb_ins256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH8; l++)
                {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum256 = _mm256_add_epi8(qlen256, eb_del256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH8; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = (int) val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;
                }
            }
            
            smithWaterman256_8(mySeq1SoA,
                               mySeq2SoA,
                               maxLen1,
                               maxLen2,
                               pairArray + i,
                               h0,
                               tid,
                               numPairs,
                               zdrop,
                               bsize,
                               qlen,
                               myband);
        }
    }

#if RDT
    st4 = __rdtsc();
#endif
    
#if SORT_PAIRS      // disbaled in bwa-mem2 (only used in separate benchmark bsw code)
    {
    // Sort the sequences according to increasing order of id
#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsId(pairArray + first, first, last - first, myTempArray);
        }
    }
    _mm_free(tempArray);
    }
#endif

#if RDT
    st5 = __rdtsc();
    setupTicks = st2 - st1;
    sort1Ticks = st3 - st2;
    swTicks = st4 - st3;
    sort2Ticks = st5 - st4;
#endif
    
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    
    return;
}


void BandedPairWiseSW::smithWaterman256_8(uint8_t seq1SoA[],
                                          uint8_t seq2SoA[],
                                          uint8_t nrow,
                                          uint8_t ncol,
                                          SeqPair *p,
                                          uint8_t h0[],
                                          uint16_t tid,
                                          int32_t numPairs,
                                          int zdrop,
                                          int32_t w,
                                          uint8_t qlen[],
                                          uint8_t myband[])
{   
    __m256i match256     = _mm256_set1_epi8(this->w_match);
    __m256i mismatch256  = _mm256_set1_epi8(this->w_mismatch);
    __m256i gapOpen256   = _mm256_set1_epi8(this->w_open);
    __m256i gapExtend256 = _mm256_set1_epi8(this->w_extend);
    __m256i gapOE256     = _mm256_set1_epi8(this->w_open + this->w_extend);
    __m256i w_ambig_256  = _mm256_set1_epi8(this->w_ambig); // ambig penalty
    __m256i five256      = _mm256_set1_epi8(5);

    __m256i e_del256    = _mm256_set1_epi8(this->e_del);
    __m256i oe_del256   = _mm256_set1_epi8(this->o_del + this->e_del);
    __m256i e_ins256    = _mm256_set1_epi8(this->e_ins);
    __m256i oe_ins256   = _mm256_set1_epi8(this->o_ins + this->e_ins);
    
    int8_t  *F  = F8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_h    = H8_ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;
    int8_t  *H_v = H8__ + tid * SIMD_WIDTH8 * MAX_SEQ_LEN8;

    int lane;
    
    int8_t i, j;

    uint8_t tlen[SIMD_WIDTH8];
    uint8_t tail[SIMD_WIDTH8] __attribute((aligned(64)));
    uint8_t head[SIMD_WIDTH8] __attribute((aligned(64)));
    
    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH8; l++) {
        tlen[l] = p[l].len1;
        qlen[l] = p[l].len2;
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m256i tlen256 = _mm256_load_si256((__m256i *) tlen);
    __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
    __m256i myband256 = _mm256_load_si256((__m256i *) myband);
    __m256i zero256 = _mm256_setzero_si256();
    __m256i one256  = _mm256_set1_epi8(1);
    __m256i two256  = _mm256_set1_epi8(2);
    __m256i max_ie256 = zero256;
    __m256i ff256 = _mm256_set1_epi8(0xFF);
        
    __m256i tail256 = qlen256, head256 = zero256;
    _mm256_store_si256((__m256i *) head, head256);
    _mm256_store_si256((__m256i *) tail, tail256);
    //__m256i ib256 = _mm256_add_epi8(qlen256, qlen256);
    // ib256 = _mm256_sub_epi8(ib256, one256);

    __m256i mlen256 = _mm256_add_epi8(qlen256, myband256);
    mlen256 = _mm256_min_epu8(mlen256, tlen256);

    uint8_t temp[SIMD_WIDTH8]  __attribute((aligned(64)));
    uint8_t temp1[SIMD_WIDTH8]  __attribute((aligned(64)));
    
    __m256i s00  = _mm256_load_si256((__m256i *)(seq1SoA));
    __m256i hval = _mm256_load_si256((__m256i *)(H_v));
    __mmask32 dmask = 0xFFFFFFFF;
    
    __m256i maxScore256 = hval;
    for(j = 0; j < ncol; j++)
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), zero256);
    
    __m256i x256 = zero256;
    __m256i y256 = zero256;
    __m256i i256 = zero256;
    __m256i gscore = _mm256_set1_epi8(-1);
    __m256i max_off256 = zero256;
    __m256i exit0 = _mm256_set1_epi8(0xFF);
    __m256i zdrop256 = _mm256_set1_epi8(zdrop);
    
    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;
    
#if RDT
    uint64_t tim = __rdtsc();
#endif
    
    for(i = 0; i < nrow; i++)
    {       
        __m256i e11 = zero256;
        __m256i h00, h11, h10;
        __m256i s10 = _mm256_load_si256((__m256i *)(seq1SoA + (i + 0) * SIMD_WIDTH8));
        
        beg = nbeg; end = nend;
        int pbeg = beg;
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;
        
        h10 = zero256;
        if (beg == 0)
            h10 = _mm256_load_si256((__m256i *)(H_v + (i+1) * SIMD_WIDTH8));
        
        __m256i j256 = zero256;
        __m256i maxRS1;
        maxRS1 = zero256;
        
        __m256i i1_256 = _mm256_set1_epi8(i+1);
        __m256i y1_256 = zero256;
        
#if RDT 
        uint64_t tim1 = __rdtsc();
#endif
        
        __m256i i256, cache256;
        __m256i phead256 = head256, ptail256 = tail256;
        i256 = _mm256_set1_epi8(i);
        cache256 = _mm256_sub_epi8(i256, myband256);
        head256 = _mm256_max_epi8(head256, cache256);
        cache256 = _mm256_add_epi8(i1_256, myband256);
        tail256 = _mm256_min_epu8(tail256, cache256);
        tail256 = _mm256_min_epu8(tail256, qlen256);
        
        // NEW, trimming.
        __m256i cmph = _mm256_cmpeq_epi8(head256, phead256);
        __m256i cmpt = _mm256_cmpeq_epi8(tail256, ptail256);
        // cmph &= cmpt;
        cmph = _mm256_and_si256(cmph, cmpt);
        __mmask32 cmp_ht = _mm256_movemask_epi8(cmph);
        
        for (int l=beg; l<end && cmp_ht != dmask; l++)
        {
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            
            __m256i pj256 = _mm256_set1_epi8(l);
            __m256i j256 = _mm256_set1_epi8(l+1);
            __m256i cmp1 = _mm256_cmpgt_epi8(head256, pj256);
            uint32_t cval = _mm256_movemask_epi8(cmp1);
            if (cval == 0x00) break;
            //__m256i cmp2 = _mm256_cmpgt_epi8(pj256, tail256);
            __m256i cmp2 = _mm256_cmpgt_epi8(j256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h256 = _mm256_blendv_epi8(h256, zero256, cmp1);
            f256 = _mm256_blendv_epi8(f256, zero256, cmp1);
            
            _mm256_store_si256((__m256i *)(F + l * SIMD_WIDTH8), f256);
            _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH8), h256);
        }
        
#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif
        // beg = nbeg; end = nend;
        //__m256i cmp256_1 = _mm256_cmpgt_epi8(i1_256, tlen256);
        
        // beg = nbeg; end = nend;
        __m256i cmp256_1 = _mm256_cmpgt_epi8(i1_256, tlen256);
        
        __m256i cmpim = _mm256_cmpgt_epi8(i1_256, mlen256);
        __m256i cmpht = _mm256_cmpeq_epi8(tail256, head256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        // NEW
        cmpht = _mm256_cmpgt_epi8(head256, tail256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        exit0 = _mm256_blendv_epi8(exit0, zero256, cmpim);
        
        
#if RDT
        tim1 = __rdtsc();
#endif
        
        j256 = _mm256_set1_epi8(beg);
#pragma unroll(4)
        for(j = beg; j < end; j++)
        {
            __m256i f11, f21, s2;
            h00 = _mm256_load_si256((__m256i *)(H_h + j * SIMD_WIDTH8));
            f11 = _mm256_load_si256((__m256i *)(F + j * SIMD_WIDTH8));
            
            s2 = _mm256_load_si256((__m256i *)(seq2SoA + (j) * SIMD_WIDTH8));
            
            __m256i pj256 = j256;
            j256 = _mm256_add_epi8(j256, one256);
            
            MAIN_CODE8(s10, s2, h00, h11, e11, f11, f21, zero256,
                       maxScore256, e_ins256, oe_ins256,
                       e_del256, oe_del256,
                       y1_256, maxRS1); //i+1

            
            // Masked writing
            __m256i cmp2 = _mm256_cmpgt_epi8(head256, pj256);
            __m256i cmp1 = _mm256_cmpgt_epi8(pj256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h10 = _mm256_blendv_epi8(h10, zero256, cmp1);
            f21 = _mm256_blendv_epi8(f21, zero256, cmp1);
            
            __m256i bmaxRS = maxRS1;                                        
            maxRS1 =_mm256_max_epi8(maxRS1, h11);                           
            __m256i cmpA = _mm256_cmpgt_epi8(maxRS1, bmaxRS);                   
            __m256i cmpB =_mm256_cmpeq_epi8(maxRS1, h11);                   
            cmpA = _mm256_or_si256(cmpA, cmpB);
            cmp1 = _mm256_cmpgt_epi8(j256, tail256); // change
            cmp1 = _mm256_or_si256(cmp1, cmp2);       // change  
            cmpA = _mm256_blendv_epi8(y1_256, j256, cmpA);
            y1_256 = _mm256_blendv_epi8(cmpA, y1_256, cmp1);
            maxRS1 = _mm256_blendv_epi8(maxRS1, bmaxRS, cmp1);                      
            
            _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), f21);
            _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH8), h10);
            
            h10 = h11;
            
            //j256 = _mm256_add_epi8(j256, one256);
            
            if (j >= minq)
            {
                __m256i cmp = _mm256_cmpeq_epi8(j256, qlen256);
                __m256i max_gh = _mm256_max_epi8(gscore, h11);
                __m256i cmp_gh = _mm256_cmpgt_epi8(gscore, h11);
                __m256i tmp256_1 = _mm256_blendv_epi8(i1_256, max_ie256, cmp_gh);
                
                tmp256_1 = _mm256_blendv_epi8(max_ie256, tmp256_1, cmp);
                tmp256_1 = _mm256_blendv_epi8(max_ie256, tmp256_1, exit0);
                
                max_gh = _mm256_blendv_epi8(gscore, max_gh, exit0);
                max_gh = _mm256_blendv_epi8(gscore, max_gh, cmp);               
                
                cmp = _mm256_cmpgt_epi8(j256, tail256); 
                max_gh = _mm256_blendv_epi8(max_gh, gscore, cmp);
                max_ie256 = _mm256_blendv_epi8(tmp256_1, max_ie256, cmp);
                gscore = max_gh;
            }
        }
        _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH8), h10);
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH8), zero256);
        
        
        /* exit due to zero score by a row */
        uint32_t cval = 0;
        __m256i bmaxScore256 = maxScore256;
        __m256i tmp = _mm256_cmpeq_epi8(maxRS1, zero256);
        cval = _mm256_movemask_epi8(tmp);
        if (cval == 0xFFFFFFFF) break;
        exit0 = _mm256_blendv_epi8(exit0, zero256,  tmp);
        
        __m256i score256 = _mm256_max_epi8(maxScore256, maxRS1);
        maxScore256 = _mm256_blendv_epi8(maxScore256, score256, exit0);
        
        __m256i cmp = _mm256_cmpgt_epi8(maxScore256, bmaxScore256);
        y256 = _mm256_blendv_epi8(y256, y1_256, cmp);
        x256 = _mm256_blendv_epi8(x256, i1_256, cmp);       
        // max_off calculations
        tmp = _mm256_sub_epi8(y1_256, i1_256);
        tmp = _mm256_abs_epi8(tmp);
        __m256i bmax_off256 = max_off256;
        tmp = _mm256_max_epi8(max_off256, tmp);
        max_off256 = _mm256_blendv_epi8(bmax_off256, tmp, cmp);
        
        // Z-score
        ZSCORE8(i1_256, y1_256);        
        
#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
        tim1 = __rdtsc();
#endif


        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val == 0xFFFFFFFF) nbeg = l;
            else
                break;
        }
        
        /* From end */
        bool flg = 1;
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val != 0xFFFFFFFF && flg)  
                break;
        }
        // int pnend =nend;
        nend = l + 2 < ncol? l + 2: ncol;
        
        __m256i tail256_ = _mm256_sub_epi8(tail256, one256);
        __m256i tmpb = ff256;
        
        __m256i exit1 = _mm256_xor_si256(exit0, ff256);
        __m256i l256 = _mm256_set1_epi8(beg);
        for (l = beg; l < end; l++)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);          
            tmp = _mm256_cmpeq_epi8(tmp, zero256);
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val == 0x00) {
                break;
            }
            tmp = _mm256_and_si256(tmp,tmpb);
            //__m256i l256 = _mm256_set1_epi8(l+1);
            l256 = _mm256_add_epi8(l256, one256);
            
            head256 = _mm256_blendv_epi8(head256, l256, tmp);
            tmpb = tmp;         
        }
        // _mm256_store_si256((__m256i *) head, head256);
        
        __m256i  index256 = tail256;
        tmpb = ff256;
        
        l256 = _mm256_set1_epi8(end);
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH8));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH8));
            
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);
            tmp = _mm256_cmpeq_epi8(tmp, zero256);          
            uint32_t val = _mm256_movemask_epi8(tmp);
            if (val == 0x00)  {
                break;
            }
            
            tmp = _mm256_and_si256(tmp,tmpb);
            l256 = _mm256_sub_epi8(l256, one256);
            
            index256 = _mm256_blendv_epi8(index256, l256, tmp);
            tmpb = tmp;
        }
        index256 = _mm256_add_epi8(index256, two256);
        tail256 = _mm256_min_epi8(index256, qlen256);
        // _mm256_store_si256((__m256i *) tail, tail256);       

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
    
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif
    
    int8_t score[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) score, maxScore256);

    int8_t maxi[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxi, x256);

    int8_t maxj[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxj, y256);

    int8_t max_off_ar[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) max_off_ar, max_off256);

    int8_t gscore_ar[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) gscore_ar, gscore);

    int8_t maxie_ar[SIMD_WIDTH8]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxie_ar, max_ie256);
    
    for(i = 0; i < SIMD_WIDTH8; i++)
    {
        p[i].score = score[i];
        p[i].tle = maxi[i];
        p[i].qle = maxj[i];
        p[i].max_off = max_off_ar[i];
        p[i].gscore = gscore_ar[i];
        p[i].gtle = maxie_ar[i];
    }
    
    return;
}

// ------------------------- AVX2 - 16 bit SIMD_LANES ---------------------------
#define PFD 2
void BandedPairWiseSW::getScores16(SeqPair *pairArray,
                                   uint8_t *seqBufRef,
                                   uint8_t *seqBufQer,
                                   int32_t numPairs,
                                   uint16_t numThreads,
                                   int32_t w)
{
    int64_t startTick, endTick;

    smithWatermanBatchWrapper16(pairArray, seqBufRef, seqBufQer, numPairs, numThreads, w);


#if MAXI
    printf("AVX2 Vecor code: Writing output..\n");
    for (int l=0; l<numPairs; l++)
    {
        fprintf(stderr, "%d (%d %d) %d %d %d\n",
                pairArray[l].score, pairArray[l].x, pairArray[l].y,
                pairArray[l].gscore, pairArray[l].max_off, pairArray[l].max_ie);

    }
    printf("Vector code: Writing output completed!!!\n\n");
#endif
    
}

void BandedPairWiseSW::smithWatermanBatchWrapper16(SeqPair *pairArray,
                                                   uint8_t *seqBufRef,
                                                   uint8_t *seqBufQer,
                                                   int32_t numPairs,
                                                   uint16_t numThreads,
                                                   int32_t w)
{
    int64_t st1, st2, st3, st4, st5;
#if RDT     
    st1 = __rdtsc();
#endif
    
    uint16_t *seq1SoA = (uint16_t *)_mm_malloc(MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);
    uint16_t *seq2SoA = (uint16_t *)_mm_malloc(MAX_SEQ_LEN16 * SIMD_WIDTH16 * numThreads * sizeof(uint16_t), 64);

    if (seq1SoA == NULL || seq2SoA == NULL) {
        fprintf(stderr, "Error! Mem not allocated!!!\n");
        exit(EXIT_FAILURE);
    }
    
    int32_t ii;
    int32_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1)/SIMD_WIDTH16 ) * SIMD_WIDTH16;
    for(ii = numPairs; ii < roundNumPairs; ii++)
    {
        pairArray[ii].id = ii;
        pairArray[ii].len1 = 0;
        pairArray[ii].len2 = 0;
    }

#if RDT 
    st2 = __rdtsc();
#endif
    
#if SORT_PAIRS      // disbaled in bwa-mem2 (only used in separate benchmark bsw code)
    // Sort the sequences according to decreasing order of lengths
    SeqPair *tempArray = (SeqPair *)_mm_malloc(SORT_BLOCK_SIZE * numThreads *
                                               sizeof(SeqPair), 64);
    int16_t *hist = (int16_t *)_mm_malloc((MAX_SEQ_LEN16 + 16) * numThreads *
                                          sizeof(int16_t), 64);
    // int16_t *histb = (int16_t *)_mm_malloc((MAX_SEQ_LEN16 + 16) * numThreads *
    //                                        sizeof(int16_t), 64);
#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;
        int16_t *myHist = hist + tid * (MAX_SEQ_LEN16 + 16);
        // int16_t *myHistb = histb + tid * (MAX_SEQ_LEN16 + 16);

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            // sortPairsLen(pairArray + first, last - first, myTempArray, myHist, myHistb);
            sortPairsLen(pairArray + first, last - first, myTempArray, myHist);
        }
    }
    _mm_free(hist);
#endif

#if RDT 
    st3 = __rdtsc();
#endif

    int eb = end_bonus;
#pragma omp parallel num_threads(numThreads)
    {
        int32_t i;
        uint16_t tid = omp_get_thread_num(); 
        uint16_t *mySeq1SoA = seq1SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint16_t *mySeq2SoA = seq2SoA + tid * MAX_SEQ_LEN16 * SIMD_WIDTH16;
        uint8_t *seq1;
        uint8_t *seq2;
        uint16_t h0[SIMD_WIDTH16]   __attribute__((aligned(64)));
        uint16_t band[SIMD_WIDTH16];        
        uint16_t qlen[SIMD_WIDTH16] __attribute__((aligned(64)));
        int32_t bsize = 0;
        
        int16_t *H1 = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
        int16_t *H2 = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
        
        __m256i zero256   = _mm256_setzero_si256();
        __m256i e_ins256  = _mm256_set1_epi16(e_ins);
        __m256i oe_ins256 = _mm256_set1_epi16(o_ins + e_ins);
        __m256i o_del256  = _mm256_set1_epi16(o_del);
        __m256i e_del256  = _mm256_set1_epi16(e_del);
        __m256i eb_ins256 = _mm256_set1_epi16(eb - o_ins);
        __m256i eb_del256 = _mm256_set1_epi16(eb - o_del);
        
        int16_t max = 0;
        if (max < w_match) max = w_match;
        if (max < w_mismatch) max = w_mismatch;
        if (max < w_ambig) max = w_ambig;
        
        int nstart = 0, nend = numPairs;

#pragma omp for schedule(dynamic, 128)
        for(i = nstart; i < nend; i+=SIMD_WIDTH16)
        {
            int32_t j, k;
            uint16_t maxLen1 = 0;
            uint16_t maxLen2 = 0;
            bsize = w;

            uint64_t tim;
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                { // prefetch block
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufRef + (int64_t)spf.idr + 64, _MM_HINT_NTA);
                }

                SeqPair sp = pairArray[i + j];
                h0[j] = sp.h0;
                seq1 = seqBufRef + (int64_t)sp.idr;
                
                for(k = 0; k < sp.len1; k++)
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = (seq1[k] == AMBIG?0xFFFF:seq1[k]);
                    H2[k * SIMD_WIDTH16 + j] = 0;
                }
                qlen[j] = sp.len2 * max;
                if(maxLen1 < sp.len1) maxLen1 = sp.len1;
            }
        
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len1; k <= maxLen1; k++) //removed "="
                {
                    mySeq1SoA[k * SIMD_WIDTH16 + j] = DUMMY1;
                    H2[k * SIMD_WIDTH16 + j] = DUMMY1;
                }
            }
//--------------------
            __m256i h0_256 = _mm256_load_si256((__m256i*) h0);
            _mm256_store_si256((__m256i *) H2, h0_256);
            __m256i tmp256 = _mm256_sub_epi16(h0_256, o_del256);

            for(k = 1; k < maxLen1; k++) {
                tmp256 = _mm256_sub_epi16(tmp256, e_del256);
                __m256i tmp256_ = _mm256_max_epi16(tmp256, zero256);
                _mm256_store_si256((__m256i *)(H2 + k* SIMD_WIDTH16), tmp256_);
            }
//-------------------
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                { // prefetch block
                    SeqPair spf = pairArray[i + j + PFD];
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq, _MM_HINT_NTA);
                    _mm_prefetch((const char*) seqBufQer + (int64_t)spf.idq + 64, _MM_HINT_NTA);
                }
                
                SeqPair sp = pairArray[i + j];
                //seq2 = seqBufQer + (int64_t)sp.id * MAX_SEQ_LEN_QER;
                seq2 = seqBufQer + (int64_t)sp.idq;             
                for(k = 0; k < sp.len2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = (seq2[k]==AMBIG?0xFFFF:seq2[k]);
                    H1[k * SIMD_WIDTH16 + j] = 0;                   
                }
                if(maxLen2 < sp.len2) maxLen2 = sp.len2;
            }
            
            for(j = 0; j < SIMD_WIDTH16; j++)
            {
                SeqPair sp = pairArray[i + j];
                for(k = sp.len2; k <= maxLen2; k++)
                {
                    mySeq2SoA[k * SIMD_WIDTH16 + j] = DUMMY2;
                    H1[k * SIMD_WIDTH16 + j] = 0;
                }
            }
//------------------------
            _mm256_store_si256((__m256i *) H1, h0_256);
            __m256i cmp256 = _mm256_cmpgt_epi16(h0_256, oe_ins256);
            tmp256 = _mm256_sub_epi16(h0_256, oe_ins256);

            tmp256 = _mm256_blendv_epi16(zero256, tmp256, cmp256);
            _mm256_store_si256((__m256i *) (H1 + SIMD_WIDTH16), tmp256);
            for(k = 2; k < maxLen2; k++)
            {
                __m256i h1_256 = tmp256;
                tmp256 = _mm256_sub_epi16(h1_256, e_ins256);
                tmp256 = _mm256_max_epi16(tmp256, zero256);
                _mm256_store_si256((__m256i *)(H1 + k*SIMD_WIDTH16), tmp256);
            }
//------------------------
            uint16_t myband[SIMD_WIDTH16] __attribute__((aligned(64)));
            uint16_t temp[SIMD_WIDTH16] __attribute__((aligned(64)));
            {
                __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
                __m256i sum256 = _mm256_add_epi16(qlen256, eb_ins256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_ins + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(bsize, max_ins);
                }
                sum256 = _mm256_add_epi16(qlen256, eb_del256);
                _mm256_store_si256((__m256i *) temp, sum256);               
                for (int l=0; l<SIMD_WIDTH16; l++) {
                    double val = temp[l]/e_del + 1.0;
                    int max_ins = val;
                    max_ins = max_ins > 1? max_ins : 1;
                    myband[l] = min_(myband[l], max_ins);
                    bsize = bsize < myband[l] ? myband[l] : bsize;                  
                }
            }

            smithWaterman256_16(mySeq1SoA,
                                mySeq2SoA,
                                maxLen1,
                                maxLen2,
                                pairArray + i,
                                h0,
                                tid,
                                numPairs,
                                zdrop,
                                bsize, 
                                qlen,
                                myband);
        }
    }

#if RDT
    st4 = __rdtsc();
#endif
    
#if SORT_PAIRS      // disbaled in bwa-mem2 (only used in separate benchmark bsw code)
    {
    // Sort the sequences according to increasing order of id
#pragma omp parallel num_threads(numThreads)
    {
        int32_t tid = omp_get_thread_num();
        SeqPair *myTempArray = tempArray + tid * SORT_BLOCK_SIZE;

#pragma omp for
        for(ii = 0; ii < roundNumPairs; ii+=SORT_BLOCK_SIZE)
        {
            int32_t first, last;
            first = ii;
            last  = ii + SORT_BLOCK_SIZE;
            if(last > roundNumPairs) last = roundNumPairs;
            sortPairsId(pairArray + first, first, last - first, myTempArray);
        }
    }
    _mm_free(tempArray);
    }
#endif

#if RDT
    st5 = __rdtsc();
    setupTicks += st2 - st1;
    sort1Ticks += st3 - st2;
    swTicks += st4 - st3;
    sort2Ticks += st5 - st4;
#endif
    
    // free mem
    _mm_free(seq1SoA);
    _mm_free(seq2SoA);
    
    return;
}

void BandedPairWiseSW::smithWaterman256_16(uint16_t seq1SoA[],
                                           uint16_t seq2SoA[],
                                           uint16_t nrow,
                                           uint16_t ncol,
                                           SeqPair *p,
                                           uint16_t h0[],
                                           uint16_t tid,
                                           int32_t numPairs,
                                           int zdrop,
                                           int32_t w,
                                           uint16_t qlen[],
                                           uint16_t myband[])
{
    __m256i match256     = _mm256_set1_epi16(this->w_match);
    __m256i mismatch256  = _mm256_set1_epi16(this->w_mismatch);
    __m256i gapOpen256   = _mm256_set1_epi16(this->w_open);
    __m256i gapExtend256 = _mm256_set1_epi16(this->w_extend);
    __m256i gapOE256     = _mm256_set1_epi16(this->w_open + this->w_extend);
    __m256i w_ambig_256  = _mm256_set1_epi16(this->w_ambig);    // ambig penalty
    __m256i five256      = _mm256_set1_epi16(5);

    __m256i e_del256    = _mm256_set1_epi16(this->e_del);
    __m256i oe_del256   = _mm256_set1_epi16(this->o_del + this->e_del);
    __m256i e_ins256    = _mm256_set1_epi16(this->e_ins);
    __m256i oe_ins256   = _mm256_set1_epi16(this->o_ins + this->e_ins);
    
    int16_t *F  = F16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_h    = H16_ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;
    int16_t *H_v = H16__ + tid * SIMD_WIDTH16 * MAX_SEQ_LEN16;

    int lane = 0;
    
    int16_t i, j;

    uint16_t tlen[SIMD_WIDTH16];
    uint16_t tail[SIMD_WIDTH16] __attribute((aligned(64)));
    uint16_t head[SIMD_WIDTH16] __attribute((aligned(64)));
    
    int32_t minq = 10000000;
    for (int l=0; l<SIMD_WIDTH16; l++) {
        tlen[l] = p[l].len1;
        qlen[l] = p[l].len2;
        if (p[l].len2 < minq) minq = p[l].len2;
    }
    minq -= 1; // for gscore

    __m256i tlen256 = _mm256_load_si256((__m256i *) tlen);
    __m256i qlen256 = _mm256_load_si256((__m256i *) qlen);
    __m256i myband256 = _mm256_load_si256((__m256i *) myband);
    __m256i zero256 = _mm256_setzero_si256();
    __m256i one256  = _mm256_set1_epi16(1);
    __m256i two256  = _mm256_set1_epi16(2);
    __m256i max_ie256 = zero256;
    __m256i ff256 = _mm256_set1_epi16(0xFFFF);
        
    __m256i tail256 = qlen256, head256 = zero256;
    _mm256_store_si256((__m256i *) head, head256);
    _mm256_store_si256((__m256i *) tail, tail256);

    __m256i mlen256 = _mm256_add_epi16(qlen256, myband256);
    mlen256 = _mm256_min_epu16(mlen256, tlen256);

    uint16_t temp[SIMD_WIDTH16]  __attribute((aligned(64)));
    uint16_t temp1[SIMD_WIDTH16]  __attribute((aligned(64)));
    
    __m256i s00  = _mm256_load_si256((__m256i *)(seq1SoA));
    __m256i hval = _mm256_load_si256((__m256i *)(H_v));
    __mmask16 dmask = 0xFFFF;
    __mmask32 dmask32 = 0xAAAAAAAA;
        
    __m256i maxScore256 = hval;
    for(j = 0; j < ncol; j++)
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), zero256);
    
    __m256i x256 = zero256;
    __m256i y256 = zero256;
    __m256i i256 = zero256;
    __m256i gscore = _mm256_set1_epi16(-1);
    __m256i max_off256 = zero256;
    __m256i exit0 = _mm256_set1_epi16(0xFFFF);
    __m256i zdrop256 = _mm256_set1_epi16(zdrop);
    
    int beg = 0, end = ncol;
    int nbeg = beg, nend = end;

#if RDT
    uint64_t tim = __rdtsc();
#endif
    
    for(i = 0; i < nrow; i++)
    {       
        __m256i e11 = zero256;
        __m256i h00, h11, h10;
        __m256i s10 = _mm256_load_si256((__m256i *)(seq1SoA + (i + 0) * SIMD_WIDTH16));

        beg = nbeg; end = nend;
        int pbeg = beg;
        if (beg < i - w) beg = i - w;
        if (end > i + w + 1) end = i + w + 1;
        if (end > ncol) end = ncol;

        h10 = zero256;
        if (beg == 0)
            h10 = _mm256_load_si256((__m256i *)(H_v + (i+1) * SIMD_WIDTH16));

        __m256i j256 = zero256;
        __m256i maxRS1;
        maxRS1 = zero256;

        __m256i i1_256 = _mm256_set1_epi16(i+1);
        __m256i y1_256 = zero256;
        
#if RDT 
        uint64_t tim1 = __rdtsc();
#endif
        
        __m256i i256, cache256;
        __m256i phead256 = head256, ptail256 = tail256;
        i256 = _mm256_set1_epi16(i);
        cache256 = _mm256_sub_epi16(i256, myband256);
        head256 = _mm256_max_epi16(head256, cache256);
        cache256 = _mm256_add_epi16(i1_256, myband256);
        tail256 = _mm256_min_epu16(tail256, cache256);
        tail256 = _mm256_min_epu16(tail256, qlen256);

        // NEW, trimming.
        __m256i cmph = _mm256_cmpeq_epi16(head256, phead256);
        __m256i cmpt = _mm256_cmpeq_epi16(tail256, ptail256);
        // cmph &= cmpt;
        cmph = _mm256_and_si256(cmph, cmpt);
        //__mmask16 cmp_ht = _mm256_movepi16_mask(cmph);
        __mmask32 cmp_ht = _mm256_movemask_epi8(cmph) & dmask32;
        
        for (int l=beg; l<end && cmp_ht != dmask32; l++)
        {
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            
            __m256i pj256 = _mm256_set1_epi16(l);
            __m256i j256 = _mm256_set1_epi16(l+1);
            __m256i cmp1 = _mm256_cmpgt_epi16(head256, pj256);
            //uint16_t cval = _mm256_movepi16_mask(cmp1);
            uint32_t cval = _mm256_movemask_epi8(cmp1) & dmask32;
            if (cval == 0x00) break;
            //__m256i cmp2 = _mm256_cmpgt_epi16(pj256, tail256);
            __m256i cmp2 = _mm256_cmpgt_epi16(j256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h256 = _mm256_blendv_epi16(h256, zero256, cmp1);
            f256 = _mm256_blendv_epi16(f256, zero256, cmp1);
            
            _mm256_store_si256((__m256i *)(F + l * SIMD_WIDTH16), f256);
            _mm256_store_si256((__m256i *)(H_h + l * SIMD_WIDTH16), h256);
        }

#if RDT
        prof[DP3][0] += __rdtsc() - tim1;
#endif

        // beg = nbeg; end = nend;
        __m256i cmp256_1 = _mm256_cmpgt_epi16(i1_256, tlen256);
        
        __m256i cmpim = _mm256_cmpgt_epi16(i1_256, mlen256);
        __m256i cmpht = _mm256_cmpeq_epi16(tail256, head256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        // NEW
        cmpht = _mm256_cmpgt_epi16(head256, tail256);
        cmpim = _mm256_or_si256(cmpim, cmpht);

        exit0 = _mm256_blendv_epi16(exit0, zero256, cmpim);

        
#if RDT
        tim1 = __rdtsc();
#endif
        
        j256 = _mm256_set1_epi16(beg);
        for(j = beg; j < end; j++)
        {
            __m256i f11, f21, s2;
            h00 = _mm256_load_si256((__m256i *)(H_h + j * SIMD_WIDTH16));
            f11 = _mm256_load_si256((__m256i *)(F + j * SIMD_WIDTH16));

            s2 = _mm256_load_si256((__m256i *)(seq2SoA + (j) * SIMD_WIDTH16));
            
            __m256i pj256 = j256;
            j256 = _mm256_add_epi16(j256, one256);

            MAIN_CODE16(s10, s2, h00, h11, e11, f11, f21, zero256,
                        maxScore256, e_ins256, oe_ins256,
                        e_del256, oe_del256,
                        y1_256, maxRS1); //i+1
            
            // Masked writing
            __m256i cmp2 = _mm256_cmpgt_epi16(head256, pj256);
            __m256i cmp1 = _mm256_cmpgt_epi16(pj256, tail256);
            cmp1 = _mm256_or_si256(cmp1, cmp2);
            h10 = _mm256_blendv_epi16(h10, zero256, cmp1);
            f21 = _mm256_blendv_epi16(f21, zero256, cmp1);
            
            __m256i bmaxRS = maxRS1;                                        
            maxRS1 =_mm256_max_epi16(maxRS1, h11);                          
            __m256i cmpA = _mm256_cmpgt_epi16(maxRS1, bmaxRS);                  
            __m256i cmpB =_mm256_cmpeq_epi16(maxRS1, h11);                  
            cmpA = _mm256_or_si256(cmpA, cmpB);
            cmp1 = _mm256_cmpgt_epi16(j256, tail256); // change
            cmp1 = _mm256_or_si256(cmp1, cmp2);         // change
            cmpA = _mm256_blendv_epi16(y1_256, j256, cmpA);
            y1_256 = _mm256_blendv_epi16(cmpA, y1_256, cmp1);
            maxRS1 = _mm256_blendv_epi16(maxRS1, bmaxRS, cmp1);                     

            _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), f21);
            _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH16), h10);

            h10 = h11;
            
            //j256 = _mm256_add_epi16(j256, one256);
            
            // gscore calculations
            if (j >= minq)
            {
                __m256i cmp = _mm256_cmpeq_epi16(j256, qlen256);
                __m256i max_gh = _mm256_max_epi16(gscore, h11);
                __m256i cmp_gh = _mm256_cmpgt_epi16(gscore, h11);
                __m256i tmp256_1 = _mm256_blendv_epi16(i1_256, max_ie256, cmp_gh);

                __m256i tmp256_t = _mm256_blendv_epi16(max_ie256, tmp256_1, cmp);
                tmp256_1 = _mm256_blendv_epi16(max_ie256, tmp256_t, exit0);             

                max_gh = _mm256_blendv_epi16(gscore, max_gh, exit0);
                max_gh = _mm256_blendv_epi16(gscore, max_gh, cmp);              

                cmp = _mm256_cmpgt_epi16(j256, tail256); 
                max_gh = _mm256_blendv_epi16(max_gh, gscore, cmp);
                max_ie256 = _mm256_blendv_epi16(tmp256_1, max_ie256, cmp);
                gscore = max_gh;            
            }
        }
        _mm256_store_si256((__m256i *)(H_h + j * SIMD_WIDTH16), h10);
        _mm256_store_si256((__m256i *)(F + j * SIMD_WIDTH16), zero256);
                
        /* exit due to zero score by a row */
        __m256i bmaxScore256 = maxScore256;
        __m256i tmp = _mm256_cmpeq_epi16(maxRS1, zero256);
        uint32_t cval = _mm256_movemask_epi8(tmp) & dmask32;
        if (cval == dmask32) break;

        exit0 = _mm256_blendv_epi16(exit0, zero256,  tmp);

        __m256i score256 = _mm256_max_epi16(maxScore256, maxRS1);
        maxScore256 = _mm256_blendv_epi16(maxScore256, score256, exit0);

        __m256i cmp = _mm256_cmpgt_epi16(maxScore256, bmaxScore256);
        y256 = _mm256_blendv_epi16(y256, y1_256, cmp);
        x256 = _mm256_blendv_epi16(x256, i1_256, cmp);      
        // max_off calculations
        tmp = _mm256_sub_epi16(y1_256, i1_256);
        tmp = _mm256_abs_epi16(tmp);
        __m256i bmax_off256 = max_off256;
        tmp = _mm256_max_epi16(max_off256, tmp);
        max_off256 = _mm256_blendv_epi16(bmax_off256, tmp, cmp);

        // Z-score
        ZSCORE16(i1_256, y1_256);       

#if RDT
        prof[DP1][0] += __rdtsc() - tim1;
        tim1 = __rdtsc();
#endif
        
        /* Narrowing of the band */
        /* From beg */
        int l;
        for (l = beg; l < end; l++) {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi16(tmp, zero256);
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val == dmask32) nbeg = l;
            else
                break;
        }
        
        /* From end */
        bool flg = 1;
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_cmpeq_epi16(tmp, zero256);
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val != dmask32 && flg)  
                break;
        }
        nend = l + 2 < ncol? l + 2: ncol;

        __m256i tail256_ = _mm256_sub_epi16(tail256, one256);
        __m256i tmpb = ff256;

        __m256i exit1 = _mm256_xor_si256(exit0, ff256);
        __m256i l256 = _mm256_set1_epi16(beg);
        for (l = beg; l < end; l++)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
    
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);          
            tmp = _mm256_cmpeq_epi16(tmp, zero256);
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val == 0x00) {
                break;
            }
            tmp = _mm256_and_si256(tmp,tmpb);
            //__m256i l256 = _mm256_set1_epi16(l+1);
            l256 = _mm256_add_epi16(l256, one256);

            head256 = _mm256_blendv_epi16(head256, l256, tmp);

            tmpb = tmp;         
        }
        // _mm256_store_si256((__m256i *) head, head256);
        
        __m256i  index256 = tail256;
        tmpb = ff256;

        l256 = _mm256_set1_epi16(end);
        for (l = end; l >= beg; l--)
        {
            __m256i f256 = _mm256_load_si256((__m256i *)(F + l * SIMD_WIDTH16));
            __m256i h256 = _mm256_load_si256((__m256i *)(H_h + l * SIMD_WIDTH16));
            
            __m256i tmp = _mm256_or_si256(f256, h256);
            tmp = _mm256_or_si256(tmp, exit1);
            tmp = _mm256_cmpeq_epi16(tmp, zero256);         
            //uint16_t val = _mm256_movepi16_mask(tmp);
            uint32_t val = _mm256_movemask_epi8(tmp) & dmask32;
            if (val == 0x00)  {
                break;
            }
            tmp = _mm256_and_si256(tmp,tmpb);
            l256 = _mm256_sub_epi16(l256, one256);

            // NEW
            index256 = _mm256_blendv_epi8(index256, l256, tmp);

            tmpb = tmp;
        }
        index256 = _mm256_add_epi16(index256, two256);
        tail256 = _mm256_min_epi16(index256, qlen256);
        // _mm256_store_si256((__m256i *) tail, tail256);       

#if RDT
        prof[DP2][0] += __rdtsc() - tim1;
#endif
    }
    
#if RDT
    prof[DP][0] += __rdtsc() - tim;
#endif
    
    int16_t score[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) score, maxScore256);

    int16_t maxi[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxi, x256);

    int16_t maxj[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxj, y256);

    int16_t max_off_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) max_off_ar, max_off256);

    int16_t gscore_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) gscore_ar, gscore);

    int16_t maxie_ar[SIMD_WIDTH16]  __attribute((aligned(64)));
    _mm256_store_si256((__m256i *) maxie_ar, max_ie256);
    
    for(i = 0; i < SIMD_WIDTH16; i++)
    {
        p[i].score = score[i];
        p[i].tle = maxi[i];
        p[i].qle = maxj[i];
        p[i].max_off = max_off_ar[i];
        p[i].gscore = gscore_ar[i];
        p[i].gtle = maxie_ar[i];
    }
    
    return;
}

#endif // AVX2