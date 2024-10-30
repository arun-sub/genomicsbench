#include <vector>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <cmath>
#include "omp.h"
#include "host_kernel.h"
#include "common.h"
// #include "minimap.h"
// #include "mmpriv.h"
// #include "kalloc.h"

#ifdef __AVX2__
    #include <immintrin.h>
#endif
#ifdef __AVX512BW__
    // #include <zmmintrin.h>
#endif
#ifdef __ARM_FEATURE_SVE
    #include <arm_sve.h>
    #define VL svcntw()
#endif

static const signed char LogTable256_dp_lib[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
    -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
        LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
        LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32_dp_lib(uint32_t v) {
    uint32_t t, tt;
    if ((tt = v >> 16)) {
        return (t = tt >> 8) ? 24 + LogTable256_dp_lib[t] : 16 + LogTable256_dp_lib[tt];
    }
    return (t = v >> 8) ? 8 + LogTable256_dp_lib[t] : LogTable256_dp_lib[v];
}

static inline int32_t ilog2_32(uint32_t v) {
    constexpr uint32_t base = 31;
    const uint32_t leading_zeros = __builtin_clz(v);

    return base - leading_zeros;
}

const int BACKSEARCH = 65;
#define MM_SEED_SEG_SHIFT  48
#define MM_SEED_SEG_MASK   (0xffULL<<(MM_SEED_SEG_SHIFT))

#if 0
void print_vector(svint64_t v) {
    int64_t _v[VL];
    svst1_s64(svptrue_b64(),_v,v);
    for (uint64_t i = 0; i < VL; i++) {
        printf("%ld,",_v[i]);
    }
    printf("\n");
}

void print_vector_f(svfloat64_t v) {
    float64_t _v[VL];
    svst1_f64(svptrue_b64(),_v,v);
    for (uint64_t i = 0; i < VL; i++) {
        printf("%f,",_v[i]);
    }
    printf("\n");
}
#endif

#ifdef __AVX512BW__
    inline __m512i get_gap_cost_vectorized_int32(__m512i dd_v, float avg_qspan, float gap_scale) {
        //Vectorized log2
        uint32_t base = 31;
        __m512i vbase = _mm512_set1_epi32(base);
        __mmask16 msk = 0xFFFF;
        __m512i vout  = _mm512_mask_lzcnt_epi32(dd_v, msk, dd_v);
        //int res = base - temp[0];
        __m512i r_v = _mm512_sub_epi32(vbase, vout);

        // log_dd = dd?ilog2:0; log_dd>>1
        __m512i zero_v = _mm512_setzero_si512();
        __mmask16 neg_mask = _mm512_cmpneq_epi32_mask(dd_v, zero_v);
        __m512i log_dd_v = _mm512_srli_epi32(_mm512_maskz_or_epi32(neg_mask, r_v, zero_v), 1);

        //dd * 0.01*avg_qspan
        float avg_qspan_val = 0.01 * avg_qspan;
        __m512 avg_qspan_v = _mm512_set1_ps(avg_qspan_val);
        __m512i cost = _mm512_cvt_roundps_epi32(_mm512_mul_ps(_mm512_cvtepi32_ps(dd_v), avg_qspan_v), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

        // gap_cost = (dd * 0.01*avg_qspan) + (log_dd>>1 + 0.499)
        __m512i gap_cost_v = (_mm512_add_epi32(cost, log_dd_v));

        return gap_cost_v;
    }
#elif __AVX2__
//#endif
    inline __m256i get_gap_cost_vectorized_int32(__m256i dd_v, float avg_qspan, float gap_scale) {

        /*
        //Vectorized log2
        uint32_t base = 31;
        __m256i vbase = _mm256_set1_epi32(base);
        __mmask8 msk = 0xFF;
        	__m256i vout  = _mm256_mask_lzcnt_epi32(dd_v, msk, dd_v);
        	//int res = base - temp[0];
        	__m256i r_v = _mm256_sub_epi32(vbase, vout);
        */
        uint32_t lg[8], dd_array[8];
        _mm256_storeu_si256((__m256i *)dd_array, dd_v);
        for (int i = 0; i < 8; i++) {
            lg[i] = ilog2_32_dp_lib(dd_array[i]);
        }
        __m256i r_v = _mm256_loadu_si256((__m256i *)lg);

        // log_dd = dd?ilog2:0; log_dd>>1
        //__mmask8 neg_mask = _mm256_cmpneq_epi32_mask(dd_v, zero_avx2_v);
        //__m256i log_dd_v = _mm256_srli_epi32(_mm256_maskz_or_epi32(neg_mask, r_v, zero_avx2_v), 1);
        __m256i zero_avx2_v = _mm256_setzero_si256();
        __m256i neg_mask = _mm256_cmpeq_epi32(dd_v, zero_avx2_v);
        __m256i log_dd_v = _mm256_srli_epi32(_mm256_andnot_si256(neg_mask, r_v), 1);

        //dd * 0.01*avg_qspan
        float avg_qspan_val = 0.01 * avg_qspan;
        __m256 avg_qspan_v = _mm256_set1_ps(avg_qspan_val);
        __m256i cost = _mm256_cvtps_epi32(_mm256_round_ps((_mm256_mul_ps(_mm256_cvtepi32_ps(dd_v), avg_qspan_v)), _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC));

        // gap_cost = (dd * 0.01*avg_qspan) + (log_dd>>1 + 0.499)
        __m256i gap_cost_v = (_mm256_add_epi32(cost, log_dd_v));

        return gap_cost_v;
    }
#endif

static void chain_dp(call_t *a, return_t *ret) {
    constexpr float gap_scale = 1.0f;
    // constexpr int max_skip = 25;
    // constexpr int is_cdna = 0;
    constexpr int max_iter = 5000;

    const auto max_dist_x = a->max_dist_x;
    const auto max_dist_y = a->max_dist_y;
    const auto bw = a->bw;

    const auto avg_qspan = a->avg_qspan;
    const float avg_qspan001 = 0.01f * avg_qspan;

    // const auto n_segs = a->n_segs;
    const auto n = a->n;

    auto *anchors_x = a->anchors_x.data() + 32;
    auto *anchors_x32 = a->anchors_x32.data() + 32;

    auto *anchors_y = a->anchors_y.data() + 32;
    auto *anchors_y32 = a->anchors_y32.data() + 32;

    auto *q_spans = a->q_spans.data() + 32;

    ret->n = n;

    // Some extra space for vectorization with intrinsics.
    ret->scores.resize(n + 64);
    ret->parents.resize(n + 64);
    ret->targets.resize(n + 64);
    ret->peak_scores.resize(n + 64);

    // Add padding before and after the data.
    auto *scores = ret->scores.data() + 32;
    auto *parents = ret->parents.data() + 32;
    auto *targets = ret->targets.data() + 32;
    auto *peak_scores = ret->peak_scores.data() + 32;

    int32_t st = 0;

#ifdef __AVX512BW__
    #pragma message("Using AVX512 version")
    //Vector code with SoA function parameters 32-bit number representation - avx512
    
    __m512i zero_v = _mm512_setzero_si512();

    int32_t dr = max_dist_x;
    int32_t dq = max_dist_y;

    __m512i dr_v = _mm512_set1_epi32((int64_t)dr);
    __m512i dq_v = _mm512_set1_epi32((int64_t)dq);
    __m512i j_idx_base = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    int32_t maxfVector_v[16];
    int32_t maxjVector_v[16];
    __m512i neg_one_v = _mm512_set1_epi32((int32_t) -1);

    for (int i = 0; i < n; i++) {


        int32_t max_j = -1;
        int32_t max_f = q_spans[i];



        //uint64_t ri = anchors_x[i];
        while (st < i && !(anchors_x[i] - anchors_x[st] <= (uint32_t)dr)) {
            ++st;    //predecessor's position is too far
        }

        // TODO: Minimap specific max_iter parameter
        if (i - st > max_iter) {
            st = i - max_iter;    //predecessor's index is too far
        }


        int j = i - 1;
        if (!(j - st <= 5)) {
            //broadcast ri and qi
            __m512i ri_v = _mm512_set1_epi32(anchors_x32[i]);
            __m512i qi_v = _mm512_set1_epi32(anchors_y32[i]);


            __m512i maxj_v = neg_one_v;
            __m512i maxf_v = _mm512_set1_epi32((int32_t)q_spans[i]);
            __m512i li_v = maxf_v;

            // 16-way vectorized
            //_mm_prefetch((const char *)(&anchors_x32[j - 30]), _MM_HINT_T2);
            //_mm_prefetch((const char *)(&anchors_y32[j - 30]), _MM_HINT_T2);

            for (j = i - 1; (j - 15) >= st; j = j - 16) {

                _mm_prefetch((const char *)(&anchors_x32[j - 60]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&anchors_y32[j - 60]), _MM_HINT_T0);

                uint32_t *rj, *qj;
                rj = &anchors_x32[j - 15];
                qj = &anchors_y32[j - 15];


                // Load rj and qj
                __m512i rj_v = _mm512_loadu_si512(rj);
                __m512i qj_v = _mm512_loadu_si512(qj);


                __m512i ddr_v = _mm512_sub_epi32(ri_v, rj_v);
                __m512i ddq_v = _mm512_sub_epi32(qi_v, qj_v);

                //TODO: Minimap2 specific continue condition
                __m512i dd_v = _mm512_abs_epi32(_mm512_sub_epi32(ddr_v, ddq_v));
                __m512i bw_v = _mm512_set1_epi32((int32_t) bw);
                __mmask16 bw_gt = _mm512_cmpgt_epi32_mask(dd_v, bw_v);
                __mmask16 mask_eq = _mm512_cmpeq_epi32_mask(ddr_v, zero_v);
                __mmask16 mask_leq = _mm512_cmple_epi32_mask(ddq_v, zero_v);
                __mmask16 mask_gt1 = _mm512_cmpgt_epi32_mask(ddq_v, dq_v);
                __mmask16 mask_gt2 = _mm512_cmpgt_epi32_mask(ddq_v, dr_v);

                __mmask16 loopContinueMask = ~(bw_gt | mask_eq | mask_leq | mask_gt1 | mask_gt2);

                // Load scores[j-8, j]
                __m512i fj_v = _mm512_loadu_si512(&scores[j - 15]);

                //Vectorized gap cost function
                __m512i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

                //---------------- Inline get_overlap_cost function -------------------
                __m512i min1 = _mm512_min_epi32(ddr_v, ddq_v);
                __m512i oc_v = _mm512_min_epi32(li_v, min1);
                //----------------------------------------------------------------------
                __m512i f_plus_oc_v = _mm512_add_epi32(fj_v, oc_v);
                __m512i sc_v = _mm512_maskz_sub_epi32(loopContinueMask, f_plus_oc_v, gc_v);

                // Update Maxf and Maxj
                __mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
                __m512i j_idx_v = _mm512_add_epi32(j_idx_base, _mm512_set1_epi32(j - 15));

                maxf_v = _mm512_max_epi32(sc_v, maxf_v);
                maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);
            }


            if (j >= st) {
                uint32_t *rj, *qj;
                rj = &anchors_x32[j - 15];
                qj = &anchors_y32[j - 15];


                // Load rj and qj
                __m512i rj_v = _mm512_loadu_si512(rj);
                __m512i qj_v = _mm512_loadu_si512(qj);


                __m512i ddr_v = _mm512_sub_epi32(ri_v, rj_v);
                __m512i ddq_v = _mm512_sub_epi32(qi_v, qj_v);

                //TODO: Minimap2 specific continue condition
                __m512i dd_v = _mm512_abs_epi32(_mm512_sub_epi32(ddr_v, ddq_v));
                __m512i bw_v = _mm512_set1_epi32((int32_t) bw);
                __mmask16 bw_gt = _mm512_cmpgt_epi32_mask(dd_v, bw_v);
                __mmask16 mask_eq = _mm512_cmpeq_epi32_mask(ddr_v, zero_v);
                __mmask16 mask_leq = _mm512_cmple_epi32_mask(ddq_v, zero_v);
                __mmask16 mask_gt1 = _mm512_cmpgt_epi32_mask(ddq_v, dq_v);
                __mmask16 mask_gt2 = _mm512_cmpgt_epi32_mask(ddq_v, dr_v);


                __mmask16 loopContinueMask = ~(bw_gt | mask_eq | mask_leq | mask_gt1 | mask_gt2);


                //Last partial vector processing mask - To enable, change loop condition to  -> j >= st

                int shift = st - (j - 15);

                loopContinueMask = loopContinueMask >> (shift);
                loopContinueMask = loopContinueMask << (shift);
                if (loopContinueMask != 0x00) {



                    // Load scores[j-8, j]
                    __m512i fj_v = _mm512_loadu_si512(&scores[j - 15]);

                    //Vectorized gap cost function
                    __m512i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

                    //---------------- Inline get_overlap_cost function -------------------
                    __m512i min1 = _mm512_min_epi32(ddr_v, ddq_v);
                    __m512i oc_v = _mm512_min_epi32(li_v, min1);
                    //----------------------------------------------------------------------
                    __m512i f_plus_oc_v = _mm512_add_epi32(fj_v, oc_v);
                    __m512i sc_v = _mm512_maskz_sub_epi32(loopContinueMask, f_plus_oc_v, gc_v);


                    // Update Maxf and Maxj
                    __mmask16 mask_max_sc = _mm512_cmpgt_epi32_mask(sc_v, maxf_v);
                    __m512i j_idx_v = _mm512_add_epi32(j_idx_base, _mm512_set1_epi32(j - 15));

                    maxf_v = _mm512_max_epi32(sc_v, maxf_v);
                    maxj_v = _mm512_mask_blend_epi32(mask_max_sc, maxj_v, j_idx_v);

                }
            }


            _mm512_store_epi32(maxfVector_v, maxf_v);
            _mm512_store_epi32(maxjVector_v, maxj_v);

            for (int iter = 15; iter >= 0; iter--) {
                if (maxfVector_v[iter] > max_f) {
                    max_f = maxfVector_v[iter];
                    max_j = maxjVector_v[iter];
                }
                else if (maxfVector_v[iter] == max_f) {
                    max_j = std::max(max_j, maxjVector_v[iter]);
                    if ((uint32_t)max_f == q_spans[i]) {
                        max_j = -1;
                    }
                }
            }


        }
        else {
            int32_t ri = anchors_x32[i], qi = anchors_y32[i];
            for (; j >= st; j--) {

                int32_t ddr, ddq;

                int32_t rj = anchors_x32[j],
                        qj = anchors_y32[j];
                ddr = ri - rj;
                ddq = qi - qj;

                if (abs(ddr - ddq) > bw) {
                    continue;
                }
                if (ddr == 0 || ddq <= 0) {
                    continue;
                }

                if (ddq > dq || ddq > dr) {
                    continue;
                }


                int32_t oc = 0;

                int32_t score = scores[j];//q_spans[i];
                oc = ddr < ddq ? ddr : ddq;
                oc = oc < (int32_t)q_spans[i] ? oc : q_spans[i];
                score += oc;

                int32_t dr = ddr;
                int32_t dq = ddq;
                int32_t dd = abs(dr - dq);//dr > dq? dr - dq : dq - dr; //dd = |dr-dq|;
                int32_t log_dd = dd ? ilog2_32_dp_lib(dd) : 0;
                int32_t gap_cost = 0;

                gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd >> 1) + 0.00; //TODO: Only multiplication should be casted to (int)
                score -= gap_cost;

                if (score > max_f) {
                    max_f = score;
                    max_j = j;
                }

            }

        }

        scores[i] = max_f;
        parents[i] = max_j;
        peak_scores[i] = max_j >= 0 && peak_scores[max_j] > max_f ? peak_scores[max_j] : max_f; // v[] keeps the peak score up to i; scores[] is the score ending at i, not always the peak
    }
#elif __AVX2__
    #pragma message("Using AVX2 version")

    __m256i zero_avx2_v = _mm256_setzero_si256();

    int32_t dr = max_dist_x;
    int32_t dq = max_dist_y;

    __m256i dr_v = _mm256_set1_epi32(dr);
    __m256i dq_v = _mm256_set1_epi32(dq);
    __m256i j_idx_base = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    int32_t maxfVector_v[8];
    int32_t maxjVector_v[8];
    __m256i neg_one_v = _mm256_set1_epi32((int32_t) -1);

    for (int i = 0; i < n; i++) {


        int32_t max_j = -1;
        int32_t max_f = q_spans[i];



        //uint64_t ri = anchors_x[i];
        while (st < i && !(anchors_x[i] - anchors_x[st] <= dr)) {
            ++st;    //predecessor's position is too far
        }

        // TODO: Minimap specific max_iter parameter
        if (i - st > max_iter) {
            st = i - max_iter;    //predecessor's index is too far
        }


        int j = i - 1;
        if (!(j - st <= 5)) {
            //broadcast ri and qi
            __m256i ri_v = _mm256_set1_epi32(anchors_x32[i]);
            __m256i qi_v = _mm256_set1_epi32(anchors_y32[i]);


            __m256i maxj_v = neg_one_v;
            __m256i maxf_v = _mm256_set1_epi32((int32_t)q_spans[i]);
            __m256i li_v = maxf_v;

            // 8-way vectorized
            for (j = i - 1; (j - 7) >= st; j = j - 8) {

                _mm_prefetch((const char *)(&anchors_x32[j - 60]), _MM_HINT_T0);
                _mm_prefetch((const char *)(&anchors_y32[j - 60]), _MM_HINT_T0);

                uint32_t *rj, *qj;
                int j_stride = j - 7;
                rj = &anchors_x32[j_stride];
                qj = &anchors_y32[j_stride];


                // Load rj and qj
                __m256i rj_v = _mm256_loadu_si256((__m256i *) rj);
                __m256i qj_v = _mm256_loadu_si256((__m256i *) qj);


                __m256i ddr_v = _mm256_sub_epi32(ri_v, rj_v);
                __m256i ddq_v = _mm256_sub_epi32(qi_v, qj_v);

                //TODO: Minimap2 specific continue condition
                __m256i dd_v = _mm256_abs_epi32(_mm256_sub_epi32(ddr_v, ddq_v));
                __m256i bw_v = _mm256_set1_epi32((int32_t) bw);

                /*
                    __mmask8 bw_gt = _mm256_cmpgt_epi32_mask(dd_v, bw_v);
                    __mmask8 mask_eq = _mm256_cmpeq_epi32_mask(ddr_v, zero_avx2_v);
                    __mmask8 mask_leq1 = _mm256_cmpgt_epi32_mask(zero_avx2_v, ddq_v);
                    __mmask8 mask_leq2 = _mm256_cmpeq_epi32_mask(ddq_v, zero_avx2_v);
                    __mmask8 mask_gt1 = _mm256_cmpgt_epi32_mask(ddq_v, dq_v);
                    __mmask8 mask_gt2 = _mm256_cmpgt_epi32_mask(ddq_v, dr_v);


                    __mmask8 loopContinueMask = ~(bw_gt | mask_eq | (mask_leq1 | mask_leq2 ) | mask_gt1 | mask_gt2);
                */

                __m256i bw_gt = _mm256_cmpgt_epi32(dd_v, bw_v);
                __m256i mask_eq = _mm256_cmpeq_epi32(ddr_v, zero_avx2_v);
                __m256i mask_leq1 = _mm256_cmpgt_epi32(zero_avx2_v, ddq_v);
                __m256i mask_leq2 = _mm256_cmpeq_epi32(ddq_v, zero_avx2_v);
                __m256i mask_gt1 = _mm256_cmpgt_epi32(ddq_v, dq_v);
                __m256i mask_gt2 = _mm256_cmpgt_epi32(ddq_v, dr_v);

                __m256i tmp1 = _mm256_or_si256(bw_gt, mask_eq);
                __m256i tmp2 = _mm256_or_si256(mask_leq1, mask_leq2);
                __m256i tmp3 = _mm256_or_si256(mask_gt1, mask_gt2);

                __m256i loopContinueMask = _mm256_or_si256(_mm256_or_si256(tmp1, tmp2), tmp3);

                // Load scores[j-8, j]
                __m256i fj_v = _mm256_loadu_si256((__m256i *) &scores[j_stride]);

                //Vectorized gap cost function
                __m256i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

                //---------------- Inline get_overlap_cost function -------------------
                __m256i min1 = _mm256_min_epi32(ddr_v, ddq_v);
                __m256i oc_v = _mm256_min_epi32(li_v, min1);
                //----------------------------------------------------------------------
                __m256i f_plus_oc_v = _mm256_add_epi32(fj_v, oc_v);
                //__m256i sc_v = _mm256_maskz_sub_epi32(loopContinueMask,f_plus_oc_v, gc_v);
                __m256i sc_v = _mm256_andnot_si256(loopContinueMask, _mm256_sub_epi32(f_plus_oc_v, gc_v));

                // Update Maxf and Maxj
                __m256i mask_max_sc = _mm256_cmpgt_epi32(sc_v, maxf_v);
                __m256i j_idx_v = _mm256_add_epi32(j_idx_base, _mm256_set1_epi32(j_stride));

                maxf_v = _mm256_max_epi32(sc_v, maxf_v);
                maxj_v = _mm256_or_si256(_mm256_andnot_si256(mask_max_sc, maxj_v), _mm256_and_si256(mask_max_sc, j_idx_v));
            }


            if (j >= st) {
                uint32_t *rj, *qj;
                int j_stride = j - 7;
                rj = &anchors_x32[j_stride];
                qj = &anchors_y32[j_stride];


                // Load rj and qj
                __m256i rj_v = _mm256_loadu_si256((__m256i *) rj);
                __m256i qj_v = _mm256_loadu_si256((__m256i *) qj);


                __m256i ddr_v = _mm256_sub_epi32(ri_v, rj_v);
                __m256i ddq_v = _mm256_sub_epi32(qi_v, qj_v);

                //TODO: Minimap2 specific continue condition
                __m256i dd_v = _mm256_abs_epi32(_mm256_sub_epi32(ddr_v, ddq_v));
                __m256i bw_v = _mm256_set1_epi32((int32_t) bw);
                /*
                                __mmask8 bw_gt = _mm256_cmpgt_epi32_mask(dd_v, bw_v);
                                __mmask8 mask_eq = _mm256_cmpeq_epi32_mask(ddr_v, zero_avx2_v);
                                __mmask8 mask_leq = _mm256_cmple_epi32_mask(ddq_v, zero_avx2_v);
                                __mmask8 mask_gt1 = _mm256_cmpgt_epi32_mask(ddq_v, dq_v);
                                __mmask8 mask_gt2 = _mm256_cmpgt_epi32_mask(ddq_v, dr_v);


                                __mmask8 loopContinueMask = ~(bw_gt | mask_eq | mask_leq | mask_gt1 | mask_gt2);
                */
                __m256i bw_gt = _mm256_cmpgt_epi32(dd_v, bw_v);
                __m256i mask_eq = _mm256_cmpeq_epi32(ddr_v, zero_avx2_v);
                __m256i mask_leq1 = _mm256_cmpgt_epi32(zero_avx2_v, ddq_v);
                __m256i mask_leq2 = _mm256_cmpeq_epi32(ddq_v, zero_avx2_v);
                __m256i mask_gt1 = _mm256_cmpgt_epi32(ddq_v, dq_v);
                __m256i mask_gt2 = _mm256_cmpgt_epi32(ddq_v, dr_v);

                __m256i tmp1 = _mm256_or_si256(bw_gt, mask_eq);
                __m256i tmp2 = _mm256_or_si256(mask_leq1, mask_leq2);
                __m256i tmp3 = _mm256_or_si256(mask_gt1, mask_gt2);

                __m256i loopContinueMask = _mm256_or_si256(_mm256_or_si256(tmp1, tmp2), tmp3);

                //Last partial vector processing mask - To enable, change loop condition to  -> j >= st

                int shift = st - (j_stride);

                int32_t msk_ar[8];
                for (int it = 0; it < 8; it++) {
                    msk_ar[it] = (it < (shift)) ? 0xFFFFFFFF : 0;
                }
                loopContinueMask = _mm256_or_si256(loopContinueMask, _mm256_loadu_si256((__m256i *)msk_ar));
                //loopContinueMask = loopContinueMask>>(shift);
                //loopContinueMask = loopContinueMask<<(shift);
                //if(loopContinueMask != 0x0)
                {



                    // Load scores[j-8, j]
                    __m256i fj_v = _mm256_loadu_si256((__m256i *) &scores[j_stride]);

                    //Vectorized gap cost function
                    __m256i gc_v = get_gap_cost_vectorized_int32(dd_v, avg_qspan, gap_scale);

                    //---------------- Inline get_overlap_cost function -------------------
                    __m256i min1 = _mm256_min_epi32(ddr_v, ddq_v);
                    __m256i oc_v = _mm256_min_epi32(li_v, min1);
                    //----------------------------------------------------------------------
                    __m256i f_plus_oc_v = _mm256_add_epi32(fj_v, oc_v);
                    //__m256i sc_v = _mm256_maskz_sub_epi32(loopContinueMask,f_plus_oc_v, gc_v);
                    __m256i sc_v = _mm256_andnot_si256(loopContinueMask, _mm256_sub_epi32(f_plus_oc_v, gc_v));


                    // Update Maxf and Maxj
                    __m256i mask_max_sc = _mm256_cmpgt_epi32(sc_v, maxf_v);
                    __m256i j_idx_v = _mm256_add_epi32(j_idx_base, _mm256_set1_epi32(j_stride));

                    maxf_v = _mm256_max_epi32(sc_v, maxf_v);
                    maxj_v = _mm256_or_si256(_mm256_andnot_si256(mask_max_sc, maxj_v), _mm256_and_si256(mask_max_sc, j_idx_v));

                }
            }


            //_mm256_store_epi32(maxfVector_v, maxf_v);
            //_mm256_store_epi32(maxjVector_v, maxj_v);
            _mm256_store_si256((__m256i *) maxfVector_v, maxf_v);
            _mm256_store_si256((__m256i *) maxjVector_v, maxj_v);

            for (int iter = 7; iter >= 0; iter--) {
                if (maxfVector_v[iter] > max_f) {
                    max_f = maxfVector_v[iter];
                    max_j = maxjVector_v[iter];
                }
                else if (maxfVector_v[iter] == max_f) {
                    max_j = std::max(max_j, maxjVector_v[iter]);
                    if (max_f == q_spans[i]) {
                        max_j = -1;
                    }
                }
            }


        }
        else {
            int32_t ri = anchors_x32[i], qi = anchors_y32[i];
            for (; j >= st; j--) {

                int32_t ddr, ddq;

                int32_t rj = anchors_x32[j],
                        qj = anchors_y32[j];
                ddr = ri - rj;
                ddq = qi - qj;

                if (abs(ddr - ddq) > bw) {
                    continue;
                }
                if (ddr == 0 || ddq <= 0) {
                    continue;
                }

                if (ddq > dq || ddq > dr) {
                    continue;
                }


                int32_t oc = 0;
                int32_t lj = q_spans[j];
                int32_t ref_overlap = rj + lj - ri;
                int32_t query_overlap = qj + lj - qi;

                int32_t score = scores[j];//q_spans[i];
                oc = ddr < ddq ? ddr : ddq;
                oc = oc < q_spans[i] ? oc : q_spans[i];
                score += oc;

                int32_t dr = ddr;
                int32_t dq = ddq;
                int32_t dd = abs(dr - dq);//dr > dq? dr - dq : dq - dr; //dd = |dr-dq|;
                int32_t log_dd = dd ? ilog2_32_dp_lib(dd) : 0;
                int32_t gap_cost = 0;

                gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd >> 1) + 0.00; //TODO: Only multiplication should be casted to (int)
                score -= gap_cost;

                bool check = score > max_f;
                if (score > max_f) {
                    max_f = score;
                    max_j = j;
                }

            }

        }

        scores[i] = max_f;
        parents[i] = max_j;
        peak_scores[i] = max_j >= 0 && peak_scores[max_j] > max_f ? peak_scores[max_j] : max_f; // v[] keeps the peak score up to i; scores[] is the score ending at i, not always the peak
    }
#elif __ARM_FEATURE_SVE
    #pragma message("Using SVE version")

    // fill the score and backtrack arrays
    for (int32_t i = 0; i < n; ++i) {
        const int32_t ri_scalar = anchors_x32[i];
        svint32_t ri = svdup_n_s32(ri_scalar);
        const int32_t qi_scalar = static_cast<int32_t>(anchors_y32[i]);
        svint32_t qi = svdup_n_s32(qi_scalar);
        const int32_t q_spani = q_spans[i];

        int32_t max_j = -1;
        int32_t max_f = q_spani;

        while (st < i && !(anchors_x[i] - anchors_x[st] <= max_dist_x)) {
            ++st; //predecessor's position is too far
        }
        if (i - st > max_iter) {
            st = i - max_iter; //predecessor's index is too far
        }

        //for (int64_t j = i - 1; j >= st; --j) {
        svbool_t ptrue = svptrue_b32();
        for (int32_t j = i - 1; j >= st; j-=VL) {
            int32_t real_j = j-VL+1;
            svbool_t valid_elements = svnot_b_z(ptrue,svwhilelt_b32_s32(real_j,st));
            //const auto rj = anchors_x[j];
            svint32_t rj = svld1_s32(valid_elements,(int32_t*)&anchors_x32[real_j]);
            //const int32_t qj = static_cast<int32_t>(anchors_y[j]);
            svint32_t qj = svld1_s32(valid_elements,(int32_t*)&anchors_y32[real_j]);
            // qj = svextw_s64_x(valid_elements,qj);

            //const int64_t dr = ri - rj;
            svint32_t dr = svsub_s32_x(valid_elements,ri,rj);
            //const int32_t dq = qi - qj;
            svint32_t dq = svsub_s32_x(valid_elements,qi,qj);

            //const int32_t dd = std::abs(dr - dq);
            svint32_t dd = svabd_s32_x(valid_elements,dr,dq);

            //if ((dr == 0 || dq <= 0) ||
            //   (dq > max_dist_y || dq > max_dist_x) ||
            //   (dd > bw)) {
            //    continue;
            //}
            svbool_t skip_anchor = svcmpeq_n_s32(valid_elements,dr,0);
            valid_elements = svbic_b_z(valid_elements,valid_elements,skip_anchor);

            skip_anchor = svcmple_n_s32(valid_elements,dq,0);
            valid_elements = svbic_b_z(valid_elements,valid_elements,skip_anchor);

            skip_anchor = svcmpgt_n_s32(valid_elements,dq,max_dist_y);
            valid_elements = svbic_b_z(valid_elements,valid_elements,skip_anchor);

            skip_anchor = svcmpgt_n_s32(valid_elements,dq,max_dist_x);
            valid_elements = svbic_b_z(valid_elements,valid_elements,skip_anchor);

            skip_anchor = svcmpgt_n_s32(valid_elements,dd,bw);
            valid_elements = svbic_b_z(valid_elements,valid_elements,skip_anchor);

            if (!svptest_any(ptrue,valid_elements)) continue;

            //const int64_t dr_dq_min = (dr < dq) ? dr : dq;
            svint32_t dr_dq_min = svmin_s32_x(valid_elements,dr,dq);
            //const int32_t oc = (dr_dq_min < q_spani) ? dr_dq_min : q_spani;
            svint32_t oc = svmin_n_s32_x(valid_elements,dr_dq_min,q_spani);

            //const int32_t log_dd = (dd) ? ilog2_32(dd) : 0;
            svbool_t valid_log = svcmpne_n_s32(valid_elements,dd,0);
            svint32_t log_dd = svreinterpret_s32(svclz_s32_x(valid_elements,dd));
            log_dd = svsubr_n_s32_z(valid_log,log_dd,31);
            
            //const int32_t gap_cost = static_cast<int>(dd * 0.01f * avg_qspan) + (log_dd >> 1);
            /*
            svint64_t gap_cost = svmul_n_s64_x(valid_elements,dd,avg_qspan001);
            log_dd = svlsr_n_s64_x(valid_elements,log_dd,1);
            gap_cost = svadd_s64_x(valid_elements,gap_cost,log_dd);
            */
            svfloat32_t gap_cost_f = svcvt_f32_s32_x(valid_elements,dd);
            //gap_cost_f = svmul_n_f64_x(valid_elements,gap_cost_f,0.01f);
            gap_cost_f = svmul_n_f32_x(valid_elements,gap_cost_f,avg_qspan001);
            svint32_t gap_cost = svcvt_s32_f32_x(valid_elements,gap_cost_f);
            log_dd = svreinterpret_s32(svlsr_n_u32_x(valid_elements,svreinterpret_u32(log_dd),1));
            gap_cost = svadd_s32_x(valid_elements,gap_cost,log_dd);

            //const int32_t score = scores[j] + oc - gap_cost;
            svint32_t score = svld1_s32(valid_elements,&scores[real_j]);
            score = svadd_s32_x(valid_elements,score,oc);
            score = svsub_s32_x(valid_elements,score,gap_cost);

            // TODO: CAN'T VECTORIZE THIS. THE COMPILER IS ONLY ABLE TO PERFORM
            // ONE REDUCTION.
            // Ideas:
            //     1. j to int32_t and: uint64_t max = (j << 32) & score.  
            //        max = ((int32_t)(max & 0xffffffff) > score) ? max = (j << 32) & score : max;
            // Reduction

            // TODO: reduction on vector register, no scalar register
            //if (score > max_f) {
            //    max_f = score;
            //    max_j = j;
            //}

            int32_t max_local = svmaxv_s32(valid_elements,score);
            if (max_local > max_f) {
                max_f = max_local;
                // WARNING
                svint32_t index = svindex_s32(real_j,1);
                svbool_t max_index = svcmpeq_n_s32(valid_elements,score,max_local);
                max_j = svlastb_s32(max_index,index);
            }
        }
        scores[i] = max_f;
        parents[i] = max_j;
        //if (max_f == 36 && max_j == 18821) {printf("ee\n");exit(0);}
        peak_scores[i] = max_j >= 0 && peak_scores[max_j] > max_f ? peak_scores[max_j] : max_f;
    }
#else // SCALAR VERSION
    #pragma message("Using SCALAR version")
    for (int32_t i = 0; i < n; ++i) {
        const uint32_t ri = anchors_x32[i];
        const int32_t qi = static_cast<int32_t>(anchors_y32[i]);
        const int32_t q_spani = q_spans[i];

        int32_t max_j = -1;
        int32_t max_f = q_spani;

        while (st < i && !(anchors_x[i] - anchors_x[st] <= max_dist_x)) {
            ++st; //predecessor's position is too far
        }
        if (i - st > max_iter) {
            st = i - max_iter; //predecessor's index is too far
        }

        // TODO: Iterate forward to vectorize the loop.
        // for (int64_t j_inv = st; j_inv < i; ++j_inv) {
        //     const int64_t j = (i - 1) - j_inv + st;
        for (int32_t j = i - 1; j >= st; --j) {
            const uint32_t rj = anchors_x32[j];
            const uint32_t qj = anchors_y32[j];

            const int32_t dr = ri - rj;
            const int32_t dq = qi - qj;

            const int32_t dd = std::abs(dr - dq);

            if ((dr == 0 || dq <= 0) ||
               (dq > max_dist_y || dq > max_dist_x) ||
               (dd > bw)) {
                continue;
            }
            // Can not vectorize "continue". Use a mask instead.
            // const bool skip = ((dr == 0 || dq <= 0) ||
            //                    (dq > max_dist_y || dq > max_dist_x) ||
            //                    (dd > bw));

            const int32_t dr_dq_min = (dr < dq) ? dr : dq;
            const int32_t oc = (dr_dq_min < q_spani) ? dr_dq_min : q_spani;

            // TODO: CAN'T VECTORIZE __builtin_clz 
            const int32_t log_dd = (dd) ? ilog2_32(dd) : 0;
            const int32_t gap_cost = static_cast<int>(dd * 0.01f * avg_qspan) + (log_dd >> 1);

            const int32_t score = scores[j] + oc - gap_cost;

            // TODO: CAN'T VECTORIZE THIS. THE COMPILER IS ONLY ABLE TO PERFORM
            // ONE REDUCTION.
            // Ideas:
            //     1. j to int32_t and: uint64_t max = (j << 32) & score.  
            //        max = ((int32_t)(max & 0xffffffff) > score) ? max = (j << 32) & score : max;
            // Reduction
            if (score > max_f) {
                max_f = score;
                max_j = j;
            }
        }
        scores[i] = max_f;
        parents[i] = max_j;
        peak_scores[i] = max_j >= 0 && peak_scores[max_j] > max_f ? peak_scores[max_j] : max_f;
    }
#endif
}

void host_chain_kernel(std::vector<call_t> &args, std::vector<return_t> &rets, int numThreads) {
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp for schedule(dynamic)
        for (size_t batch = 0; batch < args.size(); batch++) {
            call_t *arg = &args[batch];
            return_t *ret = &rets[batch];
            // fprintf(stderr, "%lld\t%f\t%d\t%d\t%d\t%d\n", arg->n, arg->avg_qspan, arg->max_dist_x, arg->max_dist_y, arg->bw, arg->n_segs);
            chain_dp(arg, ret);
        }
    }
}
