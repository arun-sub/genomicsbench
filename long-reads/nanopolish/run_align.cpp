#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <ctime>
#include "nanopolish_common.h"
#include "nanopolish_squiggle_read.h"
#include "nanopolish_pore_model_set.h"
#include "nanopolish_methyltrain.h"
#include "nanopolish_extract.h"
#include "nanopolish_raw_loader.h"
#include "nanopolish_fast5_io.h"
#include "nanopolish_fast5_loader.h"
#include <stdio.h>

extern "C" {
#include "event_detection.h"
#include "scrappie_common.h"
}

#include <fast5.hpp>

#define event_kmer_to_band(ei, ki) (ei + 1) + (ki + 1)
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) (ki) - band_lower_left[bi].kmer_idx
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

#define move_down(curr_band) { curr_band.event_idx + 1, curr_band.kmer_idx }
#define move_right(curr_band) { curr_band.event_idx, curr_band.kmer_idx + 1 }

#define ALN_BANDWIDTH 100

#define BAND_ARRAY(r, c) ( bands[((r)*(ALN_BANDWIDTH)+(c))] )
#define TRACE_ARRAY(r, c) ( trace[((r)*(ALN_BANDWIDTH)+(c))] )

inline float log_probability_match_r9(SquiggleEvent event,
                                      SquiggleEvent event_0,
                                      SquiggleScalings scalings,
                                      PoreModelStateParams state)
{
    float log_inv_sqrt_2pi = log(0.3989422804014327);
    float level_tmp =event.mean;
    float time = event.start_time - event_0.start_time;
    float level = level_tmp - time * scalings.drift;
    float gp_mean = scalings.scale * state.level_mean + scalings.shift;
    float gp_stdv = state.level_stdv * scalings.var;
    float gp_log_stdv = state.level_log_stdv + scalings.log_var;
    float a = (level - gp_mean) / gp_stdv;
    float lp = log_inv_sqrt_2pi - gp_log_stdv + (-0.5f * a * a);
    return lp;
}

//todo : can make more efficient using bit encoding
static inline uint32_t get_rank(char base) {
    if (base == 'A') { //todo: do we neeed simple alpha?
        return 0;
    } else if (base == 'C') {
        return 1;
    } else if (base == 'G') {
        return 2;
    } else if (base == 'T') {
        return 3;
    } else {
        return 0;
    }
}

// return the lexicographic rank of the kmer amongst all strings of
// length k for this alphabet
static inline uint32_t get_kmer_rank(const char* str, uint32_t k) {
    //uint32_t p = 1;
    uint32_t r = 0;

    // from last base to first
    for (uint32_t i = 0; i < k; ++i) {
        //r += rank(str[k - i - 1]) * p;
        //p *= size();
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}


std::vector<AlignedPair> adaptive_banded_simple_event_align(size_t k, std::vector<SquiggleEvent> events, 
                                    SquiggleScalings scalings, std::vector<PoreModelStateParams> states,
                                    std::string& sequence)
{
    size_t strand_idx = 0;
    // size_t k = pore_model.k;
    // const Alphabet* alphabet = pore_model.pmalphabet;
    size_t n_events = events.size();
    size_t n_kmers = sequence.size() - k + 1;

    // backtrack markers
    const uint8_t FROM_D = 0;
    const uint8_t FROM_U = 1;
    const uint8_t FROM_L = 2;
 
    // qc
    double min_average_log_emission = -5.0;
    int max_gap_threshold = 50;

    // banding
    int bandwidth = ALN_BANDWIDTH;
    int half_bandwidth = bandwidth / 2;
 
    // transition penalties
    double events_per_kmer = (double)n_events / n_kmers;
    double p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);
 
    // dp matrix
    size_t n_rows = n_events + 1;
    size_t n_cols = n_kmers + 1;
    size_t n_bands = n_rows + n_cols;
 
    // Initialize

    // Precompute k-mer ranks to avoid doing this in the inner loop
    // std::vector<size_t> kmer_ranks(n_kmers);
    size_t* kmer_ranks = (size_t*)malloc(sizeof(size_t) * n_kmers);
    for(size_t i = 0; i < n_kmers; ++i) {
        // kmer_ranks[i] = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
        char* substring = &sequence[i];
        kmer_ranks[i] = get_kmer_rank(substring, k);
    }

    float* bands = (float*)malloc(sizeof(float) * n_bands * bandwidth);
    if(bands==NULL){
        fprintf(stderr,"Memory allocation failed at %s\n",__func__);
        exit(1);
    }
    uint8_t* trace = (uint8_t*)malloc(sizeof(uint8_t) * n_bands * bandwidth);
    if(trace==NULL){
        fprintf(stderr,"Memory allocation failed at %s\n",__func__);
        exit(1);
    }
    for (size_t i = 0; i < n_bands; i++) {
        for (int j = 0; j < bandwidth; j++) {
            BAND_ARRAY(i,j) = -INFINITY;
            TRACE_ARRAY(i,j) = 0;
        }
    }

    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two bands have their coordinates initialized, the rest are computed adaptively
    struct EventKmerPair
    {
        int event_idx;
        int kmer_idx;
    };

    std::vector<EventKmerPair> band_lower_left(n_bands);
 
    // initialize range of first two bands
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    // band 0: score zero in the central cell
    int start_cell_offset = band_kmer_to_offset(0, -1);
    assert(is_offset_valid(start_cell_offset));
    assert(band_event_to_offset(0, -1) == start_cell_offset);
    BAND_ARRAY(0,start_cell_offset) = 0.0f;

    // band 1: first event is trimmed
    int first_trim_offset = band_event_to_offset(1, 0);
    assert(kmer_at_offset(1, first_trim_offset) == -1);
    assert(is_offset_valid(first_trim_offset));
    BAND_ARRAY(1,first_trim_offset) = lp_trim;
    TRACE_ARRAY(1,first_trim_offset) = FROM_U;

    int fills = 0;
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1, first_trim_offset, 0, -1, BAND_ARRAY(1,first_trim_offset));
#endif

    // fill in remaining bands
    for(int band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        float ll = BAND_ARRAY(band_idx - 1,0);
        float ur = BAND_ARRAY(band_idx - 1,bandwidth - 1);
        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;

        bool right = false;
        if(ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if(right) {
            band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1]);
        } else {
            band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1]);
        }

        // If the trim state is within the band, fill it in here
        int trim_offset = band_kmer_to_offset(band_idx, -1);
        if(is_offset_valid(trim_offset)) {
            int event_idx = event_at_offset(band_idx, trim_offset);
            if(event_idx >= 0 && event_idx < n_events) {
                BAND_ARRAY(band_idx,trim_offset) = lp_trim * (event_idx + 1);
                TRACE_ARRAY(band_idx,trim_offset) = FROM_U;
            } else {
                BAND_ARRAY(band_idx,trim_offset) = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        // We restrict the inner loop to only these values
        int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int event_max_offset = band_event_to_offset(band_idx, -1);

        int min_offset = std::max(kmer_min_offset, event_min_offset);
        min_offset = std::max(min_offset, 0);

        int max_offset = std::min(kmer_max_offset, event_max_offset);
        max_offset = std::min(max_offset, bandwidth);
        for(int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = event_at_offset(band_idx, offset);
            int kmer_idx = kmer_at_offset(band_idx, offset);

            size_t kmer_rank = kmer_ranks[kmer_idx];

            int offset_up   = band_event_to_offset(band_idx - 1, event_idx - 1);
            int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
            int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
            // verify loop conditions
            assert(kmer_idx >= 0 && kmer_idx < n_kmers);
            assert(event_idx >= 0 && event_idx < n_events);
            assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
            assert(offset_up - offset_left == 1);
            assert(offset >= 0 && offset < bandwidth);
#endif

            float up   = is_offset_valid(offset_up)   ? BAND_ARRAY(band_idx - 1,offset_up)   : -INFINITY;
            float left = is_offset_valid(offset_left) ? BAND_ARRAY(band_idx - 1,offset_left) : -INFINITY;
            float diag = is_offset_valid(offset_diag) ? BAND_ARRAY(band_idx - 2,offset_diag) : -INFINITY;
            float lp_emission = log_probability_match_r9(events[event_idx], events[0], scalings, states[kmer_rank]);
            float score_d = diag + lp_step + lp_emission;
            float score_u = up + lp_stay + lp_emission;
            float score_l = left + lp_skip;

            float max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;

#ifdef DEBUG_ADAPTIVE
            fprintf(stderr, "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n", offset_up, offset_diag, offset_left);
            fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
            fprintf(stderr, "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", band_idx, offset, event_idx, kmer_idx, max_score, from, lp_emission);
#endif
            BAND_ARRAY(band_idx,offset) = max_score;
            TRACE_ARRAY(band_idx,offset) = from;
            fills += 1;
        }
    }

    //
    // Backtrack to compute alignment
    //
    double sum_emission = 0;
    double n_aligned_events = 0;
    std::vector<AlignedPair> out;

    float max_score = -INFINITY;
    int curr_event_idx = 0;
    int curr_kmer_idx = n_kmers -1;

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for(int event_idx = 0; event_idx < n_events; ++event_idx) {
        size_t band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
        size_t offset = band_event_to_offset(band_idx, event_idx);
        if(is_offset_valid(offset)) {
            float s = BAND_ARRAY(band_idx,offset)  + (n_events - event_idx) * lp_trim;
            if(s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }

#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, max_score);
#endif

    int curr_gap = 0;
    int max_gap = 0;
    while(curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        
        // emit alignment
        out.push_back({curr_kmer_idx, curr_event_idx});
#ifdef DEBUG_ADAPTIVE
        fprintf(stderr, "[adaback] ei: %d ki: %d\n", curr_event_idx, curr_kmer_idx);
#endif
        // qc stats
        // size_t kmer_rank = alphabet->kmer_rank(sequence.substr(curr_kmer_idx, k).c_str(), k);
        char* substring = &sequence[curr_kmer_idx];
        size_t kmer_rank = get_kmer_rank(substring, k);
        sum_emission += log_probability_match_r9(events[curr_event_idx], events[0], scalings, states[kmer_rank]);
        n_aligned_events += 1;

        size_t band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        size_t offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        uint8_t from = TRACE_ARRAY(band_idx,offset);
        if(from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if(from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = std::max(curr_gap, max_gap);
        }
    }
    std::reverse(out.begin(), out.end());
    
    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    
    bool failed = false;
    if(avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
        failed = true;
        out.clear();
    }

    free(bands);
    free(trace);

    //fprintf(stderr, "ada\t%s\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", read.read_name.substr(0, 6).c_str(), failed ? "FAILED" : "OK", events_per_kmer, sequence.size(), avg_log_emission, curr_event_idx, max_gap, fills);
    return out;
}

int main() {

    std::ifstream file1;
    file1.open("align_parameter_nanopolish", std::fstream::binary);

    size_t k;
    SquiggleScalings scalings;
    file1.read((char*)&k, sizeof(k));
    file1.read((char*)&scalings, sizeof(scalings));

    uint32_t n_events;
    std::vector<SquiggleEvent> events;
    SquiggleEvent event;
    file1.read((char*)&n_events, sizeof(n_events));
    for (uint32_t i = 0; i <= n_events; i++) {
        file1.read((char*)&event, sizeof(event));
        events.push_back(event);
    }

    uint32_t n_states;
    std::vector<PoreModelStateParams> states;
    PoreModelStateParams state;
    file1.read((char*)&n_states, sizeof(n_states));
    for (uint32_t i = 0; i <= n_states; i++) {
        file1.read((char*)&state, sizeof(state));
        states.push_back(state);
    }

    std::string sequence;
    uint32_t sequence_len;
    file1.read((char*)&sequence_len, sizeof(sequence_len));
    sequence.resize(sequence_len);
    file1.read((char*)sequence.c_str(), sequence_len);
    // file1 >> sequence;

    file1.close();

    std::vector<AlignedPair> event_alignment_new = adaptive_banded_simple_event_align(k, events, scalings, states, sequence);

    std::vector<SquiggleEvent>().swap(events);
    std::vector<PoreModelStateParams>().swap(states);

    int flag = 1;

    std::ifstream file2;
    
    file2.open("align_result_nanopolish", std::fstream::binary);
    AlignedPair event_alignment;
    uint32_t n_event_alignments;
    file2.read((char*)&n_event_alignments, sizeof(n_event_alignments));
    for (uint32_t i = 0; i <= n_event_alignments; i++) {
        file2.read((char*)&event_alignment, sizeof(event_alignment));
        if (!(event_alignment == event_alignment_new[i])) {
            flag = 0; 
            printf("%d original: %d %d new: %d %d\n", i, event_alignment.ref_pos, 
                   event_alignment.read_pos, event_alignment_new[i].ref_pos, event_alignment_new[i].read_pos);
        }
    }
    file2.close();

    if (flag) std::cout << "Correct!" << std::endl;
    else std::cout << "Wrong!" << std::endl;

    return 0;
}
