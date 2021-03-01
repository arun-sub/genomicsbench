#include <vector>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include "omp.h"
#include "host_kernel.h"
#include "common.h"
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"

static const char LogTable256[256] = {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
	-1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
	LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
	LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
};

static inline int ilog2_32(uint32_t v)
{
	uint32_t t, tt;
	if ((tt = v>>16)) return (t = tt>>8) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
	return (t = v>>8) ? 8 + LogTable256[t] : LogTable256[v];
}

const int BACKSEARCH = 65;
#define MM_SEED_SEG_SHIFT  48
#define MM_SEED_SEG_MASK   (0xffULL<<(MM_SEED_SEG_SHIFT))

void chain_dp(call_t* a, return_t* ret)
{

	// TODO: make sure this works when n has more than 32 bits
	int64_t i, j, st = 0;
	int is_cdna = 0;
    const float gap_scale = 1.0f;
    const int max_iter = 5000;
    const int max_skip = 25;
    int max_dist_x = a->max_dist_x, max_dist_y = a->max_dist_y, bw = a->bw;
    float avg_qspan = a->avg_qspan;
    int n_segs = a->n_segs; 
    int64_t n = a->n;
	ret->n = n;
	ret->scores.resize(n);
	ret->parents.resize(n);
    ret->targets.resize(n);
    ret->peak_scores.resize(n);

	// fill the score and backtrack arrays
	for (i = 0; i < n; ++i) {
		uint64_t ri = a->anchors[i].x;
		int64_t max_j = -1;
		int32_t qi = (int32_t)a->anchors[i].y, q_span = a->anchors[i].y>>32&0xff; // NB: only 8 bits of span is used!!!
		int32_t max_f = q_span, n_skip = 0, min_d;
		int32_t sidi = (a->anchors[i].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
		while (st < i && ri > a->anchors[st].x + max_dist_x) ++st;
		if (i - st > max_iter) st = i - max_iter;
		for (j = i - 1; j >= st; --j) {
			int64_t dr = ri - a->anchors[j].x;
			int32_t dq = qi - (int32_t)a->anchors[j].y, dd, sc, log_dd, gap_cost;
			int32_t sidj = (a->anchors[j].y & MM_SEED_SEG_MASK) >> MM_SEED_SEG_SHIFT;
			if ((sidi == sidj && dr == 0) || dq <= 0) continue; // don't skip if an anchor is used by multiple segments; see below
			if ((sidi == sidj && dq > max_dist_y) || dq > max_dist_x) continue;
			dd = dr > dq? dr - dq : dq - dr;
			if (sidi == sidj && dd > bw) continue;
			if (n_segs > 1 && !is_cdna && sidi == sidj && dr > max_dist_y) continue;
			min_d = dq < dr? dq : dr;
			sc = min_d > q_span? q_span : dq < dr? dq : dr;
			log_dd = dd? ilog2_32(dd) : 0;
			gap_cost = 0;
			if (is_cdna || sidi != sidj) {
				int c_log, c_lin;
				c_lin = (int)(dd * .01 * avg_qspan);
				c_log = log_dd;
				if (sidi != sidj && dr == 0) ++sc; // possibly due to overlapping paired ends; give a minor bonus
				else if (dr > dq || sidi != sidj) gap_cost = c_lin < c_log? c_lin : c_log;
				else gap_cost = c_lin + (c_log>>1);
			} else gap_cost = (int)(dd * .01 * avg_qspan) + (log_dd>>1);
			sc -= (int)((double)gap_cost * gap_scale + .499);
			sc += ret->scores[j];
			if (sc > max_f) {
				max_f = sc, max_j = j;
                if (n_skip > 0) --n_skip;
			} else if (ret->targets[j] == i) {
                if (++n_skip > max_skip) {
                    break;
                }
            }
            if (ret->parents[j] >= 0) ret->targets[ret->parents[j]] = i;
		}
		ret->scores[i] = max_f, ret->parents[i] = max_j;
        ret->peak_scores[i] = max_j >= 0 && ret->peak_scores[max_j] > max_f ? ret->peak_scores[max_j] : max_f;
	}
}

void host_chain_kernel(std::vector<call_t> &args, std::vector<return_t> &rets, int numThreads)
{
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp for schedule(dynamic)
            for (size_t batch = 0; batch < args.size(); batch++) {
                call_t* arg = &args[batch];
                return_t* ret = &rets[batch];
                // fprintf(stderr, "%lld\t%f\t%d\t%d\t%d\t%d\n", arg->n, arg->avg_qspan, arg->max_dist_x, arg->max_dist_y, arg->bw, arg->n_segs);
                chain_dp(arg, ret);
            }
    }
}
