#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
#include <string.h>
#include <math.h>
#include "ksort.h"
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "ksort.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)


#ifdef VTUNE_ANALYSIS
#include <ittnotify.h>
#endif

#define mem_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv_2, bwtintv_t, mem_lt)


typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = calloc(1, sizeof(bwtintv_v));
	return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

int main(int argc, char **argv) {

#ifdef VTUNE_ANALYSIS
	 __itt_pause();
#endif

	if (argc != 4) {
		printf("Need three arguments : <index_prefix> <query.fq> <min.seed.length>\n");
		return 1;
	}
	
	int min_seed_len = atoi(argv[3]);
	
	int32_t numReads = 0;
	gzFile fp = gzopen(argv[2], "r");
	if (fp == NULL) {
		fprintf(stderr, "[E::%s] failed to open file `%s'.\n", __func__, argv[2]);
		exit(1);
	}

	printf("Loading index ...\n");
	bwaidx_t* idx = bwa_idx_load(argv[1], BWA_IDX_ALL);

	kseq_t *ks = kseq_init(fp);

	while (kseq_read(ks) >= 0) { // read one sequence
		int i, k, x = 0, old_n;
		int start_width = 1;
		int split_len = (int)(min_seed_len * 1.5 + .499);
		smem_aux_t* a = smem_aux_init();
		a->mem.n = 0;
		int len = ks->seq.l;
		char* nt_seq = ks->seq.s;
		for (i = 0; i < len; ++i) {
			nt_seq[i] = nt_seq[i] < 4? nt_seq[i] : nst_nt4_table[(int)nt_seq[i]];
		}
		uint8_t* seq = (uint8_t*) nt_seq;
		
		printf("Processing read %s\n", ks->name.s);
		
		// first pass: find all SMEMs
		while (x < len) {
			if (seq[x] < 4) {
				x = bwt_smem1(idx->bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
				for (i = 0; i < a->mem1.n; ++i) {
					bwtintv_t *p = &a->mem1.a[i];
					int slen = (uint32_t)p->info - (p->info>>32); // seed length
					if (slen >= min_seed_len)
						kv_push(bwtintv_t, a->mem, *p);
				}
			} else ++x; 
		}
		// second pass: find MEMs inside a long SMEM
		old_n = a->mem.n;
		for (k = 0; k < old_n; ++k) {
			bwtintv_t *p = &a->mem.a[k];
			int start = p->info>>32, end = (int32_t)p->info;
			if (end - start < split_len || p->x[2] > 10) continue;
			bwt_smem1(idx->bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i)
				if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= min_seed_len)
					kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
		}
		// third pass: LAST-like
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				bwtintv_t m;
				x = bwt_seed_strategy1(idx->bwt, len, seq, x, min_seed_len, 20, &m);
				if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
			} else ++x;
		}

		//sort 
		ks_introsort(mem_intv_2, a->mem.n, a->mem.a);

		for (i = 0; i < a->mem.n; ++i) {
			bwtintv_t *p = &a->mem.a[i];
			for (k = 0; k < p->x[2]; ++k) {
				int qbeg = p->info>>32;
				int64_t rbeg = bwt_sa(idx->bwt, p->x[0] + k);
				int slen = (uint32_t)p->info - (p->info>>32);
				printf("[%d,%ld,%d]\n", qbeg, rbeg, slen);
			}
		}

		smem_aux_destroy(a);

	}

	kseq_destroy(ks);
	gzclose(fp);
	bwa_idx_destroy(idx);

	return 0;

}


