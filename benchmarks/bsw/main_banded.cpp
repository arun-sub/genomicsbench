/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Vasimuddin Md, Sanchit Misra, Intel Corporation, Heng Li.

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

Authors: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>.
*****************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <limits>
#include <unistd.h>
#include <omp.h>
#include <fstream>
#include "utils.h"
#include "bandedSWA.h"

// #define VTUNE_ANALYSIS 1

#define CLMUL 8
#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

#define DEFAULT_MATCH 1
#define DEFAULT_MISMATCH 4
#define DEFAULT_OPEN 6
#define DEFAULT_EXTEND 1
#define DEFAULT_AMBIG -1

#undef MAX_SEQ_LEN_REF
#define MAX_SEQ_LEN_REF 2048
#undef MAX_SEQ_LEN_QER
#define MAX_SEQ_LEN_QER 256

// #define MAX_NUM_PAIRS 1000
// #define MATRIX_MIN_CUTOFF -100000000
// #define LOW_INIT_VALUE (INT32_MIN/2)
// #define AMBIG 52
double freq = 2.6*1e9;
int32_t w_match, w_mismatch, w_open, w_extend, w_ambig, numThreads = 1, batchSize = 0;
uint64_t SW_cells;
char *pairFileName;
FILE *pairFile;
int8_t h0 = 0;
double clock_freq;
uint64_t prof[10][112], data, SW_cells2;

void bwa_fill_scmat(int a, int b, int ambig, int8_t mat[25]) {
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = ambig; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = ambig;
}


void parseCmdLine(int argc, char *argv[])
{
	int i;
	w_match = DEFAULT_MATCH;
	w_mismatch = DEFAULT_MISMATCH;
	w_open = DEFAULT_OPEN;
	w_extend = DEFAULT_EXTEND;
	w_ambig = DEFAULT_AMBIG;

	int pairFlag = 0;
	for(i = 1; i < argc; i+=2) {
		if(strcmp(argv[i], "-match") == 0) {
			w_match = atoi(argv[i + 1]);
		}
		if(strcmp(argv[i], "-mismatch") == 0) { //penalty, +ve number
			w_mismatch = atoi(argv[i + 1]);
		}
		if(strcmp(argv[i], "-ambig") == 0) {
			w_ambig = atoi(argv[i + 1]);
		}

		if(strcmp(argv[i], "-gapo") == 0) {
			w_open = atoi(argv[i + 1]);
		}
		if(strcmp(argv[i], "-gape") == 0) {
			w_extend = atoi(argv[i + 1]);
		}
		if(strcmp(argv[i], "-pairs") == 0) {
			pairFileName = argv[i + 1];
			pairFlag = 1;
		}
		if(strcmp(argv[i], "-h0") == 0) {
			h0 = atoi(argv[i + 1]);
		}
        if(strcmp(argv[i], "-t") == 0) {
			numThreads = atoi(argv[i + 1]);
		}
        if(strcmp(argv[i], "-b") == 0) {
			batchSize = atoi(argv[i + 1]);
		}
	}
	if(pairFlag == 0) {
		fprintf(stderr, "ERROR! pairFileName not specified.\n");
		exit(EXIT_FAILURE);
	}
}

// -------------------------------------------------------------------------
// INPUT FILE FORMAT
// -------------------------------------------------------------------------
// Line 1: Seed score
// Line 2: Reference string
// Line 3: Query string
// 
// E.g:
// 19
// 01230123
// 0123

void loadPairs(SeqPair *seqPairArray, uint8_t *seqBufRef, uint8_t* seqBufQer, size_t numPairs)
{
	size_t numPairsRead = 0;
	while (numPairsRead < numPairs) {
		int32_t h0 = 0;
		char temp[10];
		fgets(temp, 10, pairFile);
		sscanf(temp, "%d", &h0);
		if (!fgets((char *)(seqBufRef + numPairsRead * (int64_t)(MAX_SEQ_LEN_REF)), MAX_SEQ_LEN_REF, pairFile)) {
			printf("WARNING! fgets returned NULL in %s. Num Pairs : %d\n", pairFileName, numPairsRead);
			break;
        }
		if (!fgets((char *)(seqBufQer + numPairsRead * (int64_t)(MAX_SEQ_LEN_QER)), MAX_SEQ_LEN_QER, pairFile)) {
			printf("WARNING! Odd number of sequences in %s\n", pairFileName);
			break;
        }
		SeqPair sp;
		sp.id = numPairsRead;
		sp.len1 = strnlen((char *)(seqBufRef + numPairsRead * MAX_SEQ_LEN_REF), MAX_SEQ_LEN_REF) - 1;
		sp.len2 = strnlen((char *)(seqBufQer + numPairsRead * MAX_SEQ_LEN_QER), MAX_SEQ_LEN_QER) - 1;
        if (sp.len1 <= 0 || sp.len2 <= 0) {
            fprintf(stderr, "%d\n", numPairsRead);
        }
        assert(sp.len1 > 0);
        assert(sp.len2 > 0);
		sp.h0 = h0;
		uint8_t *seq1 = seqBufRef + numPairsRead * MAX_SEQ_LEN_REF;
		uint8_t *seq2 = seqBufQer + numPairsRead * MAX_SEQ_LEN_QER;
		sp.idr =  numPairsRead * MAX_SEQ_LEN_REF;
		sp.idq =  numPairsRead * MAX_SEQ_LEN_QER;
		for (int l = 0; l < sp.len1; l++) {
			seq1[l] -= 48;
        }
		for (int l = 0; l < sp.len2; l++) {
			seq2[l] -= 48;
        }
		sp.seqid = sp.regid = sp.score = sp.tle = sp.gtle = sp.qle = -1;
		sp.gscore = sp.max_off = -1;
		seqPairArray[numPairsRead] = sp;
		numPairsRead++;
		// SW_cells += (sp.len1 * sp.len2);
	}
}
// profiling stats
uint64_t find_stats(uint64_t *val, int nt, double &min, double &max, double &avg) {
	min = 1e10;
	max = 0;
	avg = 0;
	for (int i=0; i<nt; i++) {
		avg += val[i];
		if (max < val[i]) max = val[i];
		if (min > val[i]) min = val[i];
	}
	avg /= nt;

	return 1;
}

int main(int argc, char *argv[])
{
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
	if (argc < 3) {
		fprintf(stderr, "usage: bsw -pairs <InSeqFile> -t <threads> -b <batch_size>\n");
		exit(EXIT_FAILURE);
	}
	
	parseCmdLine(argc, argv);

	pairFile = fopen(pairFileName, "r");	
	if (pairFile == NULL) {
		fprintf(stderr, "Could not open file: %s\n", pairFileName);
		exit(EXIT_FAILURE);
	}

    const int bufSize = 1024 * 1024;
    char* buffer = (char*)malloc(bufSize * sizeof(char));
    size_t numLines = 0;
    size_t n;
    while (n = fread(buffer, sizeof(char), bufSize, pairFile)) {
        for (int i = 0; i < n; i++) {
            if (buffer[i] == '\n') {
                numLines++;
            }
        }
    }
    free(buffer);

    // Reset file pointer back to the beginning
    fseek(pairFile, 0L, SEEK_SET);

    size_t numPairs = numLines / 3;
    size_t roundNumPairs = ((numPairs + SIMD_WIDTH16 - 1) / SIMD_WIDTH16 ) * SIMD_WIDTH16;
    if (batchSize == 0) { 
        batchSize = roundNumPairs;
    }
    printf("Number of input pairs: %ld\n", numPairs);

    size_t memAlloc = roundNumPairs * (sizeof(SeqPair) + MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_SEQ_LEN_REF * sizeof(int8_t));
	printf("Allocating %.3f GB memory for input buffers...\n", (memAlloc * 1.0)/ (1024 * 1024 * 1024));
	SeqPair *seqPairArray = (SeqPair *)_mm_malloc(roundNumPairs * sizeof(SeqPair), 64);
    uint8_t *seqBufQer = (uint8_t*) _mm_malloc(MAX_SEQ_LEN_QER * roundNumPairs * sizeof(int8_t), 64);
    uint8_t *seqBufRef = (uint8_t*) _mm_malloc(MAX_SEQ_LEN_REF * roundNumPairs * sizeof(int8_t), 64);

    int8_t mat[25];
	bwa_fill_scmat(w_match, w_mismatch, w_ambig, mat);
	int zdrop = 100, w = 100, end_bonus = 5;

	// OutScore *outScoreArray = (OutScore *)_mm_malloc(MAX_NUM_PAIRS * sizeof(OutScore), 64);
	BandedPairWiseSW *bsw[numThreads];
    for (int i = 0; i < numThreads; i++) {
        bsw[i] = new BandedPairWiseSW(w_open, w_extend, w_open, w_extend,
                    zdrop, end_bonus, mat,
                    w_match, w_mismatch, 1);
    }

	int64_t startTick, totalTicks = 0, readTim = 0;
	
	uint64_t tim = __rdtsc();

    int64_t numPairsIndex = 0;
    for (int64_t i = 0; i < roundNumPairs; i += batchSize) {
        int nPairsBatch = (numPairs - i) >= batchSize ? batchSize : numPairs - i;
        loadPairs(seqPairArray + numPairsIndex, seqBufRef + numPairsIndex * MAX_SEQ_LEN_REF, seqBufQer + numPairsIndex * MAX_SEQ_LEN_QER, nPairsBatch);
        numPairsIndex += nPairsBatch;
    }
    readTim += __rdtsc() - tim;

    startTick = __rdtsc();
    
    #ifdef VTUNE_ANALYSIS
        __itt_resume();
    #endif
    int64_t workTicks[CLMUL * numThreads];
    memset(workTicks, 0, CLMUL * numThreads * sizeof(int64_t));
#pragma omp parallel num_threads(numThreads)
{
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic, 1) 
        for (int64_t i = 0; i < roundNumPairs; i += batchSize) {
            int nPairsBatch = (numPairs - i) >= batchSize ? batchSize : numPairs - i;
            int64_t st1 = __rdtsc();
            bsw[tid]->getScores16(seqPairArray + i, seqBufRef + i * MAX_SEQ_LEN_REF, seqBufQer + i * MAX_SEQ_LEN_QER, nPairsBatch, 1, w);
            int64_t et1 = __rdtsc();
            workTicks[CLMUL * tid] += (et1 - st1);
        }
    printf("%d] workTicks = %ld\n", tid, workTicks[CLMUL * tid]);  
}
    
    #ifdef VTUNE_ANALYSIS
        __itt_pause();
    #endif
    totalTicks += __rdtsc() - startTick;
    printf("Executed AVX2 vector code...\n");

	tim = __rdtsc();
	sleep(1);
	freq = __rdtsc() - tim;
	
	printf("Processor freq: %0.2lf MHz\n", freq/1e6);

	//int64_t myTicks = bsw->getTicks();
	printf("Read time = %0.2lf s\n", readTim/freq);
	printf("Overall SW cycles = %ld, %0.2lf s\n", totalTicks, totalTicks * 1.0 / freq);
	printf("Total Pairs processed: %d\n", numPairs);
    
    int64_t sumTicks = 0;
    int64_t maxTicks = 0;
    for(int i = 0; i < numThreads; i++)
    {
        sumTicks += workTicks[CLMUL * i];
        if(workTicks[CLMUL * i] > maxTicks) maxTicks = workTicks[CLMUL * i];
    }
    double avgTicks = (sumTicks * 1.0) / numThreads;
    printf("avgTicks = %lf, maxTicks = %ld, load imbalance = %lf\n", avgTicks, maxTicks, maxTicks/avgTicks);


	// printf("SW cells(T)  = %ld\n", SW_cells);
	// printf("SW cells(||)  = %ld\n", SW_cells2);
	// printf("SW GCUPS  = %lf\n", SW_cells * freq/1e9 / myTicks);

	// printf("More stats:\n");
	// double freq = 2.3*1e9;
	// double min, max, avg;
	// find_stats(prof[1], numThreads, min, max, avg);
	// printf("Time in pre-processing: %0.2lf (%0.2lf, %0.2lf)\n",
	// avg*1.0/freq, min*1.0/freq, max*1.0/freq);
	// find_stats(prof[0], numThreads, min, max, avg);
	// printf("Time spent in smithWaterman(): %0.2lf (%0.2lf, %0.2lf)\n",
	// 	avg*1.0/freq, min*1.0/freq, max*1.0/freq);

	// printf("\nDebugging info:\n");
	// printf("Time taken for DP loop: %0.2lf\n", prof[DP][0]*1.0/freq);
	// printf("Time taken for DP loop upper part: %0.2lf\n", prof[DP3][0]*1.0/freq);	
	// printf("Time taken for DP inner loop: %0.2lf\n", prof[DP1][0]*1.0/freq);
	// printf("Time taken for DP loop lower part: %0.2lf\n", prof[DP2][0]*1.0/freq);	

	/**** free memory *****/
	_mm_free(seqPairArray);
	_mm_free(seqBufRef);
	_mm_free(seqBufQer); 
    // bsw->getTicks();
    for (int i = 0; i < numThreads; i++) {
        bsw[i]->getTicks();
	    delete bsw[i];
    }
	
	fclose(pairFile);
	return 1;
}
