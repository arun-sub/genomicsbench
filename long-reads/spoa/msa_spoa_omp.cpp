/* To complile:
   g++ -O3 msa_spoa_omp.cpp -fopenmp -std=c++11 -I include/ -L build/lib/ -lspoa -o msa_spoa
  */

#include <getopt.h>
#include <stdio.h>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <assert.h>
#include <exception>
#include <getopt.h>
#include <omp.h>
#include "spoa/spoa.hpp"
#include "spoa/graph.hpp"
#include "spoa/alignment_engine.hpp"

#define VTUNE_ANALYSIS 1

#ifdef VTUNE_ANALYSIS
    #include <ittnotify.h>
#endif

using namespace std;
using Alignment = std::vector<std::pair<std::int32_t, std::int32_t>>;

// #define DEBUG_FILE_READ
// #define PRINT_OUTPUT

typedef struct {
    vector<string> seqs;
    string consensus_seq;
} Batch;

double get_realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec*1e6 + tp.tv_usec;
}

long peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

void readFile(ifstream& in_file, vector<Batch>& batches) {
    string seq;
    if (!in_file.eof()) {
        getline(in_file, seq); // header line
    }
    while (!in_file.eof()) { // loop across batches
        if (seq[1] == '0') {
            Batch b;
            while (!in_file.eof()) { // loop across sequences in batch
                getline(in_file, seq); // sequence line
                b.seqs.push_back(seq);
                getline(in_file, seq); // header line
                if (seq[1] == '0') {
                    batches.push_back(b);
                    break;
                }
            }
            if (in_file.eof()) {
                batches.push_back(b);
            }
        }
    }

#ifdef DEBUG_FILE_READ    
    for (int i = 0; i < batches.size(); i++) {
        cout << "Batch " << i << endl;
        for (int j = 0; j < batches[i].seqs.size(); j++) {
            cout << batches[i].seqs[j] << endl;
        }
    }
#endif

}


void help() {
    std::cout <<
        "\n"
        "usage: ./msa_spoa_omp [options ...]\n"
        "\n"
        "    options:\n"
        "        -m <int>\n"
        "            default: 2\n"
        "            score for matching bases\n"
        "        -x <int>\n"
        "            default: 4\n"
        "            penalty for mismatching bases\n"
        "        -o <int(,int)>\n"
        "            default: gap_open1 = 4, gap_open2 = 24\n"
        "            gap opening penalty (must be non-negative)\n"
        "        -e <int(,int)>\n"
        "            default: gap_ext1 = 2, gap_ext2 = 1\n"
        "            gap extension penalty (must be non-negative)\n"
        "        -s <file>\n"
        "            default: seq.fa\n"
        "            the input sequence set\n"
        "        -n <int>\n"
        "            default: 10\n"
        "            number of sequences in each set\n"
        "        -t <int>\n"
        "            default: 1\n"
        "            number of CPU threads\n"
        "        -h \n"
        "            prints the usage\n";
}

int main(int argc, char** argv) {
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    string seq_file = "seq.fa";

    std::uint8_t algorithm = 1;
    std::int8_t m = 2;
    std::int8_t x = -4;
    std::int8_t o1 = -4;
    std::int8_t e1 = -2;
    std::int8_t o2 = -24;
    std::int8_t e2 = -1;

    if (argc == 1) {
        help();
        exit(EXIT_FAILURE);
    }

    char opt, *s; int n_seqs = 0, numThreads = 1;
    while ((opt = getopt(argc, argv, "l:m:x:o:n:e:q:c:s:t:h")) != -1) {
        switch (opt) {
            case 'm': m = atoi(optarg); break;
            case 'x': x = 0-atoi(optarg); break;
            case 'o': o1 = 0-strtol(optarg, &s, 10); if (*s == ',') o2 = 0-strtol(s+1, &s, 10); break;
            case 'e': e1 = 0-strtol(optarg, &s, 10); if (*s == ',') e2 = 0-strtol(s+1, &s, 10); break;
            case 'n': n_seqs = atoi(optarg); break;
            case 's': seq_file = optarg; break;
            case 't': numThreads = atoi(optarg); break;
            case 'h': help(); return 0;
            default: help(); return 1;
        }
    }

    std::int8_t oe1=o1+e1, oe2=o2+e2;

    std::unique_ptr<spoa::AlignmentEngine> alignment_engine[numThreads];
    for (int i = 0; i < numThreads; i++) {
        try {
            alignment_engine[i] = spoa::createAlignmentEngine(
                    static_cast<spoa::AlignmentType>(algorithm), m, x, oe1, e1, oe2, e2);
        } catch(std::invalid_argument& exception) {
            std::cerr << exception.what() << std::endl;
            return 1;
        }
    }

    ifstream fp_seq;
    fp_seq.open(seq_file, ios::in);
    assert(fp_seq.is_open());

    vector<Batch> batches;

    struct timeval start_time, end_time, t_start, t_end;
    double runtime = 0; int seq_i; string seq;
    double realtime = 0, real_start, real_end;
    double graphCreationTime = 0, alignTime = 0, addToGraphTime = 0, generateConsensusTime = 0;

#pragma omp parallel num_threads(numThreads) 
{
    int tid = omp_get_thread_num();
    if (tid == 0) {
        fprintf(stderr, "Running with threads: %d\n", numThreads);
    }
}

    readFile(fp_seq, batches);

    gettimeofday(&start_time, NULL); real_start = get_realtime();

#ifdef VTUNE_ANALYSIS
    __itt_resume();
#endif
#pragma omp parallel num_threads(numThreads)
{
    int tid = omp_get_thread_num();
    #pragma omp for schedule(dynamic, 1)
        for (int i = 0; i < batches.size(); i++) {
            // gettimeofday(&t_start, NULL);
            auto graph = spoa::createGraph();
            // gettimeofday(&t_end, NULL);
            // graphCreationTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
            for (int j = 0; j < batches[i].seqs.size(); j++) {
                // gettimeofday(&t_start, NULL);
                auto alignment = alignment_engine[tid]->align(batches[i].seqs[j], graph);
                // gettimeofday(&t_end, NULL);
                // alignTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;

                // gettimeofday(&t_start, NULL);
                graph->add_alignment(alignment, batches[i].seqs[j]);
                // gettimeofday(&t_end, NULL);
                // addToGraphTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
            }
            // gettimeofday(&t_start, NULL);
            batches[i].consensus_seq = graph->generate_consensus();
            // gettimeofday(&t_end, NULL);
            // generateConsensusTime += (t_end.tv_sec - t_start.tv_sec)*1e6 + t_end.tv_usec - t_start.tv_usec;
        }   
}
#ifdef VTUNE_ANALYSIS
    __itt_pause();
#endif
    gettimeofday(&end_time, NULL); real_end = get_realtime();
    runtime += (end_time.tv_sec - start_time.tv_sec)*1e6 + end_time.tv_usec - start_time.tv_usec;
    realtime += (real_end-real_start);

#ifdef PRINT_OUTPUT
    for (int i = 0; i < batches.size(); i++) {
        cout << ">Consensus_sequence" << endl;
        cout << batches[i].consensus_seq.c_str() << endl;
    }
#endif

    fprintf(stderr, "Runtime: %.2f, GraphCreate: %.2f, Align: %.2f, AddSeqGraph: %.2f, Consensus %.2f %.2f %.3f \n", runtime*1e-6, graphCreationTime*1e-6, alignTime*1e-6, addToGraphTime*1e-6, generateConsensusTime*1e-6, realtime*1e-6, peakrss()/1024.0/1024.0);

    fp_seq.close();
    return 0;
}
