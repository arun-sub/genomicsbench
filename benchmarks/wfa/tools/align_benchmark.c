/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Wavefront Alignment benchmarking tool
 */

#if VTUNE_ANALYSIS
  #include <ittnotify.h>
#endif

#include "utils/commons.h"
#include "gap_affine/affine_wavefront.h"
#include "gap_affine/affine_wavefront_align.h"

#include "omp.h"

/*
 * Parameters
 */
#define MAX_SEQUENCE_LENGTH 100000

/*
 * Generic parameters
 */
typedef struct {
  // I/O
  char *input;
  char *output;
  // Penalties
  affine_penalties_t affine_penalties;
  int min_wavefront_length;
  int max_distance_threshold;
  // Misc
  int nthreads;
  int progress;
  bool verbose;
} benchmark_args;
benchmark_args parameters = {
  // Input
  .input=NULL,
  .output=NULL,
  // Penalties
  .affine_penalties = {
      .match = 0,
      .mismatch = 4,
      .gap_opening = 6,
      .gap_extension = 2,
  },
  .min_wavefront_length = -1, // 10,
  .max_distance_threshold = -1, //50,
  // Misc
  .nthreads = 1,
  .progress = 10000,
  .verbose = false
};

/*
 * Parsing Input
 */

typedef struct {
  int id;
  char* pattern;
  int pattern_length;
  char* text;
  int text_length;
} input_pair_sequences_t;
input_pair_sequences_t** parse_input_sequences(
    FILE* const input_file,
    int** const total_sequences) {

  *total_sequences = malloc(parameters.nthreads * sizeof((*total_sequences)[0]));
  input_pair_sequences_t **input_buffers = malloc(parameters.nthreads * sizeof(input_pair_sequences_t*));

  char more_seqs = 1;
  int id = 0;
  #pragma omp parallel num_threads(parameters.nthreads)
  {
    // Init allocation
    int thread_id = omp_get_thread_num();
    int total_allocated = 10000;
    input_pair_sequences_t* input_buffer = malloc(total_allocated * sizeof(input_pair_sequences_t));

    // Input reading loop
    int total_parsed = 0;
    while(more_seqs) {
      // Avoid that some threads exit the while loop if another thread enters the
      // critical section before the threads check the while condition.
      #pragma omp barrier
      #pragma omp for schedule(static, 1)
      for (size_t i = 0; i < omp_get_num_threads(); ++i) {
        #pragma omp critical
        {
          char *line1 = NULL; 
          char *line2 = NULL;

          size_t line1_allocated = 0;
          size_t line2_allocated = 0;

          // Read queries
          int line1_length = getline(&line1, &line1_allocated, input_file);
          int line2_length = getline(&line2, &line2_allocated, input_file);
    
          if (line1_length == -1 || line2_length == -1) {
            free(line1);
            free(line2);

            // OMP: Implicit flush of more_seqs at the end of the for loop.
            more_seqs = 0;
          }
          else {
            // Process input
            char* const pattern = line1 + 1;
            const int pattern_length = line1_length - 2;
            pattern[pattern_length] = '\0';

            char* const text = line2 + 1;
            const int text_length = line2_length - 2;
            text[text_length] = '\0';

            // Allocate
            if (total_parsed + 1 >= total_allocated) {
              total_allocated += 10000;
              input_buffer = realloc(input_buffer, total_allocated * sizeof(input_pair_sequences_t));
            }

            // Copy Pattern
            input_buffer[total_parsed].id = id;
            ++id; // OMP: Implicit flush of id on entry and exit of the critical section.
            input_buffer[total_parsed].pattern = malloc(pattern_length + 1);
            strncpy(input_buffer[total_parsed].pattern, pattern, pattern_length + 1);
            input_buffer[total_parsed].pattern_length = pattern_length;

            // Copy Text
            input_buffer[total_parsed].text = malloc(text_length + 1);
            strncpy(input_buffer[total_parsed].text, text, text_length + 1);
            input_buffer[total_parsed].text_length = text_length;

            free(line1);
            free(line2);

            ++total_parsed;
          }
        }
      }
    }
    input_buffers[thread_id] = input_buffer;
    (*total_sequences)[thread_id] = total_parsed;
  }

  return input_buffers;
}
/*
 * Generic Menu
 */
void usage() {
  fprintf(stderr,
      "USE: ./align_benchmark -i <input> [-o <output>]                      \n"
      "      Options::                                                      \n"
      "        [I/O]                                                       \n"
      "          --input|i <File>                                           \n"
      "          --output|o <File>                                          \n"
      "        [Penalties]                                                  \n"
      "          --affine-penalties|g M,X,O,E                               \n"
      "        [Specifics]                                                  \n"
      "          --minimum-wavefront-length <INT>                           \n"
      "          --maximum-difference-distance <INT>                        \n"
      "        [Misc]                                                       \n"
      "          --nthreads|t <integer> [default=1]                         \n"
      "          --progress|P <integer>                                     \n"
      "          --verbose|v                                                \n"
      "          --help|h                                                   \n");
}
void parse_arguments(int argc,char** argv) {
  struct option long_options[] = {
    /* I/O */
    { "input", required_argument, 0, 'i' },
    { "output", required_argument, 0, 'o' },
    /* Penalties */
    { "affine-penalties", required_argument, 0, 'p' },
    /* Specifics */
    { "minimum-wavefront-length", required_argument, 0, 1000 },
    { "maximum-difference-distance", required_argument, 0, 1001 },
    /* Misc */
    { "nthreads", required_argument, 0, 't' },
    { "progress", required_argument, 0, 'P' },
    { "verbose", no_argument, 0, 'v' },
    { "help", no_argument, 0, 'h' },
    { 0, 0, 0, 0 } };
  int c,option_index;
  if (argc <= 1) {
    usage();
    exit(0);
  }
  while (1) {
    c=getopt_long(argc,argv,"i:o:p:t:P:vh",long_options,&option_index);
    if (c==-1) break;
    switch (c) {
    /*
     * I/O
     */
    case 'i':
      parameters.input = optarg;
      break;
    case 'o':
      parameters.output = optarg;
      break;
    /*
     * Penalties
     */
    case 'p': { // --affine-penalties
      char* sentinel = strtok(optarg,",");
      parameters.affine_penalties.match = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.mismatch = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_opening = atoi(sentinel);
      sentinel = strtok(NULL,",");
      parameters.affine_penalties.gap_extension = atoi(sentinel);
      break;
    }
    /*
     * Specific parameters
     */
    case 1000: // --minimum-wavefront-length
      parameters.min_wavefront_length = atoi(optarg);
      break;
    case 1001: // --maximum-difference-distance
      parameters.max_distance_threshold = atoi(optarg);
      break;
    /*
     * Misc
     */
    case 't':
      parameters.nthreads = atoi(optarg);
      break;
    case 'P':
      parameters.progress = atoi(optarg);
      break;
    case 'v':
      parameters.verbose = true;
      break;
    case 'h':
      usage();
      exit(1);
    // Other
    case '?': default:
      fprintf(stderr,"Option not recognized \n");
      exit(1);
    }
  }
  // Checks
  if (parameters.input == NULL) {
    fprintf(stderr,"Option --input is required \n");
    exit(1);
  }
}
int main(int argc,char* argv[]) {
#if VTUNE_ANALYSIS
  __itt_pause();
#endif
  // Parsing command-line options
  parse_arguments(argc,argv);
  // Parameters
  FILE *input_file = NULL;
  FILE *output_file = NULL;
  // Init I/O files
  input_file = fopen(parameters.input, "r");
  if (input_file == NULL) {
    fprintf(stderr,"Input file '%s' couldn't be opened\n",parameters.input);
    exit(1);
  }
  if (parameters.output != NULL) {
    output_file = fopen(parameters.output, "w");
    if (output_file == NULL) {
      fprintf(stderr,"Output file '%s' couldn't be opened\n",parameters.output);
      exit(1);
    }
  }

  struct timeval benchmark_start;
  gettimeofday(&benchmark_start, NULL);

  // Parse input file
  int *total_sequences;
  input_pair_sequences_t** const input_buffers =
      parse_input_sequences(input_file, &total_sequences);

  struct timeval alignment_start;
  struct timeval alignment_end;

  int progress_mod = 0;
  #pragma omp parallel num_threads(parameters.nthreads)
  {
    // Init MM-allocator
    mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);

    // Init wavefront
    affine_wavefronts_t* affine_wavefronts;
    if (parameters.min_wavefront_length < 0) {
      affine_wavefronts = affine_wavefronts_new_complete(
          MAX_SEQUENCE_LENGTH,MAX_SEQUENCE_LENGTH,
          &parameters.affine_penalties,mm_allocator);
    } else {
      affine_wavefronts = affine_wavefronts_new_reduced(
          MAX_SEQUENCE_LENGTH,MAX_SEQUENCE_LENGTH,
          &parameters.affine_penalties,parameters.min_wavefront_length,
          parameters.max_distance_threshold,mm_allocator);
    }

    #pragma omp barrier
    #pragma omp master
    {
#if VTUNE_ANALYSIS
      __itt_resume();
#endif
      gettimeofday(&alignment_start, NULL);
    }

    // Pointer to thread private data.
    int thread_id = omp_get_thread_num();
    input_pair_sequences_t *input_buffer = input_buffers[thread_id];
    int thread_total_sequences = total_sequences[thread_id];

    edit_cigar_t *edit_cigars = malloc(thread_total_sequences * sizeof(edit_cigars[0])); 

    // Read-align loop
    for (int i = 0; i < thread_total_sequences; ++i) {
      // Align
      affine_wavefronts_clear(affine_wavefronts);
      affine_wavefronts_align(affine_wavefronts,
          input_buffer[i].pattern, input_buffer[i].pattern_length,
          input_buffer[i].text, input_buffer[i].text_length);

      // Store output
      // Alloc new operations array.
      int noperations = 
        affine_wavefronts->edit_cigar.end_offset - 
        affine_wavefronts->edit_cigar.begin_offset;
      edit_cigars[i].operations = mm_allocator_malloc(mm_allocator, noperations);
      // Copy operations.
      memcpy(edit_cigars[i].operations, 
          &affine_wavefronts->edit_cigar.operations[affine_wavefronts->edit_cigar.begin_offset], 
          noperations * sizeof(edit_cigars->operations[0]));
      
      // Copy the rest of the fields.
      edit_cigars[i].max_operations = affine_wavefronts->edit_cigar.max_operations;
      edit_cigars[i].score = affine_wavefronts->edit_cigar.score;
      edit_cigars[i].begin_offset = 0;
      edit_cigars[i].end_offset = noperations;

      // Update progress
      if (parameters.verbose) {
        #pragma omp critical
        {
          ++progress_mod;
          if (progress_mod % parameters.progress == 0) {
            fprintf(stderr,"...processed %d reads \n", progress_mod);
          }
        }
      }
    }

    #pragma omp barrier
    #pragma omp master
    {
      gettimeofday(&alignment_end, NULL);
#if VTUNE_ANALYSIS
      __itt_pause();
#endif
    }

    // Print the output.
    #pragma omp critical
    {
      if (!parameters.verbose) {
        progress_mod += thread_total_sequences;
      }
      if (output_file != NULL) {
        for (int i = 0; i < thread_total_sequences; ++i) {
          fprintf(output_file, "id=%d ", input_buffer[i].id);
          edit_cigar_print(output_file, &edit_cigars[i]);
          fprintf(output_file, "\n");
        }
      }
    }

    // Free
    for (int i = 0; i < thread_total_sequences; ++i) {
      mm_allocator_free(mm_allocator, edit_cigars[i].operations);
      free(input_buffer[i].pattern);
      free(input_buffer[i].text);
    }
    free(edit_cigars);

    free(input_buffer);

    affine_wavefronts_delete(affine_wavefronts);
    mm_allocator_delete(mm_allocator);
  }
  free(total_sequences);
  free(input_buffers);
  fclose(input_file);
  if (output_file != NULL) fclose(output_file);

  struct timeval benchmark_end;
  gettimeofday(&benchmark_end, NULL);

  printf("Total.reads: %d\n", progress_mod);
  printf("Time.Benchmark: %f s\n", (benchmark_end.tv_sec - benchmark_start.tv_sec) + 
      (benchmark_end.tv_usec - benchmark_start.tv_usec) * 1E-6);
  printf("Time.Alignment: %f s\n", (alignment_end.tv_sec - alignment_start.tv_sec) + 
      (alignment_end.tv_usec - alignment_start.tv_usec) * 1E-6);
}
