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
 * DESCRIPTION: Wavefront alignment algorithm for pairwise gap-affine
 *   alignment (Main module)
 */

#include "affine_wavefront.h"
#include "affine_wavefront_backtrace.h"
#include "affine_wavefront_utils.h"

/*
 * Setup
 */
void affine_wavefronts_allocate_wavefront_null(
    affine_wavefronts_t* const affine_wavefronts) {
  // Allocate null wavefront
  const int wavefront_length = affine_wavefronts->pattern_length + affine_wavefronts->text_length + 1;
  awf_offset_t* const offsets_null = mm_allocator_calloc(
      affine_wavefronts->mm_allocator,wavefront_length,awf_offset_t,false);
  // Initialize
  affine_wavefronts->wavefront_null.null = true;
  affine_wavefronts->wavefront_null.lo =  1;
  affine_wavefronts->wavefront_null.hi = -1;
  affine_wavefronts->wavefront_null.lo_base =  1;
  affine_wavefronts->wavefront_null.hi_base = -1;
  affine_wavefronts->wavefront_null.offsets = offsets_null + affine_wavefronts->pattern_length; // Center at k=0
  int i;
  for (i=0;i<wavefront_length;++i) {
    offsets_null[i] = AFFINE_WAVEFRONT_OFFSET_NULL;
  }
}
void affine_wavefronts_allocate_wavefront_components(
    affine_wavefronts_t* const affine_wavefronts) {
  // Parameters
  mm_allocator_t* const mm_allocator = affine_wavefronts->mm_allocator;
  // Initialize wavefronts
  const int num_wavefronts = affine_wavefronts->num_wavefronts;
  affine_wavefronts->mwavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine_wavefront_t*,true);
  affine_wavefronts->iwavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine_wavefront_t*,true);
  affine_wavefronts->dwavefronts =
      mm_allocator_calloc(mm_allocator,num_wavefronts,affine_wavefront_t*,true);
  // Allocate bulk-memory (for all wavefronts)
  affine_wavefront_t* const wavefronts_mem =
      mm_allocator_calloc(mm_allocator,3*num_wavefronts,affine_wavefront_t,false);
  affine_wavefronts->wavefronts_mem = wavefronts_mem;
  affine_wavefronts->wavefronts_current = wavefronts_mem;
}
affine_wavefronts_t* affine_wavefronts_new(
    const int pattern_length,
    const int text_length,
    affine_penalties_t* const penalties,
    const wavefronts_penalties_strategy penalties_strategy,
    mm_allocator_t* const mm_allocator) {
  // Handler
  affine_wavefronts_t* const affine_wavefronts = mm_allocator_alloc(mm_allocator,affine_wavefronts_t);
  // MM
  affine_wavefronts->mm_allocator = mm_allocator;
  affine_wavefronts->mm_stack = mm_stack_new(BUFFER_SIZE_8M);
  // Dimensions
  const int max_score_misms = MIN(pattern_length,text_length) * penalties->mismatch;
  const int max_score_indel = penalties->gap_opening + ABS(pattern_length-text_length) * penalties->gap_extension;
  const int num_wavefronts = max_score_misms + max_score_indel;
  affine_wavefronts->pattern_length = pattern_length;
  affine_wavefronts->text_length = text_length;
  affine_wavefronts->num_wavefronts = num_wavefronts;
  affine_wavefronts->max_allocated_wavefront = 0;
  // Limits
  const int single_gap_penalty = penalties->gap_opening + penalties->gap_extension;
  const int max_penalty = MAX(penalties->mismatch,single_gap_penalty);
  affine_wavefronts->max_penalty = max_penalty;
  // Penalties
  affine_wavefronts_penalties_init(&affine_wavefronts->penalties,penalties,penalties_strategy);
  // Allocate wavefronts
  affine_wavefronts_allocate_wavefront_components(affine_wavefronts);
  affine_wavefronts_allocate_wavefront_null(affine_wavefronts);
  // CIGAR
  edit_cigar_allocate(&affine_wavefronts->edit_cigar,pattern_length,text_length,mm_allocator);
  // Return
  return affine_wavefronts;
}
void affine_wavefronts_clear(
    affine_wavefronts_t* const affine_wavefronts) {
  // Clear wavefronts memory
  const int num_wavefronts = MIN(affine_wavefronts->max_allocated_wavefront,affine_wavefronts->num_wavefronts);
  memset(affine_wavefronts->mwavefronts,0,num_wavefronts*sizeof(affine_wavefront_t*));
  memset(affine_wavefronts->iwavefronts,0,num_wavefronts*sizeof(affine_wavefront_t*));
  memset(affine_wavefronts->dwavefronts,0,num_wavefronts*sizeof(affine_wavefront_t*));
  mm_stack_clear(affine_wavefronts->mm_stack);
  // Clear CIGAR
  edit_cigar_clear(&affine_wavefronts->edit_cigar);
  // Clear wavefronts-ptr
  affine_wavefronts->wavefronts_current = affine_wavefronts->wavefronts_mem;
}
void affine_wavefronts_delete(
    affine_wavefronts_t* const affine_wavefronts) {
  // Parameters
  mm_allocator_t* const mm_allocator = affine_wavefronts->mm_allocator;
  // Free MID-Wavefronts
  mm_allocator_free(mm_allocator,affine_wavefronts->mwavefronts);
  mm_allocator_free(mm_allocator,affine_wavefronts->iwavefronts);
  mm_allocator_free(mm_allocator,affine_wavefronts->dwavefronts);
  mm_allocator_free(mm_allocator,affine_wavefronts->wavefront_null.offsets - affine_wavefronts->pattern_length);
  // Free bulk memory
  mm_allocator_free(mm_allocator,affine_wavefronts->wavefronts_mem);
  // CIGAR
  edit_cigar_free(&affine_wavefronts->edit_cigar,mm_allocator);
  // MM
  mm_stack_delete(affine_wavefronts->mm_stack);
  // Handler
  mm_allocator_free(mm_allocator,affine_wavefronts);
}
/*
 * Setup WF-modes
 */
affine_wavefronts_t* affine_wavefronts_new_complete(
    const int pattern_length,
    const int text_length,
    affine_penalties_t* const penalties,
    mm_allocator_t* const mm_allocator) {
  // Create new
  affine_wavefronts_t* const affine_wavefronts =
      affine_wavefronts_new(
          pattern_length,text_length,
          penalties,wavefronts_penalties_force_zero_match,mm_allocator);
  // Limits
  affine_wavefronts->max_k = text_length;
  affine_wavefronts->min_k = -pattern_length;
  // Reduction
  affine_wavefronts_reduction_set_none(&affine_wavefronts->reduction);
  // Return
  return affine_wavefronts;
}
affine_wavefronts_t* affine_wavefronts_new_reduced(
    const int pattern_length,
    const int text_length,
    affine_penalties_t* const penalties,
    const int min_wavefront_length,
    const int max_distance_threshold,
    mm_allocator_t* const mm_allocator) {
  // Create new
  affine_wavefronts_t* const affine_wavefronts =
      affine_wavefronts_new(
          pattern_length,text_length,
          penalties,wavefronts_penalties_force_zero_match,mm_allocator);
  // Limits
  affine_wavefronts->max_k = text_length;
  affine_wavefronts->min_k = -pattern_length;
  // Reduction
  affine_wavefronts_reduction_set_dynamic(
      &affine_wavefronts->reduction,min_wavefront_length,max_distance_threshold);
  // Return
  return affine_wavefronts;
}
/*
 * Allocate individual wavefront
 */
affine_wavefront_t* affine_wavefronts_allocate_wavefront(
    affine_wavefronts_t* const affine_wavefronts,
    const int lo_base,
    const int hi_base) {
  // Compute limits
  const int wavefront_length = hi_base - lo_base + 2; // (+1) for k=0
  // Allocate wavefront
  affine_wavefront_t* const wavefront = affine_wavefronts->wavefronts_current;
  ++(affine_wavefronts->wavefronts_current); // Next
  // Configure offsets
  wavefront->null = false;
  wavefront->lo = lo_base;
  wavefront->hi = hi_base;
  wavefront->lo_base = lo_base;
  wavefront->hi_base = hi_base;
  // Allocate offsets
  awf_offset_t* const offsets_mem = mm_stack_calloc(
      affine_wavefronts->mm_stack,wavefront_length,awf_offset_t,false);
  awf_offset_t* const offsets = offsets_mem - lo_base; // Center at k=0
  wavefront->offsets = offsets;
  // Return
  return wavefront;
}

