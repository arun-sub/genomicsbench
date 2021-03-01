#ifndef SHACC_PAIRHMM_H
#define SHACC_PAIRHMM_H

#ifdef __APPLE__
#define WEAK __attribute__((weak_import))
#else
#define WEAK __attribute__((weak))
#endif

namespace shacc_pairhmm {

  struct Read {
    int length;
    const char* bases;
    const char* q;
    const char* i;
    const char* d;
    const char* c;
  };
  
  struct Haplotype {
    int length;
    const char* bases;
  };
  
  struct Batch {
    int id;
    int num_reads;
    int num_haps;
    long num_cells;
    Read* reads;
    Haplotype* haps;
    double* results;

    bool operator < (const Batch& b) const
    {
        return (num_cells < b.num_cells);
    }
  };

    struct SortByCells
    {
        bool operator()( const Batch& lx, const Batch& rx ) const {
            return lx.num_cells < rx.num_cells;
        }
    };

    struct SortById
    {
        bool operator()( const Batch& lx, const Batch& rx ) const {
            return lx.id < rx.id;
        }
    };

  extern WEAK bool calculate(Batch& batch);
}

#endif
