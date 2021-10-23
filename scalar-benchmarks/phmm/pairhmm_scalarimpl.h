#ifndef PAIRHMM_SCALARIMPL_H
#define PAIRHMM_SCALARIMPL_H

#include <xmmintrin.h>
#include "pairhmm_impl.h"

template <class PRECISION>
class PairhmmScalarImpl: public PairhmmImpl<PRECISION, Diagonals3<PRECISION>, Constants<PRECISION>, 1>  {
  using Base =  PairhmmImpl<PRECISION,Diagonals3<PRECISION>, Constants<PRECISION>, 1>;
public:
  PairhmmScalarImpl(const size_t initial_size = Base::INITIAL_SIZE): Base {initial_size} { }
  virtual ~PairhmmScalarImpl() { }
protected:
  double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {
    const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
    const auto rl = read.original_length;       // read original length (unpadded)
    const auto rows = rl + read.left_padding;  // number of rows in the diagonals (padded read length)
    const auto mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
    const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    auto result = 0.l;                          // result accumulator
    auto &diags = this->m_diagonals;
    auto &consts = this->m_constants;
    for (auto d = fd; d != mrl + hl - 1; ++d) { // d for diagonal
      const auto hap_offset = mrl+hl-1;
      for (auto r = 1u; r != rows; ++r) {       // r for row
        const auto read_base = read.bases[r];
        const auto hap_base = haplotype.bases[hap_offset+r-d];
        const auto base_qual = read.base_quals[r];
        const auto prior = ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ?  static_cast<PRECISION>(1) - base_qual : base_qual;
        diags.m[r] = prior * ((diags.mpp[r-1] * consts.mm[r]) + (consts.gm[r] * (diags.xpp[r-1] + diags.ypp[r-1])));
        diags.x[r] = diags.mp[r-1] * consts.mx[r] + diags.xp[r-1] * consts.xx[r];
        diags.y[r] = diags.mp[r] * consts.my[r] + diags.yp[r] * consts.yy[r];
      }
      result += diags.m[rows-1] + diags.x[rows-1];
      diags.rotate();
    }
    return result < this->MIN_ACCEPTED ?
      this->FAILED_RUN_RESULT : // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
  }

};

#endif
