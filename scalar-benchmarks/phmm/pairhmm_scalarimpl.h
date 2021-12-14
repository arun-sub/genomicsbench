#ifndef PAIRHMM_SCALARIMPL_H
#define PAIRHMM_SCALARIMPL_H

#include <xmmintrin.h>
#include "pairhmm_impl.h"

#define TRISTATE_CORRECTION 3
#define doNotUseTristateCorrection 0

template <class PRECISION>
class PairhmmScalarImpl: public PairhmmImpl<PRECISION, Diagonals3<PRECISION>, Constants<PRECISION>, 1>  {
  using Base =  PairhmmImpl<PRECISION,Diagonals3<PRECISION>, Constants<PRECISION>, 1>;
public:
  PairhmmScalarImpl(const size_t initial_size = Base::INITIAL_SIZE): Base {initial_size} { }
  virtual ~PairhmmScalarImpl() { }
protected:

  // // Original
  // double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {

  //   const auto hl = haplotype.original_length;  // haplotype original length (unpadded)
  //   const auto rl = read.original_length;       // read original length (unpadded)
  //   const auto rows = rl + read.left_padding;  // number of rows in the diagonals (padded read length)
  //   const auto mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
  //   const auto fd = mrl - rl;                   // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
  //   auto result = 0.l;                          // result accumulator
  //   auto &diags = this->m_diagonals;
  //   auto &consts = this->m_constants;

  //   // for (auto d = fd; d != mrl + hl - 1; ++d) { // d for diagonal
  //   for (int d = 0; d != mrl + rl - 1; ++d) {
  //     const auto hap_offset = mrl+hl-1;
  //     for (auto r = 1u; r != rows; ++r) {       // r for row
  //       const auto read_base = read.bases[r];
  //       const auto hap_base = haplotype.bases[hap_offset+r-d];
  //       const auto base_qual = read.base_quals[r];
  //       const auto prior = ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ?  static_cast<PRECISION>(1) - base_qual : base_qual;
  //       diags.m[r] = prior * ((diags.mpp[r-1] * consts.mm[r]) + (consts.gm[r] * (diags.xpp[r-1] + diags.ypp[r-1])));
  //       diags.x[r] = diags.mp[r-1] * consts.mx[r] + diags.xp[r-1] * consts.xx[r];
  //       diags.y[r] = diags.mp[r] * consts.my[r] + diags.yp[r] * consts.yy[r];
  //     }
  //     result += diags.m[rows-1] + diags.x[rows-1];
  //     diags.rotate();
  //   }
  //   result = result < this->MIN_ACCEPTED ?
  //     this->FAILED_RUN_RESULT : // if we underflowed return failed constant to rerun with higher precision if desired
  //     log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
  //   printf("%lf\n", result);
  //   return result;
  // }




  // Reorder
  double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) override {

    int r, h;
    const int hl = haplotype.original_length;  // haplotype original length (unpadded)
    const int rl = read.original_length;       // read original length (unpadded)
    const int mrl = this->max_original_read_length();  // alias for max original read length for readability in the code below (max read length in the testcase)
    float result = 0.l;                          // result accumulator
    auto &consts = this->m_constants;

    std::vector<std::vector<PRECISION> > m;
    std::vector<std::vector<PRECISION> > x;
    std::vector<std::vector<PRECISION> > y;

    m.resize(rl+1, std::vector<PRECISION>(hl+1, 0));
    x.resize(rl+1, std::vector<PRECISION>(hl+1, 0));
    y.resize(rl+1, std::vector<PRECISION>(hl+1, 0));

    for (h = 0; h < hl+1; h++) {
      y[0][h] = constants_with_precision::INITIAL_CONSTANT_WITH_PRECISION<PRECISION>() / hl;
    }

    for (r = 1; r < rl+1; r++) {
      for (h = 1; h < hl+1; h++) {
        const int read_base = read.bases[r];
        const int hap_base = haplotype.bases[mrl+hl+1-h];
        const float base_qual = read.base_quals[r];
        const float prior = ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ?  static_cast<PRECISION>(1) - base_qual : base_qual / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION);
        m[r][h] = prior * (m[r-1][h-1] * consts.mm[r] + (x[r-1][h-1] + y[r-1][h-1]) * consts.gm[r]);
        
        x[r][h] = m[r-1][h] * consts.mx[r] + x[r-1][h] * consts.xx[r];
        y[r][h] = m[r][h-1] * consts.my[r] + y[r][h-1] * consts.yy[r];
      }
    }

    for (h = 1; h < hl+1; h++) {
      result += m[rl][h] + x[rl][h];
    }
    
    m.clear();
    x.clear();
    y.clear();


    // const int rows = rl + read.left_padding;
    // auto &diags = this->m_diagonals;
    // for (int d = 0; d != mrl + rl - 1; ++d) {
    //   const int hap_offset = mrl+hl-1;
    //   // printf("%d %f %f %f\n", d, diags.mpp[0], diags.xpp[0], diags.ypp[0]);
    //   // if (d < rows) printf("%f %f %f\n", diags.mpp[d], diags.xpp[d], diags.ypp[d]);
    //   for (r = 1u; r != rows; ++r) {       // r for row
    //     const int read_base = read.bases[r];
    //     const int hap_base = haplotype.bases[hap_offset+r-d];
    //     const PRECISION base_qual = read.base_quals[r];
    //     const PRECISION prior = ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ?  static_cast<PRECISION>(1) - base_qual : base_qual / (doNotUseTristateCorrection ? 1.0 : TRISTATE_CORRECTION);
    //     diags.m[r] = prior * ((diags.mpp[r-1] * consts.mm[r]) + (consts.gm[r] * (diags.xpp[r-1] + diags.ypp[r-1])));
    //     diags.x[r] = diags.mp[r-1] * consts.mx[r] + diags.xp[r-1] * consts.xx[r];
    //     diags.y[r] = diags.mp[r] * consts.my[r] + diags.yp[r] * consts.yy[r];
    //   }
    //   result += diags.m[rows-1] + diags.x[rows-1];
    //   diags.rotate();
    // }



    // const int rows = rl + read.left_padding;    // number of rows in the diagonals (padded read length)
    // // const int fd = mrl - rl;                    // first diagonal to compute (saves compute of all-0 diagonals when read is shorter than the padding - which will be the maximum read length in the testcase)
    // auto &diags = this->m_diagonals;
    // // for (int d = fd; d != mrl + hl - 1; ++d) { // d for diagonal
    // for (int d = 0; d != mrl + rl - 1; ++d) { // d for diagonal
    //   const int hap_offset = mrl+hl-1;
    //   // printf("%lf %lf %lf\n", diags.mpp[0], diags.xpp[0], diags.ypp[0]);
    //   // if (d < rows) printf("%lf %lf %lf\n", diags.mpp[d], diags.xpp[d], diags.ypp[d]);
    //   h = d-rl+2;
    //   for (int r = 1u; r != rows; ++r) {       // r for row
    //   // printf("%d %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", r, diags.mpp[r], diags.xpp[r], diags.ypp[r], diags.mp[r], diags.xp[r], diags.yp[r], diags.m[r], diags.x[r], diags.y[r]);
    //     const int read_base = read.bases[r];
    //     const int hap_base = haplotype.bases[hap_offset+r-d];
    //     const double base_qual = read.base_quals[r];
    //     const double prior = ((read_base == hap_base) || (read_base == 'N') || (hap_base == 'N')) ? static_cast<PRECISION>(1) - base_qual : base_qual;
    //     diags.m[r] = prior * ((diags.mpp[r-1] * consts.mm[r]) + (consts.gm[r] * (diags.xpp[r-1] + diags.ypp[r-1])));
    //     diags.x[r] = diags.mp[r-1] * consts.mx[r] + diags.xp[r-1] * consts.xx[r];
    //     diags.y[r] = diags.mp[r] * consts.my[r] + diags.yp[r] * consts.yy[r];
    //     if ((h > 0) && (h < hl+1) && (r > 0) && (r < rl+1)) {
    //       // printf("h=%d r=%d\n", h, r);
    //       m[r][d-r+2] = diags.m[r];
    //       x[r][d-r+2] = diags.x[r];
    //       y[r][d-r+2] = diags.y[r];          
    //       // printf("h=%d r=%d m[%d][%d]\n", hap_offset+r-d, r, r, d-r+2);
    //     }
    //     // if ((mrl < hap_offset+r-d) && (hap_offset+r-d <= mrl+hl)) printf("h=%d r=%d %lf %lf\n", hap_offset+r-d, r, diags.m[r], diags.x[r]);
    //   }
    //   // if (d+1 < rows) printf("h=1 %lf\n", diags.m[d+1]);
    //   result += diags.m[rl] + diags.x[rl];
    //   // if (d >= rl - 1) printf("h=%d r=%d %lf %lf %lf\n", mrl+hl+1-h, rl, result, diags.m[rl], diags.x[rl]);
    //   // if (d >= rl - 1) printf("h=%d r=%d %lf %lf %lf\n", mrl+hl+1-h, rl, result, m[rl][h], x[rl][h]);
    //   // if (d >= rl - 1) printf("h=%d r=%d %lf %lf %lf\n", h, rl, result, m[rl][h], x[rl][h]);
    //   diags.rotate();
    // }
    // // printf("\n");

    // // for (h = 1; h < hl+1; h++) {
    // //   for (r = 1; r < rl+1; r++) {
    // //     // printf("h=%d r=%d %lf %lf\n", mrl+hl+1-h, r, m[r][h], x[r][h]);
    // //     printf("h=%d r=%d %lf %lf\n", h, r, m[r][h], x[r][h]);
    // //   }
    // // }



    // printf("%lf %lf ", log10(static_cast<double>(result)), log10(static_cast<double>(this->INITIAL_CONSTANT)));
    result = result < this->MIN_ACCEPTED ?
      this->FAILED_RUN_RESULT : // if we underflowed return failed constant to rerun with higher precision if desired
      log10(static_cast<double>(result)) - log10(static_cast<double>(this->INITIAL_CONSTANT));
    // printf("%lf\n", result);
    return result;
  }

};

#endif
