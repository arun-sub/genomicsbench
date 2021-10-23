#ifndef PAIRHMM_IMPL_H
#define PAIRHMM_IMPL_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ostream>
#include <cassert>

#include "testcase.h"
#include "constants.h"
#include "diagonals.h"

#define P(x) x

namespace constants_with_precision {
  template <class T>
  static constexpr T INITIAL_CONSTANT_WITH_PRECISION();

  template <>
  constexpr double INITIAL_CONSTANT_WITH_PRECISION<double>() {return std::numeric_limits<double>::max() / 16;}

  template <>
  constexpr float INITIAL_CONSTANT_WITH_PRECISION<float>() {return std::numeric_limits<float>::max() / 16;}

  template <class T>
  static constexpr T MIN_ACCEPTED_WITH_PRECISION();

  template <>
  constexpr double MIN_ACCEPTED_WITH_PRECISION<double>() {return 0.l;}

  template <>
  constexpr float MIN_ACCEPTED_WITH_PRECISION<float>() {return 1e-28f;}
}

template <class PRECISION, class DIAGONALS, class CONSTANTS, int VECSIZE = 1>
class PairhmmImpl {
  size_t m_max_original_read_length = 0;  // updated during calculate() -- inside pad_haplotypes

public:
  size_t max_original_read_length(void) const { return m_max_original_read_length; }

  PairhmmImpl(const size_t initial_size = INITIAL_SIZE) :
    m_constants {initial_size},
    m_diagonals {initial_size},
    m_ph2pr{} {
    m_ph2pr.reserve(MAX_PH2PR_INDEX);
    for (auto i=static_cast<PRECISION>(0); i!=MAX_PH2PR_INDEX; ++i)
      m_ph2pr.push_back(pow(static_cast<PRECISION>(10.0), (-i) / static_cast<PRECISION>(10.0)));
  }

  virtual ~PairhmmImpl() { }

  std::vector<double> calculate (const Testcase& testcase) {
    const auto padded_haplotypes = pad_haplotypes(testcase); // updates m_max_original_read_length (important!)
    const auto padded_reads = pad_reads(testcase.reads);
    const auto max_padded_read_length = m_max_original_read_length + total_read_padding();
    m_constants.reserve(max_padded_read_length+(VECSIZE-1));
    m_diagonals.reserve(max_padded_read_length+(VECSIZE-1));
    return calculate(padded_reads, padded_haplotypes);
  }

  void recalculate (const Testcase& testcase, std::vector<double>& results) {
    m_max_original_read_length = calculate_max_read_length(testcase.reads);
    const auto max_padded_read_length = m_max_original_read_length + total_read_padding();
    m_constants.reserve(max_padded_read_length+(VECSIZE-1));
    m_diagonals.reserve(max_padded_read_length+(VECSIZE-1));
    calculate_failed(testcase, results);
  }
 protected:
  CONSTANTS m_constants;
  DIAGONALS m_diagonals;
  std::vector<PRECISION> m_ph2pr;

  static constexpr auto READ_LEFT_PADDING  = 1;
  static constexpr auto READ_RIGHT_PADDING = VECSIZE-1;
  static constexpr auto MAX_PH2PR_INDEX    = 128;
  static constexpr auto INITIAL_CONSTANT   = constants_with_precision::INITIAL_CONSTANT_WITH_PRECISION<PRECISION>();
  static constexpr auto INITIAL_SIZE       = 1;
  static constexpr auto MIN_ACCEPTED       = constants_with_precision::MIN_ACCEPTED_WITH_PRECISION<PRECISION>();
  static constexpr auto FAILED_RUN_RESULT  = std::numeric_limits<double>::min();

  size_t total_read_padding() {
    return READ_LEFT_PADDING + READ_RIGHT_PADDING;
  }

  template<class T>
    void pad (T& v, size_t padding) const {
      v.insert(v.end(), padding, 0);
    }

  Haplotype<PRECISION> pad_and_reverse_haplotype(const Haplotype<uint8_t>& haplotype, const size_t left_padding, const size_t right_padding) const {
    const auto padded_length = left_padding + haplotype.bases.size() + right_padding;
    auto p = Haplotype<PRECISION>{};
    p.original_length = haplotype.bases.size();
    p.bases.reserve(padded_length);
    pad(p.bases, left_padding);
    p.bases.insert(p.bases.end(), haplotype.bases.rbegin(), haplotype.bases.rend());
    pad(p.bases, right_padding);
    return p;
  }

  const Haplotype<PRECISION> pad_haplotype(const Haplotype<uint8_t>& haplotype) const {
    const auto left_padding = m_max_original_read_length + 1;
    const auto right_padding = m_max_original_read_length + (VECSIZE-1);
    return pad_and_reverse_haplotype(haplotype, left_padding, right_padding);
  }

  template<class Result_Type, const bool convert = true>
    std::vector<Result_Type> pad_and_convert_qual(const std::vector<uint8_t>& qual, const size_t left_padding, const size_t right_padding) const {
      const auto padded_length = left_padding + qual.size() + right_padding;
      auto p = std::vector<Result_Type>{};
      p.reserve(padded_length);
      pad(p, left_padding);
      for (const auto q : qual)
        p.push_back(convert ? m_ph2pr[q] : q);
      pad(p, right_padding);
      return p;
    }

  Read<PRECISION,PRECISION> pad_read(const Read<uint8_t,uint8_t>& read) const {
    auto padded_read = Read<PRECISION,PRECISION>{};
    padded_read.left_padding = READ_LEFT_PADDING;
    padded_read.original_length = read.bases.size();
    padded_read.right_padding = READ_RIGHT_PADDING;
    padded_read.bases      = pad_and_convert_qual<PRECISION, false>(read.bases, READ_LEFT_PADDING, READ_RIGHT_PADDING);
    padded_read.base_quals = pad_and_convert_qual<PRECISION>(read.base_quals, READ_LEFT_PADDING, READ_RIGHT_PADDING);
    padded_read.ins_quals  = pad_and_convert_qual<PRECISION>(read.ins_quals, READ_LEFT_PADDING, READ_RIGHT_PADDING);
    padded_read.del_quals  = pad_and_convert_qual<PRECISION>(read.del_quals, READ_LEFT_PADDING, READ_RIGHT_PADDING);
    padded_read.gcp_quals  = pad_and_convert_qual<PRECISION>(read.gcp_quals, READ_LEFT_PADDING, READ_RIGHT_PADDING);
    return padded_read;
  }

  std::vector<Read<PRECISION,PRECISION>> pad_reads(const std::vector<Read<uint8_t, uint8_t>>& reads) const {
    auto padded_reads = std::vector<Read<PRECISION,PRECISION>>{};
    padded_reads.reserve(reads.size());
    for (auto& read : reads)
      padded_reads.push_back(pad_read(read));
    return padded_reads;
  }

  std::vector<Haplotype<PRECISION>> pad_haplotypes(const Testcase& testcase) {
    auto padded_haplotypes = std::vector<Haplotype<PRECISION>>{};
    m_max_original_read_length = calculate_max_read_length(testcase.reads);
    padded_haplotypes.reserve(testcase.haplotypes.size());
    for (const auto& haplotype : testcase.haplotypes)
      padded_haplotypes.push_back(pad_haplotype(haplotype));
    return padded_haplotypes;
  }

  template <class T, class U>
    size_t calculate_max_read_length(const std::vector<Read<T,U>>& reads) const {
      auto max_read_length = 0u;
      for (const auto& read : reads)
        max_read_length = max_read_length >= read.bases.size() ? max_read_length : read.bases.size();
      return max_read_length;
    }

  double compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) {
    return do_compute_full_prob(read, haplotype);
  }

  std::vector<double> calculate (const std::vector<Read<PRECISION,PRECISION>>& padded_reads, const std::vector<Haplotype<PRECISION>>& padded_haplotypes) {
    auto results = std::vector<double>{};
    results.reserve(padded_reads.size() * padded_haplotypes.size());
    for (const auto& read : padded_reads) {
      m_constants.update(read);
      for (const auto& hap : padded_haplotypes) {
        m_diagonals.update(INITIAL_CONSTANT / hap.original_length);
        results.push_back(compute_full_prob(read, hap));
      }
    }
    return results;
  }

  void calculate_failed (const Testcase& testcase, std::vector<double>& results) {
    auto master_idx = 0;
    auto padded_read = Read<PRECISION,PRECISION>{};
    for (auto read : testcase.reads) {
      auto has_padded_read = false;
      for (auto hap_idx = 0u; hap_idx != testcase.haplotypes.size(); ++hap_idx) {
        if (results[master_idx] == FAILED_RUN_RESULT) {
          if (!has_padded_read) {
            padded_read = pad_read(read); // making a copy here so I can keep it across haplotype iterations
            m_constants.update(padded_read);
            has_padded_read = true;
          }
          auto padded_haplotype = pad_haplotype(testcase.haplotypes[hap_idx]); // can try to optimize this later, but only after benchmarking
          m_diagonals.update(INITIAL_CONSTANT/padded_haplotype.original_length);
          results[master_idx] = compute_full_prob(padded_read, padded_haplotype);
        }
        ++master_idx;
      }
    }
  }


protected:

  virtual double do_compute_full_prob(const Read<PRECISION,PRECISION>& read, const Haplotype<PRECISION>& haplotype) = 0;

};

#endif
