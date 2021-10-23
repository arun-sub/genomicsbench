#ifndef HAPLOTYPE_H
#define HAPLOTYPE_H

#include <ostream>
#include <vector>
#include <string>

#include "utils.h"

template <class QUAL_TYPE>
struct Haplotype {
  size_t original_length;
  std::vector<QUAL_TYPE> bases;

  Haplotype() = default;

  explicit Haplotype(const std::vector<QUAL_TYPE>& bases_, const size_t original_length_ = 0) : 
    original_length{original_length_},
    bases {bases_}
  {}

  explicit Haplotype(const std::string& bases_, const size_t original_length_ = 0) :
    original_length{original_length_},
    bases {convert_bytes<std::vector<QUAL_TYPE>>(bases_)}
  {}
};

template <class QUAL_TYPE>
std::ostream& operator<<(std::ostream& out, const Haplotype<QUAL_TYPE>& haplotype) { 
  out << convert_bytes<std::string>(haplotype.bases);
  return out;
}

#endif
