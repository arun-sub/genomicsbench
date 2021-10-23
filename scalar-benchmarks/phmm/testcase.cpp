#include <ostream>

#include "testcase.h"

std::ostream& operator<<(std::ostream& out, const Testcase& testcase) {
  out << "Reads:" << std::endl;
  for (const auto& read : testcase.reads)
    out << " " << read << std::endl;
  out << "Haplotypes:" << std::endl;
  for (const auto& haplotype : testcase.haplotypes)
    out << " " << haplotype << std::endl;
  return out;
}
