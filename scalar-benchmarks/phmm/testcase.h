#ifndef TESTCASE_H
#define TESTCASE_H

#include <ostream>
#include <vector>

#include "read.h"
#include "haplotype.h"

struct Testcase {
  std::vector<Read<uint8_t,uint8_t>> reads;
  std::vector<Haplotype<uint8_t>> haplotypes;

  Testcase() {}
  Testcase(std::vector<Read<uint8_t,uint8_t>>& reads_, std::vector<Haplotype<uint8_t>>& haplotypes_) = delete;
  Testcase(std::vector<Read<uint8_t,uint8_t>>&& reads_, std::vector<Haplotype<uint8_t>>&& haplotypes_) :
    reads {std::move(reads_)},
    haplotypes {std::move(haplotypes_)}
  {}

  size_t size() const {return haplotypes.size() * reads.size();}
};

std::ostream& operator<<(std::ostream& out, const Testcase& testcase); 

#endif
