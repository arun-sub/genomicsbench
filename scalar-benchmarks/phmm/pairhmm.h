#ifndef PAIRHMM_H
#define PAIRHMM_H

#include <vector>

#include "pairhmm_impl.h"
#include "testcase.h"

template <typename FAST, typename CORRECT>
class Pairhmm {
 private:
  FAST fast;
  CORRECT correct;
 public:
  std::vector<double> calculate(const Testcase& testcase) {
    auto results = fast.calculate(testcase); // calculate all the testcases using the float precision pairhmm implementation
    correct.recalculate(testcase, results);  // recalculate only the ones that failed with double precision (results = DOUBLE_MIN for failed ones)
    return results;
  }
};

#endif
