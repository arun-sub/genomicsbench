#include <iostream>
#include <vector>
#include <cmath>

#ifndef DENORMALS
#include <xmmintrin.h>
#include <immintrin.h>
#endif

#include "input_reader.h"
#include "testcase_iterator.h"
#include "testcase.h"
#include "pairhmm.h"
#include "pairhmm_scalarimpl.h"
#include "chronos.h"

using namespace std;

int main (const int argc, char const * const argv[]) {
#ifndef DENORMALS
    _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
    _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
#endif
  auto pairhmm = Pairhmm<
    PairhmmScalarImpl<float>,
    PairhmmScalarImpl<double>
  >{};
  InputReader<TestcaseIterator> reader {};
  if (argc == 2)
    reader.from_file(argv[1]);
  double computation_time = 0.f;
  Chronos time;
  for (auto& testcase : reader) {
    time.reset();
    auto results = pairhmm.calculate(testcase);
    computation_time += time.elapsed();
    for (auto x : results)
      cout << x << endl;
  }
  std::cerr << "done in " << computation_time << "ms\n";
  return 0;
}
