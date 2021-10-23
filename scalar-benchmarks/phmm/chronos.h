#ifndef CHRONOS_H
#define CHRONOS_H

#include <chrono>     // std::chrono for timing benchmarks

class Chronos {
  public:
    Chronos(void): m_start(now()) { }

    void reset(void) {
        m_start = now();
    }

    double elapsed(void) const {
      using namespace std::chrono;
      double t = duration_cast<microseconds>(now()-m_start).count();
      return t/1000.;
    }
  private:
    using clock = std::chrono::high_resolution_clock;
    using time_point = std::chrono::time_point<clock>;
    static time_point now() { return clock::now(); }
    time_point m_start;
};

#endif
