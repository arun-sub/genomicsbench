#ifndef DIAGONALS_H
#define DIAGONALS_H

#include <vector>
#include <iostream>

template<class PRECISION, class ALLOCATOR = std::allocator<PRECISION>>
struct Diagonals3 {
  std::vector<PRECISION, ALLOCATOR> m;
  std::vector<PRECISION, ALLOCATOR> mp;
  std::vector<PRECISION, ALLOCATOR> mpp;
  std::vector<PRECISION, ALLOCATOR> x;
  std::vector<PRECISION, ALLOCATOR> xp;
  std::vector<PRECISION, ALLOCATOR> xpp;
  std::vector<PRECISION, ALLOCATOR> y;
  std::vector<PRECISION, ALLOCATOR> yp;
  std::vector<PRECISION, ALLOCATOR> ypp;

  Diagonals3() = default;
  Diagonals3(const size_t read_length) :
    m   (read_length, static_cast<PRECISION>(0)),
    mp  (read_length, static_cast<PRECISION>(0)),
    mpp (read_length, static_cast<PRECISION>(0)),
    x   (read_length, static_cast<PRECISION>(0)),
    xp  (read_length, static_cast<PRECISION>(0)),
    xpp (read_length, static_cast<PRECISION>(0)),
    y   (read_length, static_cast<PRECISION>(0)),
    yp  (read_length, static_cast<PRECISION>(0)),
    ypp (read_length, static_cast<PRECISION>(0))
  {}

  void reserve(const size_t new_size) {
    // printf("%d\n", new_size);
    m.resize(new_size); mp.resize(new_size); mpp.resize(new_size);
    x.resize(new_size); xp.resize(new_size); xpp.resize(new_size);
    y.resize(new_size); yp.resize(new_size); ypp.resize(new_size);
    // m.reserve(new_size); mp.reserve(new_size); mpp.reserve(new_size);
    // x.reserve(new_size); xp.reserve(new_size); xpp.reserve(new_size);
    // y.reserve(new_size); yp.reserve(new_size); ypp.reserve(new_size);
    // printf("%ld %ld %ld %ld %ld %ld %ld %ld %ld\n", mpp.size(), xpp.size(), ypp.size(), mp.size(), xp.size(), yp.size(), m.size(), x.size(), y.size());
  }

  void rotate() {
    std::swap(mpp, mp);
    std::swap(xpp, xp);
    std::swap(ypp, yp);
    std::swap(mp, m);
    std::swap(xp, x);
    std::swap(yp, y);
  }

  void update(const PRECISION first_row_value) {
    std::fill(m.begin(), m.end(), static_cast<PRECISION>(0));
    std::fill(mp.begin(), mp.end(), static_cast<PRECISION>(0));
    std::fill(mpp.begin(), mpp.end(), static_cast<PRECISION>(0));
    std::fill(x.begin(), x.end(), static_cast<PRECISION>(0));
    std::fill(xp.begin(), xp.end(), static_cast<PRECISION>(0)); 
    std::fill(xpp.begin(), xpp.end(), static_cast<PRECISION>(0));
    std::fill(y.begin(), y.end(), static_cast<PRECISION>(0)); 
    std::fill(yp.begin(), yp.end(), static_cast<PRECISION>(0)); 
    std::fill(ypp.begin(), ypp.end(), static_cast<PRECISION>(0));
    
    // for (int i = 0; i < new_size; i++){
    //   m[i] = static_cast<PRECISION>(0);
    //   mp[i] = static_cast<PRECISION>(0);
    //   mpp[i] = static_cast<PRECISION>(0);
    //   x[i] = static_cast<PRECISION>(0);
    //   xp[i] = static_cast<PRECISION>(0);
    //   xpp[i] = static_cast<PRECISION>(0);
    //   y[i] = static_cast<PRECISION>(0);
    //   yp[i] = static_cast<PRECISION>(0);
    //   ypp[i] = static_cast<PRECISION>(0);
    // }

    ypp[0] = yp[0] = y[0] = first_row_value;

  }

};

template<class PRECISION, class ALLOCATOR = std::allocator<PRECISION>>
struct Diagonals2 {
  std::vector<PRECISION, ALLOCATOR> m;
  std::vector<PRECISION, ALLOCATOR> mp;
  std::vector<PRECISION, ALLOCATOR> x;
  std::vector<PRECISION, ALLOCATOR> xp;
  std::vector<PRECISION, ALLOCATOR> y;
  std::vector<PRECISION, ALLOCATOR> yp;

  Diagonals2() = default;
  Diagonals2(const size_t read_length) :
    m   (read_length, static_cast<PRECISION>(0)),
    mp  (read_length, static_cast<PRECISION>(0)),
    x   (read_length, static_cast<PRECISION>(0)),
    xp  (read_length, static_cast<PRECISION>(0)),
    y   (read_length, static_cast<PRECISION>(0)),
    yp  (read_length, static_cast<PRECISION>(0))
  {}

  void reserve(const size_t new_size) {
    m.reserve(new_size); mp.reserve(new_size);
    x.reserve(new_size); xp.reserve(new_size);
    y.reserve(new_size); yp.reserve(new_size);
  }

  void rotate() {
    std::swap(mp, m);
    std::swap(xp, x);
    std::swap(yp, y);
  }

  void update(const PRECISION first_row_value) {
    std::fill(m.begin(), m.end(), static_cast<PRECISION>(0));
    std::fill(mp.begin(), mp.end(), static_cast<PRECISION>(0));
    std::fill(x.begin(), x.end(), static_cast<PRECISION>(0));
    std::fill(xp.begin(), xp.end(), static_cast<PRECISION>(0));
    std::fill(y.begin(), y.end(), static_cast<PRECISION>(0));
    std::fill(yp.begin(), yp.end(), static_cast<PRECISION>(0));
    yp[0] = y[0] = first_row_value;
  }

};


template<class PRECISION>
std::ostream& operator<<(std::ostream& out, const Diagonals3<PRECISION>& diagonals);

#endif
