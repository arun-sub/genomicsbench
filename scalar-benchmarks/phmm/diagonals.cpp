#include <ostream>

#include "diagonals.h"

using namespace std;

template<class PRECISION>
ostream& operator<<(ostream& out, const Diagonals3<PRECISION>& diagonals) {
  out << "M   -> "; for (auto x : diagonals.m) out << x << ","; out << endl;
  out << "Mp  -> "; for (auto x : diagonals.mp) out << x << ","; out << endl;
  out << "Mpp -> "; for (auto x : diagonals.mpp) out << x << ","; out << endl;
  out << "X   -> "; for (auto x : diagonals.x) out << x << ","; out << endl;
  out << "Xp  -> "; for (auto x : diagonals.xp) out << x << ","; out << endl;
  out << "Xpp -> "; for (auto x : diagonals.xpp) out << x << ","; out << endl;
  out << "Y   -> "; for (auto x : diagonals.y) out << x << ","; out << endl;
  out << "Yp  -> "; for (auto x : diagonals.yp) out << x << ","; out << endl;
  out << "Ypp -> "; for (auto x : diagonals.ypp) out << x << ","; out << endl;
  return out;
}
