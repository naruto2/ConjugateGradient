#ifndef SIMPLE_H
#define SIMPLE_H
#include <map>
#include <vector>

namespace simple {

  using namespace std;

  template<typename Real> class matrix : public vector< map<long, Real> > {
  public:
    matrix(long n) : vector< map<long, Real> >(n){}
  };
#include "printmatrix.h"
}
#endif
