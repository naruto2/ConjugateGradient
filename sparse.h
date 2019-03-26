#ifndef SPARSE_H
#define SPARSE_H
#include <map>
#include <vector>

namespace sparse {

  using namespace std;

  template<typename Real> class customap {
    map<long, Real> mp;
  public:
    auto& operator []( const long i ) {
      for ( auto it : mp )
	if (it.second == 0.0) mp.erase(it.first);
      return mp[i];
    }
    auto begin() { return mp.begin();}
    auto end()   { return mp.end();}
  };

  template<typename Real> class matrix : public vector< customap<Real> > {
  public:
    matrix(long n) : vector< customap<Real> >(n){}
  };

  template<typename Real> class device_matrix : public matrix<Real> {
  public:
    device_matrix(long n) : matrix<Real>(n){}
  };
  
}
#endif
