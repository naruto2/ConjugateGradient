#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <unistd.h>
#include <iterator>

//#define MYMAP map<long,T>
#define MYMAP mymap<T>

using namespace std;

namespace sparse {

  template<typename Real>
  class mymap : public map<long, Real> {

    map<long,Real> mp;
  public:
    mymap(){}

    Real& operator [] ( const long i )  // 書き込み用
    {
      return mp[i];
    }
    map<long,double>::iterator begin(){ return mp.begin();}
    map<long,double>::iterator end(){ return mp.end();}
  };

  
  template<typename T>
  class matrix : public vector< MYMAP > {

  public:
    matrix(){}
    matrix(long n) : vector< MYMAP >(n){};
  };

  template<typename T>
  void printmatrix(matrix<T>&A)
  {
    long i, j, n;
    n = A.size();
    for(i=0;i<n;i++){
      for(j=0;j<n;j++) cout<<A[i][j]<<" ";
      cout << endl;
    }
  };
}
  
#endif
