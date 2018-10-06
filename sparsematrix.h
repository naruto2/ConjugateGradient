#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H
#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <unistd.h>

using namespace std;

namespace sparse {

  template<typename T>
  class matrix : public vector< map<long, T> > {

  public:
    matrix(){}
    matrix(long n) : vector< map<long, T> >(n){};
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
