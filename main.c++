#ifndef noGPU
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <thrust/device_vector.h>
using namespace thrust;
#endif

#include <iostream>
#include <cmath>
#include "sparse.h"
using namespace std;
using namespace sparse;
#include "crs.h"
#include "operator.h"
#include "ConjugateGradient.h"
#include "cgs.h"

int main()
{
#ifndef noGPU
  device_matrix<double> A; device_vector<double> x, b;
#else
  matrix<double> A; vector<double> x, b;
#endif      
  getProb(A,b);

  // 解ベクトル(初期値)の設定
  long i, n = A.size();
  x.resize(n);
  for (i = 0; i < n; i++) x[i] = 1.0;
  
  CGS(A,x,b);

  printv(x);
  return 0;
}
