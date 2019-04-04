#include "solver.h"

int main(int argc, char**argv)
{
#ifndef noGPU
  device_matrix<double> A; device_vector<double> x, b;
#else
  matrix<double> A; vector<double> x, b;
#endif      
  //getProb(A,x,b);
  //PSC98init(A,x,b);

  static double *val; static int *I, *J;
  int M, N, nz, ret;

  ret = mm_read_unsymmetric_sparse(argv[1],&M,&N,&nz,&val,&I,&J);

  if ( ret != 0 ) printf("mm_read_unsymmetric_sparse()=%d %s\n",ret,
                         argv[1]);
  if ( M != N ) return 0;

  A.resize(N); x.resize(N); b.resize(N);

  for (int k=0;k<nz;k++)
    if ( 0<=I[k]&&I[k]<=N&&0<=J[k]&&J[k]<=N) A[I[k]][J[k]] = val[k];

  for (int k=0;k<N;k++) x[k] = 1.0;
  CRSinit(A);
  b = A*x;
  CRSdestory(A);
  for (int k=0;k<N;k++) x[k] = 0.0;
  M_init(A);
  solver(A,x,b);
  printv(x);
  //PSC98check(x);
  return 0;
}
