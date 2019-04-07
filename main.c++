#include "solver.h"

int main(int argc, char**argv)
{
#ifndef noGPU
  device_matrix<double> A; device_vector<double> x, b;
#else
  matrix<double> A; vector<double> x, b;
#endif      
  MatrixMarket(argv[1],A,x,b);
  solver(A,x,b);
  printv(x);
  return 0;
}
