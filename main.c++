#include "solver.h"

int main()
{
#ifndef noGPU
  device_matrix<double> A; device_vector<double> x, b;
#else
  matrix<double> A; vector<double> x, b;
#endif      
  getProb(A,x,b);
  ConjugateGradient(A,x,b);
  printv(x);
  return 0;
}
