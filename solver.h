#ifndef SOLVER_H
#define SOLVER_H

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
#define SOLVERROR_NONE      0
#define SOLVERROR_MAXIT     1
#define SOLVERROR_BREAKDOWN 2
using namespace std;
using namespace sparse;
int progress(string str, long i, double res);
#include "crs.h"
#include "operator.h"
#include "ConjugateGradient.h"
#include "cgs.h"
#include "bicgstab.h"
#include "gmres.h"
#include "bicg.h"
#include "qmr.h"

template < class Matrix, class Vector >
int
solver(Matrix& A, Vector& x, const Vector&b)
{
  if( QMR(A,x,b) == 0) return 0;
  return GMRES(A,x,b);
}

int progress(string str, long i, double res)
{
  cout<<str<<" i= "<<i<<" res= "<<res<<endl;
  if ( i > 4000 ) return 1;
  return 0;
}

#endif
