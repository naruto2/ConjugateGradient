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
#include "psc98.h"


template < class Matrix, class Vector >
double
resi(Matrix& A, Vector& x, const Vector& b)
{
  Vector r = y_ax(-1.0*(A*x), 1.0, b);
  return nrm2(r)/nrm2(b);
}

template < class Matrix, class Vector >
int
solver(Matrix& A, Vector& x, const Vector&b)
{
  Vector y = x; int ret = 0;

  ret = QMR(A,x,b);
  if ( ret == 0 ) return 0;
  CRSinit(A);
  if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
  CRSdestory(A);
  
  ret = BiCG(A,x,b); 
  if ( ret == 0 ) return 0;
  CRSinit(A);
  if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
  CRSdestory(A);
  
  ret = BiCGSTAB(A,x,b);
  if ( ret == 0 ) return 0; 
  CRSinit(A);
  if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
  CRSdestory(A);
  
  ret = CGS(A,x,b); 
  if ( ret == 0 ) return 0;
  CRSinit(A);
  if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
  CRSdestory(A);
  
  ret = ConjugateGradient(A,x,b);
  if ( ret == 0 ) return 0;
  CRSinit(A);
  if (resi(A,x,b)>resi(A,y,b)) x=y; else y=x;
  CRSdestory(A);
  
  return 1;
}

int progress(string str, long i, double res)
{
  cout<<str<<" i= "<<i<<" res= "<<res<<endl;
  if ( str == "QMR" ) {
    if ( res < 0.1) return 1;
    if ( i > 10000) return 2;
  }
  if ( str == "BiCG" ) {
    if (res < 0.00001) return 1;
    if ( i > 20000) return 2;
  }
  if ( str == "BiCGSTAB" ) {
    if (res < 0.0000001) return 1;
    if ( i > 20000) return 2;
  }
  if ( str == "CGS" ) {
    if (res < 0.00000001) return 1;
    if ( i > 20000) return 2;
  }
  return 0;
}

#endif
