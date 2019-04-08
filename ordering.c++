#include "ordering.h"

long count_kern(exp_matrix<double>& A)
{
  long k;
  for(k=1; A[k][k]; k++) ;
  return k;
}


void gen_minimatrix(exp_matrix<double>& A0, exp_matrix<double>& A)
{
  long k =  count_kern(A0);
  printf("k= %d\n",k);

  long n = A0.size()-k+1;
  A.resize(n);
  A.ROW.resize(n), A.COL.resize(n), A.nROW.resize(n), A.nCOL.resize(n);

  long i, j;

  for(i=1; i<n; i++) for(j=1; j<n; j++) A[i][j] = A0[i+k-1][j+k-1];
  
  for(i=1; i<n; i++) {
    A.nROW[i] = A0.nROW[i+k-1];
    A.nCOL[i] = A0.nCOL[i+k-1];
    A.ROW[i] = A0.ROW[i+k-1];
    A.COL[i] = A0.COL[i+k-1];
  }
    
  printm(A);
}



int main()
{
  exp_matrix<double> A0, A;
  
  setval(A0);
  alg_single(A0);
  printm(A0);

  gen_minimatrix(A0,A);
  
  
  return 0;
}

