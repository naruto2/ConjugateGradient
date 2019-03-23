#include <iostream>

template<typename Real> void printmatrix(matrix<Real> &A) {
  long i, j, n;
  n = A.size();
  for(i=0; i<n; i++){
    for(j=0; j<n; j++) cout << A[i][j] << " ";
    cout << endl;
  }
}
