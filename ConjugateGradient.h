#define dev(someDev) thrust::raw_pointer_cast(&(someDev)[0])

using namespace std;
using namespace thrust;


long countelements(matrix<double> &A) {
  long N=A.size();
  long n=0;
  for (long i=0; i<N; i++) for ( auto it : A[i] ) n++;
  return n;
}


void matrix2CRS(matrix<double> &A,
		device_vector<double>&  valDev,
		device_vector<int>& col_indDev,
		device_vector<int>& row_ptrDev,
		CRSdata &CRS) {
  long N = A.size();
  long n = countelements(A);

  double*  val = new double[n];;
  int*     col_ind = new int[n];
  int*     row_ptr = new int[N+1];

  if ( val == NULL || col_ind == NULL || row_ptr == NULL ) {
    fprintf(stderr,"Out of memory\n");
    exit(1);
  }
  
  n=0;
  for (long i=0; i<N; i++){
    row_ptr[i] = n;
    for ( auto it : A[i] ) {
      long j=it.first;
      val[n] = A[i][j];
      col_ind[n] = j;
      n++;
    }
  }
  row_ptr[N] = n;

  copy_n(val,       n, valDev.begin());
  copy_n(col_ind,   n, col_indDev.begin());
  copy_n(row_ptr, N+1, row_ptrDev.begin());
  delete[] val;
  delete[] col_ind;
  delete[] row_ptr;

  CRS.N                = N;
  CRS.n                = n;
  CRS.val              = dev(valDev);
  CRS.col_ind          = dev(col_indDev);
  CRS.row_ptr          = dev(row_ptrDev);
}

#include "d.h"
#include "operator.h"

int ConjugateGradient(matrix<double>& A, double* x, double* b)
{
  long N = A.size();
  long n = countelements(A);
  
  double alpha, beta, normr0, normr, rhop, rho, tol = 0.000000001;
  long maxit = 2*N;
  CRSdata CRS;
  if ( cublas_cusparse_init(CRS) != 0 )
    return EXIT_FAILURE;

  device_vector<double> val(n);
  device_vector<int>    col_ind(n);
  device_vector<int>    row_ptr(N+1);
  matrix2CRS(A, val, col_ind, row_ptr, CRS);

  /**********************************/
  /********** 入力値の転送 **********/
  /**********************************/
  // GPU側の配列を確保
  // （ポインタ管理が面倒なのでthrust使うと便利！）

  device_vector<double> xnDev(N);  
  device_vector<double> xDev(N);
  device_vector<double> bDev(N);
  device_vector<double> zDev(N);
  device_vector<double> pDev(N);
  device_vector<double> qDev(N);
  device_vector<double> rDev(N);


  // 解ベクトルを設定
  double* xn = new double[N];
  for (int i = 0; i < N; i++) xn[i] = i * 0.1;
  
  // GPU側配列へ入力値（行列とベクトル）を複製
  copy_n(xn, N, xnDev.begin());
  delete[] xn;

  copy_n(x,  N,  xDev.begin());
  copy_n(b,  N,  bDev.begin());

  /******************************************/
  /********** 行列×ベクトルの計算 **********/
  /******************************************/
  // thrust配列からCUDA用ポインタに変換

  // CRSmv（CRS形式行列とベクトルの積）を実行
  // y = Ax
  Dclass D(CRS);
  gCRS = CRS;

  long i;
  for (i = 0; i< N; i++) if ( b[i] != 0.0 ) break;
  if ( i >= N ) bDev =  A * xnDev;

  rDev = A * xDev;

  rDev = -1.0 * rDev;

  rDev = rDev + bDev;
  
  normr0 = nrm2(rDev);

  if ( normr0 == 0.0 ) return 0;

  for (long i=0; i<maxit; i++){
    // 本来は Mz = r
    zDev = rDev;

    rhop = rho;

    rho = dot(rDev, zDev);

    if ( i == 0 ){
      pDev = zDev;
    } else {
      beta = rho/rhop;

      zDev = zDev + beta * pDev;

      pDev = zDev;
    }

    qDev = A * pDev;

    alpha = rho/dot(pDev,qDev);

    xDev = xDev + alpha * pDev;

    rDev = rDev - alpha * qDev;

    normr = nrm2(rDev);

    if ( normr/normr0 < tol ) break;
  }
  CRSdataDestory(CRS);
  copy_n(xDev.begin(), N, x);
  return 0;
}
