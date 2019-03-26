#ifndef CSR_H
#define CSR_H
#define dev(someDev) thrust::raw_pointer_cast(&(someDev)[0])

class CRSdata {
 public:
  cublasHandle_t     cublas;
  cusparseHandle_t   cusparse;
  cusparseMatDescr_t matDescr;
  long N;
  long n;
  double* val;
  int*    col_ind;
  int*    row_ptr;
};


static CRSdata gCRS;

int cublas_cusparse_init(CRSdata& CRS)
{
  /************************************/
  /********** cuBLASの準備 ************/
  /************************************/
  // cuBLASハンドルを作成
  cublasStatus_t stat =cublasCreate(&CRS.cublas);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    cout<<"cuBLAS initialization failed"<<endl;
    return EXIT_FAILURE;
  }
  /************************************/
  /********** cuSPARSEの準備 **********/
  /************************************/
  // cuSPARSEハンドルを作成
  cusparseStatus_t cusparsestat = cusparseCreate(&CRS.cusparse);
  if (cusparsestat != CUSPARSE_STATUS_SUCCESS) {
    cout<<"cuSPARSE initialization failed"<<endl;
    return EXIT_FAILURE;
  }
  // 行列形式を作成
  // * 一般的な形式
  // * 番号は0から開始
  cusparsestat =   cusparseCreateMatDescr(&CRS.matDescr);
  if (cusparsestat != CUSPARSE_STATUS_SUCCESS) {
    cout<<"cuSPARSE initialization failed(cusparseCreateMatDescr)"<<endl;
    return EXIT_FAILURE;
  }
  cusparseSetMatType(CRS.matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(CRS.matDescr, CUSPARSE_INDEX_BASE_ZERO);

  return 0;
}



long countelements(device_matrix<double> &A) {
  long N=A.size();
  long n=0;
  for (long i=0; i<N; i++) for ( auto it : A[i] ) n++;
  return n;
}


void matrix2CRS(device_matrix<double> &A, CRSdata &CRS) {
  long N = A.size();
  long n = countelements(A);

  static device_vector<double>     valDev(n);
  static device_vector<int>    col_indDev(n);
  static device_vector<int>    row_ptrDev(N+1);

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

void CRSinit(matrix<double>& A)
{
}

void CRSinit(device_matrix<double>& A)
{
  if ( cublas_cusparse_init(gCRS) != 0 )
    exit(EXIT_FAILURE);
  matrix2CRS(A, gCRS);
}


void CRSdataDestory(device_matrix<double>& A) {
  cublasDestroy(gCRS.cublas);
  cusparseDestroy(gCRS.cusparse);
  cusparseDestroyMatDescr(gCRS.matDescr);
}

void CRSdataDestory(matrix<double>& A) {
}

void y_Ax(double* resultPtr, CRSdata CRS, double* vectorPtr)
{
  // CRSmv（CRS形式行列とベクトルの積）を実行
  // y = α*Ax + β*y;
  const double ALPHA = 1;
  const double BETA = 0;
  cusparseDcsrmv(CRS.cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
		 CRS.N, CRS.N, CRS.n,
		 &ALPHA, CRS.matDescr, CRS.val, CRS.row_ptr, CRS.col_ind,
		 vectorPtr,
		 &BETA, resultPtr);
}


void printv(long N, double* result)
{
  // 結果の表示
  for(long i = 0; i < N; i++)
    {
      std::cout << result[i] << std::endl;
    }
}


void printDev(long N, thrust::device_vector<double> someDev)
{
  double* some = new double[N];
  thrust::copy_n(someDev.begin(), N, some);
  printv(N, some);
  delete[] some;
}
#endif
