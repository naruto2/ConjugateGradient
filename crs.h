#ifndef CSRMATRIX_H
#define CSRMATRIX_H
using namespace std;

typedef struct {
  cublasHandle_t     cublas;
  cusparseHandle_t   cusparse;
  cusparseMatDescr_t matDescr;
  long N;
  long n;
  double* val;
  int*    col_ind;
  int*    row_ptr;
} CRSdata;


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


void CRSdataDestory(CRSdata& CRS) {
  cublasDestroy(CRS.cublas);
  cusparseDestroy(CRS.cusparse);
  cusparseDestroyMatDescr(CRS.matDescr);
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
