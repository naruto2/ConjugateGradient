#ifndef CSRMATRIX_H
#define CSRMATRIX_H
using namespace std;

typedef struct {
  cublasHandle_t     cublas;
  cusparseHandle_t   cusparse;
  cusparseMatDescr_t matDescr;
  long N;
  long nonZeroCount;
  double* elementsPtr;
  int* columnIndecesPtr;
  int* rowOffsetsPtr;
} CsrMatrix;


int cublas_cusparse_init(CsrMatrix* Ap)
{
  /************************************/
  /********** cuBLASの準備 ************/
  /************************************/
  // cuBLASハンドルを作成
  cublasStatus_t stat =cublasCreate(&Ap->cublas);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    cout<<"cuBLAS initialization failed"<<endl;
    return EXIT_FAILURE;
  }
  /************************************/
  /********** cuSPARSEの準備 **********/
  /************************************/
  // cuSPARSEハンドルを作成
  cusparseStatus_t cusparsestat = cusparseCreate(&Ap->cusparse);
  if (cusparsestat != CUSPARSE_STATUS_SUCCESS) {
    cout<<"cuSPARSE initialization failed"<<endl;
    return EXIT_FAILURE;
  }
  // 行列形式を作成
  // * 一般的な形式
  // * 番号は0から開始
  cusparsestat =   cusparseCreateMatDescr(&Ap->matDescr);
  if (cusparsestat != CUSPARSE_STATUS_SUCCESS) {
    cout<<"cuSPARSE initialization failed(cusparseCreateMatDescr)"<<endl;
    return EXIT_FAILURE;
  }
  cusparseSetMatType(Ap->matDescr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(Ap->matDescr, CUSPARSE_INDEX_BASE_ZERO);

  return 0;
}

void y_Ax(double* resultPtr, CsrMatrix A, double* vectorPtr)
{
  // Csrmv（CSR形式行列とベクトルの積）を実行
  // y = α*Ax + β*y;
  const double ALPHA = 1;
  const double BETA = 0;
  ::cusparseDcsrmv(A.cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE,
		   A.N, A.N, A.nonZeroCount,
		   &ALPHA, A.matDescr, A.elementsPtr, A.rowOffsetsPtr, A.columnIndecesPtr,
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
