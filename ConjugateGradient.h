#define dev(someDev) thrust::raw_pointer_cast(&(someDev)[0])
using namespace std;
using namespace thrust;

int ConjugateGradient(sparse::matrix<double> spA, double* x, double* b)
{
  long N = spA.size();
  double alpha, beta, normr0, normr, rhop, rho, temp, tol = 0.000000001;
  long maxit = 2*N;
  CsrMatrix A;
  if ( cublas_cusparse_init(&A) != 0 ) return EXIT_FAILURE;
  
  /**********************************/
  /********** 入力値の転送 **********/
  /**********************************/
  // GPU側の配列を確保
  // （ポインタ管理が面倒なのでthrust使うと便利！）
  long n=0;
  for (long i=0; i<N; i++) for ( auto it : spA[i] ) n++;

  device_vector<double> elementsDev(n);
  device_vector<int>    columnIndecesDev(n);
  device_vector<int>    rowOffsetsDev(N+1);
  device_vector<double> xnDev(N);  
  device_vector<double> xDev(N);
  device_vector<double> bDev(N);
  device_vector<double> zDev(N);
  device_vector<double> pDev(N);
  device_vector<double> qDev(N);
  device_vector<double> rDev(N);

  double*  elements      = new double[n];
  int*     columnIndeces = new int[n];
  int*     rowOffsets    = new int[N+1];
  
  n=0;
  for (long i=0; i<N; i++){
    rowOffsets[i] = n;
    for ( auto it : spA[i] ) {
      long j=it.first;
      elements[n] = spA[i][j];
      columnIndeces[n] = j;
      n++;
    }
  }
  rowOffsets[N] = n;
  
  // 解ベクトルを設定
  double* xn = new double[N];
  for (int i = 0; i < N; i++) xn[i] = i * 0.1;
  
  // GPU側配列へ入力値（行列とベクトル）を複製
  copy_n(elements,      rowOffsets[N], elementsDev.begin());
  copy_n(columnIndeces, rowOffsets[N], columnIndecesDev.begin());
  copy_n(rowOffsets,    N+1, rowOffsetsDev.begin());
  copy_n(xn, N, xnDev.begin());
  copy_n(x,  N,  xDev.begin());
  copy_n(b,  N,  bDev.begin());

  /******************************************/
  /********** 行列×ベクトルの計算 **********/
  /******************************************/
  // thrust配列からCUDA用ポインタに変換
  A.N                = N;
  A.nonZeroCount     = rowOffsets[N];
  A.elementsPtr      = dev(elementsDev);
  A.columnIndecesPtr = dev(columnIndecesDev);
  A.rowOffsetsPtr    = dev(rowOffsetsDev);

  // Csrmv（CSR形式行列とベクトルの積）を実行
  // y = Ax
  { long i;
    for (i = 0; i< N; i++) if ( b[i] != 0.0 ) break;
    if ( i >= N ) y_Ax(dev(bDev), A, dev(xnDev));
  }

  y_Ax(dev(rDev), A, dev(xDev));

  alpha = -1.0;
  cublasDscal(A.cublas, N, &alpha, dev(rDev), 1);

  alpha = 1.0;
  cublasDaxpy(A.cublas,N, &alpha, dev(bDev), 1, dev(rDev), 1);

  cublasDnrm2(A.cublas, N, dev(rDev), 1, &normr0);
  if ( normr0 == 0.0 ) return 0;

  for (long i=0; i<maxit; i++){
    cublasDcopy(A.cublas, N, dev(rDev), 1, dev(zDev), 1);
    // 本来は Mz = r

    //4: \rho = r^{T} z
    rhop = rho;
    cublasDdot(A.cublas, N, dev(rDev), 1, dev(zDev), 1, &rho);

    if ( i == 0 ){
      //6: p = z
      cublasDcopy(A.cublas, N, dev(zDev), 1, dev(pDev), 1);
    } else {
      //8: \beta = rho_{i} / \rho_{i-1}
      beta = rho/rhop;
      //9: p = z + \beta p
      cublasDaxpy(A.cublas, N, &beta, dev(pDev), 1, dev(zDev), 1);
      cublasDcopy(A.cublas, N, dev(zDev), 1, dev(pDev), 1);
    }

    //11: Compute q = Ap (sparse matrix-vector multiplication)
    y_Ax(dev(qDev), A, dev(pDev));

    //12: \alpha = \rho_{i} / (p^{T} q)
    cublasDdot(A.cublas, N, dev(pDev), 1, dev(qDev), 1, &temp);

    //13: x = x + \alpha p
    alpha = rho/temp;
    cublasDaxpy(A.cublas, N, &alpha, dev(pDev), 1, dev(xDev), 1);

    //14: r = r - \alpha q
    alpha = -alpha;
    cublasDaxpy(A.cublas, N, &alpha, dev(qDev), 1, dev(rDev), 1);
    
    //check for convergence
    cublasDnrm2(A.cublas, N, dev(rDev), 1, &normr);
    //cout << normr/normr0 <<endl;
    if ( normr/normr0 < tol ) break;
  }
  cublasDestroy(A.cublas);
  cusparseDestroy(A.cusparse);
  cusparseDestroyMatDescr(A.matDescr);
  copy_n(xDev.begin(), N, x);
  delete[] elements;
  delete[] columnIndeces;
  delete[] rowOffsets;
  delete[] xn;
  return 0;
}
