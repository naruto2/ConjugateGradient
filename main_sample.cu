#include <iostream>
#include <cuda_runtime_api.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <thrust/device_vector.h>
#include "csrmatrix.h"
#include "sparsematrix.h"
#include "ConjugateGradient.h"


int main()
{
      /**********************************/
      /********** 入力値の準備 **********/
      /**********************************/
      long N = 1024;
      sparse::matrix<double> A(N);

	// 中央差分行列を準備する
	//（対角項が2でその隣が1になる、↓こんなやつ）
	// | 2 1 0 0 0 0 0 0 ・・・ 0 0 0|
	// | 1 2 1 0 0 0 0 0 ・・・ 0 0 0|
	// | 0 1 2 1 0 0 0 0 ・・・ 0 0 0|
	// | 0 0 1 2 1 0 0 0 ・・・ 0 0 0|
	// | 0 0 0 1 2 1 0 0 ・・・ 0 0 0|
	// | 0 0 0 0 1 2 1 0 ・・・ 0 0 0|
	// | 0 0 0 0 0 1 2 1 ・・・ 0 0 0|
	// | 0 0 0 0 0 0 1 2 ・・・ 0 0 0|
	// | 0 0 0 0 0 0 0 0 ・・・ 2 1 0|
	// | 0 0 0 0 0 0 0 0 ・・・ 1 2 1|
	// | 0 0 0 0 0 0 0 0 ・・・ 0 1 2|

      for (long i = 0; i < N; i++) {
	  A[i][i] = 2;
	  if(i > 0) A[i][i-1] = 1;
	  if(i < N-1) A[i][i+1] = 1;
      }

      // 解ベクトル(初期値)の設定
      double* x = new double[N];
      for (long i = 0; i < N; i++) x[i] = 1.0;

      // 右辺ベクトルを設定 (b==0ならばbを自動生成)
      double* b = new double[N];
      for (long i = 0; i < N; i++) b[i] = 0.0;

      ConjugateGradient(A,x,b);
      for (long i = 0; i < N; i++) cout << x[i] << endl;
      delete[] x;
      delete[] b;
      return 0;
}
