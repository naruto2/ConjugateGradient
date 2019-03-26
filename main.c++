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
using namespace std;
using namespace sparse;
#include "crs.h"
#include "operator.h"
#include "ConjugateGradient.h"

int main()
{
      /**********************************/
      /********** 入力値の準備 **********/
      /**********************************/
      long i, n = 1024;
#ifndef noGPU
      device_matrix<double> A(n);      
      device_vector<double> x(n), b(n);
#else
      matrix<double> A(n);       
      vector<double> x(n), b(n);
#endif      
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

      for (i = 0; i < n; i++) {
	  A[i][i] = 2;
	  if(i > 0) A[i][i-1] = 1;
	  if(i < n-1) A[i][i+1] = 1;
      }

      // 解ベクトル(初期値)の設定
      for (i = 0; i < n; i++) x[i] = 1.0;

      // 右辺ベクトルを設定 (b==0ならばbを自動生成)
      for (i = 0; i < n; i++) b[i] = 0.0;

      ConjugateGradient(A,x,b);

      for (i = 0; i < n; i++) cout << x[i] << endl;
      return 0;
}
