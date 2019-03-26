
#ifndef noGPU
double dot(const device_vector<double>&x, const device_vector<double>&y){
  double rho;
  long N = x.size();
  cublasDdot(gCRS.cublas, N, dev(x), 1, dev(y), 1, &rho);
  return rho;
}

double nrm2(const device_vector<double>& r) {
  double norm;
  long N = r.size();
  cublasDnrm2(gCRS.cublas, N, dev(r), 1, &norm);
  return norm;
}


device_vector<double>& operator*(const device_matrix<double>& A,
				 device_vector<double>&b) {
  long N = b.size();
  static device_vector<double> x(N);
  y_Ax(dev(x), gCRS, dev(b));
  return x;
}



device_vector<double>& operator*(double a, device_vector<double>&x){
  int n = x.size();
  cublasDscal(gCRS.cublas, n, &a, dev(x), 1);
  return x;
}



device_vector<double>& operator+(device_vector<double>&y,
				 device_vector<double>&x){
  int n = x.size();
  double alpha = 1.0;
  cublasDaxpy(gCRS.cublas, n, &alpha, dev(x), 1, dev(y), 1);
  return y;
}


device_vector<double>& operator-(device_vector<double>&y,
				 device_vector<double>&x){
  int n = x.size();
  double alpha = -1.0;
  cublasDaxpy(gCRS.cublas, n, &alpha, dev(x), 1, dev(y), 1);
  return y;
}

#endif

vector<double>& operator+(vector<double>& y, vector<double>& x){
  int i, n = x.size();
  static vector<double> z(n);
  for ( i = 0; i < n; i++ ) z[i] = y[i] + x[i];
  return z;
}



vector<double>& operator-(vector<double>& y,  vector<double>& x) {
  int i, n = x.size();
  static vector<double> z(n);
  for( i = 0; i < n; i++ ) z[i] = y[i] - x[i];
  return z;
}


double dot(const vector<double>&x, const vector<double>&y){
  double rho = 0.0;
  long i, n = x.size();
  for ( i = 0; i < n; i++ ) rho += x[i]*y[i];
  return rho;
}

double nrm2(const vector<double>& r) {
  double norm = 0.0;
  long i, n = r.size();
  for ( i = 0; i < n; i++ ) norm += r[i]*r[i];
  return sqrt(norm);
}


vector<double>& operator*(const matrix<double>& A, vector<double>& b) {
  long n = b.size();
  static vector<double> x(n);
  for (long i=0; i<n; i++ ) {
    auto Ai = A[i];
    auto j = Ai.begin();
    for ( x[i]=0.0; j != Ai.end(); j++ )
      x[i] += j->second * b[j->first];
  }
  return x;
}
vector<double>& operator*(double a, vector<double>&x){
  long i, n = x.size();
  static vector<double> y(n);
  for ( i = 0; i < n; i++ ) y[i] = a*x[i];
  return y;
}

