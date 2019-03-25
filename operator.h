static CRSdata gCRS;


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

device_vector<double>& operator*(const matrix<double>& A,
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

