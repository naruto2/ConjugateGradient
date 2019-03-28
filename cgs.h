template < class Vector >
void printv(Vector &x)
{
  long i, n = x.size();
  for ( i = 0; i < n; i++ ) {
    cout<< x[i] << endl;
  }
  exit(1);
}

#ifndef noGPU
void cp(device_vector<double> &x, device_vector<double> &y){
  int I = 1, n = x.size();
  cublasDcopy(gCRS.cublas, n, dev(x), I, dev(y), I);
}
#endif
	  
void cp(vector<double> &x, vector<double> &y){
  long i, n = x.size();
  for ( i = 0; i < n; i++ ) y[i] = x[i];
}


#define M_solve(p) p
template < class Matrix, class Vector>
int CGS(Matrix &A, Vector &x, Vector &b)
{
  CRSinit(A);

  long i, n = A.size();
  Vector x0(n);
  for (i = 0; i < n; i++) x0[i] = i*0.1;

  for (i = 0; i< n; i++)
    if ( b[i] != 0.0 ) break;
  if ( i >= n ) b =  A * x0;

  double resid, tol = 0.000000000001, rho_1, rho_2, alpha, beta, normb=nrm2(b);
  Vector p(n), phat(n), q(n), qhat(n), vhat(n), u(n), uhat(n), r(n), rtilde(n),
    tmp;
  long max_iter = 2*n;

  tmp = A * x;
  cp(tmp,r);
  
  r = -1.0*r;
  r = r + b;
  cp(r,rtilde);

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = nrm2(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }
  for (int i = 1; i <= max_iter; i++) {
    rho_1 = dot(rtilde, r);
    if (rho_1 == 0) {
      tol = nrm2(r) / normb;
      return 2;
    }
    
    if (i == 1) {
      cp(r,u);
      cp(r,p);
    } else {
      beta = rho_1 / rho_2;
      cp(r,u); cp(q,tmp); u = u + beta * tmp;
               cp(p,tmp); q = q + beta * tmp;
      cp(u,p); cp(q,tmp); p = p + beta * tmp;
    }
    cp(M_solve(p), phat);

    cp(A*phat, vhat);
    
    alpha = rho_1 / dot(rtilde, vhat);

    cp(u,uhat);
    cp(vhat, tmp);
    u = u - alpha * tmp;
    cp(u, q);

    uhat = M_solve(uhat + q);

    cp(uhat, tmp); x = x + alpha * tmp;
    
    cp(A*uhat, qhat);

    r = r - alpha * qhat;

    rho_2 = rho_1;
    if ((resid = nrm2(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  CRSdestory(A);
  return 1;
}
