int ConjugateGradient(matrix<double>& A, double* x1, double* b1)
{
  CRSinit(A);

  long i, n = A.size();
  device_vector<double> x0(n), x(n), b(n), z(n), p(n), q(n), r(n);
  copy_n(x1, n, x.begin()); copy_n(b1, n, b.begin());
  double alpha, beta, normr0, normr, rhop, rho, tol = 0.000000001;
  long maxit = 2*n;

  for (i = 0; i < n; i++) x0[i] = i*0.1;
  
  for (i = 0; i< n; i++)
    if ( b[i] != 0.0 ) break;
  if ( i >= n ) b =  A * x0;
  
  r = A * x;
  
  r = -1.0 * r;
  
  r = r + b;
  
  normr0 = nrm2(r);
  
  if ( normr0 == 0.0 ) goto end;
  
  for (long i=0; i<maxit; i++){

    z = r;    // 本来は Mz = r

    rhop = rho;

    rho = dot(r, z);

    if ( i == 0 ){
      p = z;
    } else {
      beta = rho/rhop;

      z = z + beta * p;

      p = z;
    }

    q = A * p;

    alpha = rho/dot(p,q);

    x = x + alpha * p;

    r = r - alpha * q;
    
    normr = nrm2(r);
    
    if ( normr/normr0 < tol ) break;
  }

 end:
  CRSdataDestory();
  copy_n(x.begin(), n, x1);
  return 0;
}
