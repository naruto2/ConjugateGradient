template<typename T> 
int ConjugateGradient(matrix<double>& A, T& x, T& b)
{
  CRSinit(A);

  long i, n = A.size();
  T x0(n), z(n), p(n), q(n), r(n);
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
  
  for (i = 0; i<maxit; i++){

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
  return 0;
}
