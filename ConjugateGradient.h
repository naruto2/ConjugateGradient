template<typename Matrix, typename Vector> 
int ConjugateGradient(Matrix& A, Vector& x, Vector& b)
{
  CRSinit(A);

  long i, n=A.size(), maxit=2*n, ret=0;
  double alpha, beta, normr0, rhop, rho, tol=0.00000000001;
  Vector p(n), q(n), r(n), z(n);  

  r = b - A*x;

  normr0 = nrm2(r);
  
  if ( normr0 == 0.0 ) goto end;
  
  for (i = 0; i < maxit; i++){
    
    z = M_solve(r);
    
    rhop = rho;
    
    rho = dot(r, z);
    
    
    if ( i == 0 ){
      p = z;
    } else {
      beta = rho/rhop;

      axpy(n, &beta, p, 1, z, 1);

      p = z;
    }

    q = A * p;

    alpha = rho/dot(p,q);

    axpy(n, &alpha, p, 1, x, 1);

    alpha = -alpha;

    axpy(n, &alpha, q, 1, r, 1);
    
    if ( nrm2(r)/normr0 < tol ) goto end;
  }
  ret=1;

 end:
  CRSdestory(A);
  return ret;
}
