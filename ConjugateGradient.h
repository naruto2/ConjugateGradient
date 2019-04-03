template<typename Matrix, typename Vector> 
int ConjugateGradient(Matrix& A, Vector& x, const Vector& b)
{
  CRSinit(A);

  long i, n=A.size(), maxit=2*n, ret=SOLVERROR_NONE;
  double alpha, beta, normr0, rhop, rho, tol=0.00000000001;
  Vector p, q, r, z;  

  r = A*x;
  r = -1.0*r;
  y_ax(r, 1.0, b);
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

      y_ax(z, beta, p);
      p = z;
    }

    q = A * p;

    alpha = rho/dot(p,q);

    y_ax(x,  alpha, p);
    y_ax(r, -alpha, q);
    
    if ( nrm2(r)/normr0 < tol ) goto end;
  }
  ret=SOLVERROR_MAXIT;
 end:
  CRSdestory(A);
  return ret;
}
