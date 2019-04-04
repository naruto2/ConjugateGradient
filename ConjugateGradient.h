template<typename Matrix, typename Vector> 
int ConjugateGradient(Matrix& A, Vector& x, const Vector& b)
{
  CRSinit(A);

  long i, maxit=100000, ret=SOLVERROR_NONE;
  double alpha, beta, normb, rhop, rho, tol=0.00000000001;
  Vector p, q, r, z;  

  r = A*x;
  r = -1.0*r;
  y_ax(r, 1.0, b);
  normb = nrm2(b);
  
  if ( normb == 0.0 ) goto end;
  
  for (i = 0; i < maxit; i++){
    if ( i % 100 == 0 ) if(progress("ConjugateGradient",i,nrm2(r)/normb)!=0) {
	goto end;
      }
    
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
    
    if ( nrm2(r)/normb < tol ) goto end;
  }
  ret=SOLVERROR_MAXIT;
 end:
  CRSdestory(A);
  return ret;
}
