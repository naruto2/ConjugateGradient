template < class Matrix, class Vector>
int CGS(Matrix &A, Vector &x, const Vector &b)
{
  CRSinit(A);

  long i, n=A.size(), maxit=2*n, ret=SOLVERROR_NONE;
  double alpha, beta, rho_1, rho_2, normb=nrm2(b), tol=0.000000000001;
  Vector p, phat, q, qhat, r, rtilde, u, uhat, vhat;

  r = y_ax(-1.0*(A*x), 1.0, b);
  rtilde = r;
  
  if (normb == 0.0)  normb = 1;
  
  if ( nrm2(r) / normb <= tol)  goto end;

  for (i = 0; i <= maxit; i++) {
    if ( i % 100 == 0 ) if(progress("CGS",i,nrm2(r)/normb)!=0) {
	ret = 1;
	goto end;
      }
    rho_1 = dot(rtilde, r);

    if (rho_1 == 0) {
      ret = SOLVERROR_BREAKDOWN; goto end;
    }
    
    if (i == 0) {
      u = r;
      p = u;
    } else {
      beta = rho_1 / rho_2;
      u = r; 
      y_ax(u, beta, q);
      y_ax(q, beta, p);
      p = u;
      y_ax(p, beta, q);
    }
    phat = M_solve(p);
    
    vhat = A*phat;
    
    alpha = rho_1 / dot(rtilde, vhat);

    uhat = u;

    y_ax(u, -alpha, vhat);
    
    q = u;
    
    uhat = M_solve(uhat + q);

    y_ax(x, alpha, uhat);
    
    qhat = A*uhat;

    y_ax(r, -alpha, qhat);
    
    rho_2 = rho_1;

    if ( nrm2(r) / normb < tol) goto end;
  }
  ret = SOLVERROR_MAXIT;
 end:
  CRSdestory(A);
  return ret;
}
