template < class Matrix, class Vector>
int 
BiCGSTAB(Matrix &A, Vector &x, const Vector &b)
{
  CRSinit(A);

  long i, n=A.size(), maxit=10*n, ret=SOLVERROR_NONE;
  double alpha, beta, omega=1.0, rho_1, rho_2=1.0, tol=0.000000000001,
    normb = nrm2(b);
  Vector p, phat, s, shat, t, v, r, rtilde;
  
  r = y_ax(-1.0*(A*x),1.0,b);
  
  rtilde = r;

  if (normb == 0.0) normb = 1;
  
  if (nrm2(r) / normb <= tol)  goto end;

  for (i = 0; i <= maxit; i++) {
    rho_1 = dot(rtilde, r);
    if (rho_1 == 0) {
      ret = SOLVERROR_BREAKDOWN;
      goto end;
    }
    if (i == 0)
      p = r;
    else {
      beta = (rho_1/rho_2) * (alpha/omega);
      p = r;
      t = p;
      y_ax(t,-omega,v);
      y_ax(p, beta, t);
    }
    phat = M_solve(p);

    v = A*phat;

    alpha =  rho_1 / dot(rtilde, v);
    s = r;
    y_ax(s, -alpha, v);

    if (nrm2(s)/normb < tol) {
      y_ax(x, alpha, phat);
      goto end;
    }
    shat = M_solve(s);
    t = A*shat;
    omega = dot(t,s) / dot(t,t);

    y_ax(x,alpha,phat);
    y_ax(x,omega,shat);
    
    r = s;
    y_ax(r, -omega, t);
    
    rho_2 = rho_1;

    if (nrm2(r) / normb < tol)  goto end;

    if (omega == 0) {
      ret = SOLVERROR_BREAKDOWN;
      goto end;
    }
  }
  ret = SOLVERROR_MAXIT;
 end:
  CRSdestory(A);
  return ret;
}
