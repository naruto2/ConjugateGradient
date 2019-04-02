template < class Matrix, class Vector>
int CGS(Matrix &A, Vector &x, Vector &b)
{
  CRSinit(A);

  long i, n=A.size(), maxit=2*n, ret=0;
  double alpha, beta, rho_1, rho_2, normb=nrm2(b), tol=0.000000000001;
  Vector p(n), phat(n), q(n), qhat(n), r(n), rtilde(n), u(n), uhat(n), vhat(n);

  r = b - A*x;

  rtilde = r;
  
  if (normb == 0.0)  normb = 1;
  
  if ( nrm2(r) / normb <= tol) {
    ret = 1; goto end;
  }
  for (i = 0; i <= maxit; i++) {

    rho_1 = dot(rtilde, r);

    if (rho_1 == 0) {
      tol = nrm2(r) / normb;
      ret = 2; goto end;
    }
    
    if (i == 0) {
      u = r;
      p = u;
    } else {
      beta = rho_1 / rho_2;
      u = r; 
      axpy(n, &beta, q, 1, u, 1);
      axpy(n, &beta, p, 1, q, 1);
      p = u;
      axpy(n, &beta, q, 1, p, 1);
    }
    phat = M_solve(p);
    
    vhat = A*phat;
    
    alpha = rho_1 / dot(rtilde, vhat);

    uhat = u;

    alpha = -alpha;

    axpy(n, &alpha, vhat, 1, u, 1);

    alpha = -alpha;

    q = u;
    
    uhat = M_solve(uhat + q);

    axpy(n, &alpha, uhat, 1, x, 1);
    
    qhat = A*uhat;

    alpha = -alpha;

    axpy(n, &alpha, qhat, 1, r, 1);

    rho_2 = rho_1;

    if ( nrm2(r) / normb < tol) {
      goto end;
    }
  }
  ret = 3;
 end:
  CRSdestory(A);
  return ret;
}
