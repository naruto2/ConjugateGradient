template < class Matrix, class Vector>
int 
BiCGSTAB(Matrix &A, Vector &x, const Vector &b)
{
  CRSinit(A);

  long i, n=A.size(), max_iter=2*n;
  double resid, tol=0.000000000001;
  double rho_1, rho_2=1.0, alpha, beta, omega=1.0;
  Vector p, phat, s, shat, t, v, r, rtilde;
  double normb = nrm2(b);

  r = A*x;
  r = -1.0*r;
  alpha = 1.0;
  axpy(n, &alpha, b, 1, r, 1);

  rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = nrm2(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (i = 0; i <= max_iter; i++) {
    rho_1 = dot(rtilde, r);
    if (rho_1 == 0) {
      tol = nrm2(r) / normb;
      return 2;
    }
    if (i == 0)
      p = r;
    else {
      beta = (rho_1/rho_2) * (alpha/omega);

      omega = -omega;
      p = r;
      t = p;
      axpy(n, &omega, v, 1, t, 1);
      axpy(n, &beta,  t, 1, p, 1);
      
    }
    phat = M_solve(p);

    v = A * phat;

    alpha = - rho_1 / dot(rtilde, v);
    s = r;
    axpy(n, &alpha, v, 1, s, 1);
    alpha = -alpha;

    if ((resid = nrm2(s)/normb) < tol) {
      x = x + alpha * phat;
      tol = resid;
      return 0;
    }
    shat = M_solve(s);
    t = A * shat;
    omega = dot(t,s) / dot(t,t);

    axpy(n, &alpha, phat, 1, x, 1);

    axpy(n, &omega, shat, 1, x, 1);
    r = s;
    omega = -omega;
    axpy(n, &omega, t, 1, r, 1);
    
    rho_2 = rho_1;
    if ((resid = nrm2(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega == 0) {
      tol = nrm2(r) / normb;
      return 3;
    }
  }

  tol = resid;
  CRSdestory(A);
  return 1;
}
