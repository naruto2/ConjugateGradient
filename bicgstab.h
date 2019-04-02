template < class Matrix, class Vector>
int 
BiCGSTAB(Matrix &A, Vector &x, Vector &b)
{
  CRSinit(A);

  long n=A.size(), max_iter=2*n;
  double resid, tol=0.000000000001;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1), omega(1);
  Vector p, phat, s, shat, t, v, r(n), rtilde(n);

  double normb = nrm2(b);
  r = b - A*x;
  rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = nrm2(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (int i = 1; i <= max_iter; i++) {
    rho_1[0] = dot(rtilde, r);
    if (rho_1[0] == 0) {
      tol = nrm2(r) / normb;
      return 2;
    }
    if (i == 1)
      p = r;
    else {
      beta[0] = (rho_1[0]/rho_2[0]) * (alpha[0]/omega[0]);

      p = r + beta[0] * (p - omega[0] * v);
    }
    phat = M_solve(p);
    v = A * phat;
    alpha[0] = rho_1[0] / dot(rtilde, v);

    s = r - alpha[0] * v;
    if ((resid = nrm2(s)/normb) < tol) {
      x = x + alpha[0] * phat;
      tol = resid;
      return 0;
    }
    shat = M_solve(s);
    t = A * shat;
    omega[0] = dot(t,s) / dot(t,t);


    Vector x1, x2;
    x1 = alpha[0] * phat;
    x2 = omega[0] * shat;
    x = x + x1;
    x = x + x2;
    //x = x + alpha[0] * phat + omega[0] * shat;


    r = s - omega[0] * t;
    rho_2[0] = rho_1[0];
    if ((resid = nrm2(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
    if (omega[0] == 0) {
      tol = nrm2(r) / normb;
      return 3;
    }
  }

  tol = resid;
  CRSdestory(A);
  return 1;
}
