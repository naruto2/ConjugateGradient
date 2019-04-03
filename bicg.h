#define trans_mult(A,x) A*x
#define M_trans_solve(x) x

template < class Matrix, class Vector>
int 
BiCG(Matrix &A, Vector &x, const Vector &b)
{
  CRSinit(A);

  long   i, n=A.size(), max_iter=10*n;
  double resid, tol=0.0000000001;
  Vector rho_1(1), rho_2(1), alpha(1), beta(1);
  Vector z, ztilde, p, ptilde, q, qtilde;

  double normb = nrm2(b);
  Vector r = y_ax(-1.0*(A*x), 1.0, b);
  Vector rtilde = r;

  if (normb == 0.0)
    normb = 1;
  
  if ((resid = nrm2(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  for (i = 0; i < max_iter; i++) {

    z = M_solve(r);
    ztilde = M_trans_solve(rtilde);
    rho_1[0] = dot(z, rtilde);
    if (rho_1[0] == 0) { 
      tol = nrm2(r) / normb;
      max_iter = i;
      return 2;
    }
    if (i == 0) {
      p = z;
      ptilde = ztilde;
    } else {
      beta[0] = rho_1[0] / rho_2[0];
      Vector t;
      t = z;
      y_ax(t, beta[0], p);
      p = t;
      t = ztilde;
      y_ax(t, beta[0], ptilde);
      ptilde = t;
    }
    q = A * p;
    qtilde = trans_mult(A,ptilde);
    alpha[0] = rho_1[0] / dot(ptilde, q);

    y_ax(x, alpha[0], p);

    y_ax(r,-alpha[0], q);

    y_ax(rtilde,-alpha[0],qtilde);
    
    rho_2[0] = rho_1[0];
    if ((resid = nrm2(r) / normb) < tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;
}
