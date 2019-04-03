#include <math.h>
#define M1_solve(r) r
#define M1_trans_solve(r) r
#define M2_solve(r) r
#define M2_trans_solve(r) r

template < class Matrix, class Vector>
  int
  QMR(Matrix &A, Vector &x, const Vector &b)
{
  CRSinit(A);
  long   n=A.size(), max_iter = 10*n;
  double resid, tol=0.00000000001;
  
  Vector rho(1), rho_1(1), xi(1), gamma(1), gamma_1(1), theta(1), theta_1(1);
  Vector eta(1), delta(1), ep(1), beta(1);

  Vector r, v_tld, y, w_tld, z;
  Vector v, w, y_tld, z_tld;
  Vector p, q, p_tld, d, s, t;

  double normb = nrm2(b);

  r = y_ax(-1.0*(A*x), 1.0, b);

  if (normb == 0.0)
    normb = 1;

  if ((resid = nrm2(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  v_tld = r;
  y = M1_solve(v_tld);
  rho[0] = nrm2(y);

  w_tld = r;
  z = M2_trans_solve(w_tld);
  xi[0] = nrm2(z);

  gamma[0] = 1.0;
  eta[0] = -1.0;
  theta[0] = 0.0;

  for (int i = 1; i <= max_iter; i++) {

    if (rho[0] == 0.0)
      return 2;                        // return on breakdown

    if (xi[0] == 0.0)
      return 7;                        // return on breakdown
    v.clear();
    v.resize(n);
    y_ax(v,(1./rho[0]), v_tld);
    y = (1. / rho[0]) * y;

    t = w_tld;
    w = (1./xi[0]) * t;
    z = (1./xi[0]) * z;

    delta[0] = dot(z, y);
    if (delta[0] == 0.0)
      return 5;                        // return on breakdown

    y_tld = M2_solve(y);               // apply preconditioners
    z_tld = M1_trans_solve(z);

    if (i > 1) {
      t = y_tld;
      y_ax(t, -(xi[0]*delta[0]/ep[0]), p);
      p = t;
      t = z_tld;
      y_ax(t, -(rho[0]*delta[0]/ep[0]), q);
      q = t;
    } else {
      p = y_tld;
      q = z_tld;
    }

    p_tld = A * p;
    ep[0] = dot(q, p_tld);
    if (ep[0] == 0.0)
      return 6;                        // return on breakdown

    beta[0] = ep[0] / delta[0];
    if (beta[0] == 0.0)
      return 3;                        // return on breakdown
    t = p_tld;
    v_tld = y_ax(t, -beta[0], v);

    y = M1_solve(v_tld);

    rho_1[0] = rho[0];
    rho[0] = nrm2(y);

    t = trans_mult(A,q);
    w_tld = y_ax(t, -beta[0], w);

    z = M2_trans_solve(w_tld);

    xi[0] = nrm2(z);

    gamma_1[0] = gamma[0];
    theta_1[0] = theta[0];

    theta[0] = rho[0] / (gamma_1[0] * beta[0]);
    gamma[0] = 1.0 / sqrt(1.0 + theta[0] * theta[0]);

    if (gamma[0] == 0.0)
      return 4;                        // return on breakdown

    eta[0] = -eta[0] * rho_1[0] * gamma[0] * gamma[0] / 
      (beta[0] * gamma_1[0] * gamma_1[0]);

    if (i > 1) {
      Vector d1, d2, s1, s2;
      t = p;
      d1 = eta[0] * t;
      d2 = (theta_1[0] * theta_1[0] * gamma[0] * gamma[0]) * d;
      d = d1 + d2;
      s1 = eta[0] * p_tld;
      s2 = (theta_1[0] * theta_1[0] * gamma[0] * gamma[0]) * s;
      s = s1 + s2;
    } else {
      d = eta[0] * p;
      s = eta[0] * p_tld;
    }

    x = x + d;                            // update approximation vector
    r = r - s;                            // compute residual

    if ((resid = nrm2(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;                            // no convergence
}
