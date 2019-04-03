#define H(i,j) H[i][j]
#define h(i,j) h[i][j]

template < class Matrix, class Vector >
void 
Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[])
{
  Vector y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y[j] -= h(j,i) * y[i];
  }

  for (int j = 0; j <= k; j++)
    x = x + y[j] * v[j];
}


template < class Real >
Real 
abs(Real x)
{
  return (x > 0 ? x : -x);
}


#include <math.h> 


template<class Real> 
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
    Real temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    Real temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}


template<class Real>
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}


template < class Matrix, class Vector>
int 
GMRES(Matrix &A, Vector &x, const Vector &b)
{
  CRSinit(A);
  
  long m = 8, n=A.size(), max_iter=2*n;
  double resid, tol=0.0000000001;
  long i, j = 1, k, jj;
  Vector s(m+1), cs(m+1), sn(m+1), w;
  Matrix Hmatrix(n), H = Hmatrix;
  double normb = nrm2(M_solve(b));

  Vector r = A*x;
  r = -1.0*r;
  y_ax(r, 1.0, b);
  r = M_solve(r);
  double beta = nrm2(r);
  
  if (normb == 0.0)
    normb = 1;
  
  if ((resid = nrm2(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  Vector *v = new Vector[m+1];

  while (j <= max_iter) {
    v[0] = r;
    v[0] =  (1.0 / beta)*v[0];    
    for (jj=0; jj<(long)s.size(); jj++) s[jj] = 0.0;
    s[0] = beta;
    
    for (i = 0; i < m && j <= max_iter; i++, j++) {

      w = M_solve(A * v[i]);
      for (k = 0; k <= i; k++) {
        H(k, i) = dot(w, v[k]);
	Vector t = v[k];
        w = w -  H(k, i) * t;
      }
      H(i+1, i) = nrm2(w);
      v[i+1] = (1.0 / H(i+1, i)) * w;

      for (k = 0; k < i; k++){
	double csk = cs[k], snk = sn[k];
        ApplyPlaneRotation(H(k,i), H(k+1,i), csk, snk);
	cs[k] = csk, sn[k] = snk;
      }
      
      {
	double csi = cs[i], sni = sn[i], si=s[i], si1=s[i+1];
	GeneratePlaneRotation(H(i,i), H(i+1,i), csi, sni);
	ApplyPlaneRotation(H(i,i), H(i+1,i), csi, sni);
	ApplyPlaneRotation(si, si1, csi, sni);
	cs[i] = csi, sn[i] = sni, s[i] = si, s[i+1] = si1;
      }
      
      double si1 = s[i+1];
      if ((resid = abs(si1) / normb) < tol) {
        Update(x, i, H, s, v);
        tol = resid;
        max_iter = j;
        delete [] v;
        return 0;
      }
    }
    Update(x, m - 1, H, s, v);


    r = A*x;
    r = -1.0*r;
    y_ax(r, 1.0, b);
    r = M_solve(r);

    beta = nrm2(r);
    if ((resid = beta / normb) < tol) {
      tol = resid;
      max_iter = j;
      delete [] v;
      return 0;
    }
  }
  
  tol = resid;
  delete [] v;
  CRSdestory(A);
  return 1;
}
