template < class Matrix, class Vector >
void 
Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[])
{
  Vector y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y[i] /= h[i][i];
    for (int j = i - 1; j >= 0; j--)
      y[j] -= h[j][i] * y[i];
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
  printf("hello\n");
  long i, j = 1, k, jj, m = 100, n=A.size(), maxit=20000,
    ret = SOLVERROR_NONE;
  double beta, tol=0.000000000001, normb;
  Vector s(m+1), cs(m+1), sn(m+1), r, w, *v = new Vector[m+1], t;
  Matrix H(n);
  t = b;
  normb = nrm2(M_solveGMRES(t));

  r = M_solveGMRES(y_ax(-1.0*(A*x), 1.0, b));

  beta = nrm2(r);
  printf("world\n");
  if (normb == 0.0) normb = 1;
  printf("nrm2(r)/normb = %f\n",nrm2(r)/normb);
  if (nrm2(r)/normb < tol) goto end;
  printf("earth\n");
  while (j <  maxit) {
    if ( 1  ) if(progress("GMRES",j,nrm2(r)/normb)!=0) {
	ret = 1;
	goto end;
      }
    v[0] = r;
    v[0] = (1.0/beta)*v[0];    

    for (jj=0; jj<(long)s.size(); jj++) s[jj] = 0.0;
    s[0] = beta;
    
    for (i = 0; i < m && j <= maxit; i++, j++) {
      t = A*v[i];
      w = M_solveGMRES(t);
      for (k = 0; k <= i; k++) {
        H[k][i] = dot(w, v[k]);
	y_ax(w, -H[k][i], v[k]);
      }
      H[i+1][i] = nrm2(w);
      v[i+1] = (1.0/H[i+1][i])*w;

      for (k = 0; k < i; k++){
	double csk = cs[k], snk = sn[k];
        ApplyPlaneRotation(H[k][i], H[k+1][i], csk, snk);
	cs[k] = csk, sn[k] = snk;
      }
      
      double csi = cs[i], sni = sn[i], si=s[i], si1=s[i+1];
      GeneratePlaneRotation(H[i][i], H[i+1][i], csi, sni);
      ApplyPlaneRotation(H[i][i], H[i+1][i], csi, sni);
      ApplyPlaneRotation(si, si1, csi, sni);
      cs[i] = csi, sn[i] = sni, s[i] = si, s[i+1] = si1;
      
      si1 = s[i+1];

      if ( abs(si1)/normb < tol) {
        Update(x, i, H, s, v);
	goto end;
      }
    }
    Update(x, m-1, H, s, v);
    t = y_ax(-1.0*(A*x), 1.0, b);
    r = M_solveGMRES(t);

    beta = nrm2(r);

    if ( beta/normb < tol) goto end;
  }
  ret = SOLVERROR_MAXIT;
 end:
  delete [] v;
  CRSdestory(A);
  return ret;
}
