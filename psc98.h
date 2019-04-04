extern "C" {
  void genmat(int,int*,double*,double*);
  void chkval(FILE*,int,double*);
}

template < class Matrix, class Vector >
void PSC98init(Matrix& A, Vector& x, Vector& b)
{
  vector<double> AA(10);
  double B;
  int    JA[10], i, j, n, w;

  genmat(-1,&JA[0],&AA[0],&B);
  n = JA[0]; w = JA[2];

  A.resize(n);  x.resize(n); b.resize(n); 

  for ( i=1; i<=n; i++) {
    for (j=0;j<=w-1;j++) {
      JA[j] =  -1;
      AA[j] = 0.0;
    }
    genmat(i,&JA[0],&AA[0],&B);
    for ( j=0; j<w; j++) if (JA[j] != -1){
	if( JA[j] <=0 || n < JA[j] ) {
	  ;
	} else {
	  A[i-1][JA[j]-1] = AA[j];
	  b[i-1] = B;
	}
      }
  }
}


template < class Vector>
void PSC98check(Vector& x)
{
  vector<double> x1(x.size());
  for (unsigned int i=0; i<x.size(); i++) x1[i] = x[i];
  chkval(stdout,x1.size(),&x1[0]);
}