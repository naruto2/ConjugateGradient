namespace single {
  
void swapCOL(long s0, long s1, matrix<double>& A,
	     vector<long>& ROW, vector<long>& COL,
	     vector<long>&nROW, vector<long>&nCOL)
{
  double x;
  long i;

  for ( i=0; i<A.size(); i++){
    x = A[i][s0];
    A[i][s0] = A[i][s1];
    A[i][s1] = x;
  }
  x = COL[s0];
  COL[s0] = COL[s1];
  COL[s1] = x;
  
  x = nCOL[s0];
  nCOL[s0] = nCOL[s1];
  nCOL[s1] = x;
}


void swapROW(long s0, long s1, matrix<double>& A,
	     vector<long>& ROW, vector<long>& COL,
	     vector<long>&nROW, vector<long>&nCOL)
{
  double x;
  long j;

  for ( j=0; j<A.size(); j++){
    x = A[s0][j];
    A[s0][j] = A[s1][j];
    A[s1][j] = x;
  }
  x = ROW[s0];
  ROW[s0] = ROW[s1];
  ROW[s1] = x;
  
  x = nROW[s0];
  nROW[s0] = nROW[s1];
  nROW[s1] = x;
}


void procedure0(exp_matrix<double>& A)
{
  long i, j, c;
  A.ROW.resize(A.size());
  A.COL.resize(A.size());
  for (i=1; i<A.size(); i++){
    A.ROW[i] = 0;
    for (j=1; j<A.size(); j++) if( A[i][j] != 0) A.ROW[i]++;
  }
  for (j=1; j<A.size(); j++){
    A.COL[j] = 0;
    for (i=1; i<A.size(); i++) if( A[i][j] != 0) A.COL[j]++;
  }
}


int p1step1(matrix<double>& A, vector<long>& ROW,long& s,long& t)
{
  long i,j, k, m = A.size()-1;
  for(s=0, k=1; k<=m; k++) {
    if(ROW[k] == 1) s = k;
  }
  if ( s != 0 ){
    i = s;
    for (j=1; j<m; j++) if (A[i][j] !=0) t =j;
  } else return 1;
  return 0;
}


void p1step2(matrix<double>& A, vector<long>& ROW, long& t)
{
  long i, m= A.size()-1;
  for(i=1; i<=m; i++) if(A[i][t] !=0) ROW[i]--;
}


void shiftCOL(long I, long v, matrix<double>& A,
	 vector<long>& ROW, vector<long>& COL,
	 vector<long>&nROW, vector<long>&nCOL)
{
  for(;I<v;I++) swapCOL(I,v,A,ROW,COL,nROW,nCOL);
}

void shiftROW(long I, long v, matrix<double>& A,
	 vector<long>& ROW, vector<long>& COL,
	 vector<long>&nROW, vector<long>&nCOL)
{
  for(;I<v;I++) swapROW(I,v,A,ROW,COL,nROW,nCOL);
}

void procedure1(exp_matrix<double>& A)
{
  long i, j, k, m = A.size()-1;
  long s, t;  
  for (k=1; k<=m; k++){
    if(p1step1(A,A.ROW,s,t)!=0) break;
    p1step2(A,A.ROW,t);
    shiftCOL(k,t,A,A.ROW,A.COL,A.nROW,A.nCOL);
    shiftROW(k,s,A,A.ROW,A.COL,A.nROW,A.nCOL);
    return;
  }

}



void procedure2(exp_matrix<double>& A)
{
  vector<long> ROW=A.ROW, COL=A.COL, nROW=A.nROW, nCOL=A.nCOL;
  long I, i, ii, j, k, m = A.size()-1, u, v;
  I=1;
  while(!COL[I++]);

 start:
  for(j=I; j<A.size(); j++) {
    if ( COL[j] == 1 ){
      v = j;
      for(i=I;i<A.size();i++) {
	if (A[i][j] == 0 ) continue;
	break;
      }
      u = i;
    } else continue;
    shiftCOL(I,v,A,ROW,COL,nROW,nCOL);
    shiftROW(I,u,A,ROW,COL,nROW,nCOL);
    I++;
    for(ii = I; ii <A.size(); ii++) if(A[I-1][ii] !=0){
	COL[ii]--;
      }
    goto start;
  }
 end:
  A.ROW=ROW, A.COL=COL, A.nROW=nROW, A.nCOL=nCOL;
}
}

void alg_single(exp_matrix<double>& A)
{
  single::procedure0(A);
  single::procedure1(A);
  single::procedure2(A);
}    

