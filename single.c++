#define noGPU
#include "/home/tsukud-y/ConjugateGradient/solver.h"

typedef struct {
  matrix<double> A;
  vector<long> ROW, COL, nROW, nCOL;
} EnhancedMatrix;



void setvalE(EnhancedMatrix& D)
{
  D.A.resize(17);

  matrix<double> A =D.A;
  
  A[1][4] = A[1][5] = A[1][8] = 1;
  A[2][1] = A[2][5] = A[2][7] = A[2][8] = A[2][11] = A[2][16] = 1;
  A[3][7] = 1;
  A[4][2] = A[4][5] = A[4][6] = A[4][7] = A[4][9] = A[4][15] = A[4][16] = 1;
  A[5][2] = A[5][4] = A[5][5] = 1;
  A[6][5] = A[6][10]= A[6][13]= A[6][15]= 1;
  A[7][4] = A[7][10]= 1;
  A[8][4] = A[8][10]= 1;
  A[9][3] = A[9][5] = A[9][6] = A[9][8] = A[9][12] = 1;
  A[10][1]= A[10][2]= A[10][5]= A[10][16]=1;
  A[11][1]= A[11][3]= A[11][12]=1;
  A[12][3]= A[12][6]= A[12][12]=A[12][15]=1;
  A[13][1]= A[13][16]=1;
  A[14][7]= A[14][10]=A[14][13]=1;
  A[15][1]= A[15][9] =A[15][10]=A[15][13]=A[15][14]=A[15][15]=1;
  A[16][1]= A[16][2] =A[16][8] =A[16][11]=1;

  D.nROW.resize(A.size());
  D.nCOL.resize(A.size());

  long k;
  for ( k=0; k<A.size(); k++) {
    D.nROW[k] = k;
    D.nCOL[k] = k;
  }
  D.A = A;
}


void setval(matrix<double>& A, vector<long>&nROW, vector<long>&nCOL)
{
  A.resize(17);
  
  A[1][4] = A[1][5] = A[1][8] = 1;
  A[2][1] = A[2][5] = A[2][7] = A[2][8] = A[2][11] = A[2][16] = 1;
  A[3][7] = 1;
  A[4][2] = A[4][5] = A[4][6] = A[4][7] = A[4][9] = A[4][15] = A[4][16] = 1;
  A[5][2] = A[5][4] = A[5][5] = 1;
  A[6][5] = A[6][10]= A[6][13]= A[6][15]= 1;
  A[7][4] = A[7][10]= 1;
  A[8][4] = A[8][10]= 1;
  A[9][3] = A[9][5] = A[9][6] = A[9][8] = A[9][12] = 1;
  A[10][1]= A[10][2]= A[10][5]= A[10][16]=1;
  A[11][1]= A[11][3]= A[11][12]=1;
  A[12][3]= A[12][6]= A[12][12]=A[12][15]=1;
  A[13][1]= A[13][16]=1;
  A[14][7]= A[14][10]=A[14][13]=1;
  A[15][1]= A[15][9] =A[15][10]=A[15][13]=A[15][14]=A[15][15]=1;
  A[16][1]= A[16][2] =A[16][8] =A[16][11]=1;

  nROW.resize(A.size());
  nCOL.resize(A.size());

  long k;
  for ( k=0; k<A.size(); k++) {
    nROW[k] = k;
    nCOL[k] = k;
  }
}


void printm(matrix<double>& A,
	    vector<long>& ROW, vector<long>& COL,
	    vector<long>&nROW, vector<long>&nCOL)
{
  long i, j,n;
  ROW.resize(A.size());
  COL.resize(A.size());
  
  
  cout <<"  ";
  for ( i=1; i<A.size(); i++) cout<<" "<<(nCOL[i]/10)%10;
  cout<<endl;
  cout <<"  ";
  for ( i=1; i<A.size(); i++) cout<<" "<<nCOL[i]%10;
  cout<<endl;

  for (i = 1; i <A.size(); i++){
    printf("%2d ",nROW[i]);;
    for(j = 1; j<A.size(); j++) {
      if ( A[i][j] != 0 )cout<<"* ";
      else cout<<"  ";
    }
    cout<<ROW[i]<<endl;
  }
  cout<<"  ";
  for (j= 1; j< A.size(); j++){
    cout<<" "<<COL[j];
  }
  cout<<endl<<endl;
}

void printmE(EnhancedMatrix& D)
{
  matrix<double> A = D.A;
  long i, j,n;
  D.ROW.resize(A.size());
  D.COL.resize(A.size());
  
  
  cout <<"  ";
  for ( i=1; i<A.size(); i++) cout<<" "<<(D.nCOL[i]/10)%10;
  cout<<endl;
  cout <<"  ";
  for ( i=1; i<A.size(); i++) cout<<" "<<D.nCOL[i]%10;
  cout<<endl;

  for (i = 1; i <A.size(); i++){
    printf("%2d ",D.nROW[i]);;
    for(j = 1; j<A.size(); j++) {
      if ( A[i][j] != 0 )cout<<"* ";
      else cout<<"  ";
    }
    cout<<D.ROW[i]<<endl;
  }
  cout<<"  ";
  for (j= 1; j< A.size(); j++){
    cout<<" "<<D.COL[j];
  }
  cout<<endl<<endl;
}


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



void procedure0(matrix<double>&A, vector<long>& ROW, vector<long>& COL)
{
  long i, j, c;
  ROW.resize(A.size());
  COL.resize(A.size());
  for (i=1; i<A.size(); i++){
    ROW[i] = 0;
    for (j=1; j<A.size(); j++) if( A[i][j] != 0) ROW[i]++;
  }
  for (j=1; j<A.size(); j++){
    COL[j] = 0;
    for (i=1; i<A.size(); i++) if( A[i][j] != 0) COL[j]++;
  }
}

void procedure0E(EnhancedMatrix& D)
{
  long i, j, c;
  D.ROW.resize(D.A.size());
  D.COL.resize(D.A.size());
  for (i=1; i<D.A.size(); i++){
    D.ROW[i] = 0;
    for (j=1; j<D.A.size(); j++) if( D.A[i][j] != 0) D.ROW[i]++;
  }
  for (j=1; j<D.A.size(); j++){
    D.COL[j] = 0;
    for (i=1; i<D.A.size(); i++) if( D.A[i][j] != 0) D.COL[j]++;
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

void procedure1E(EnhancedMatrix& D)
{
  long i, j, k, m = D.A.size()-1;
  long s, t;  
  for (k=1; k<=m; k++){
    if(p1step1(D.A,D.ROW,s,t)!=0) break;
    p1step2(D.A,D.ROW,t);
    printf("s,t,k= %d, %d, %d\n",s,t,k);
    shiftCOL(k,t,D.A,D.ROW,D.COL,D.nROW,D.nCOL);
    shiftROW(k,s,D.A,D.ROW,D.COL,D.nROW,D.nCOL);
    return;
  }

}
void procedure1(matrix<double>&A,
		vector<long>& ROW, vector<long>& COL,
		vector<long>&nROW, vector<long>&nCOL)
{
  long i, j, k, m = A.size()-1;
  long s, t;
  
  for (k=1; k<=m; k++){
    if(p1step1(A,ROW,s,t)!=0) break;
    p1step2(A,ROW,t);
    printf("s,t= %d, %d\n",s,t);
    shiftCOL(k,t,A,ROW,COL,nROW,nCOL);
    shiftROW(k,s,A,ROW,COL,nROW,nCOL);
  }

}


void procedure2(matrix<double>&A,
		vector<long>& ROW, vector<long>& COL,
		vector<long>&nROW, vector<long>&nCOL)
{
  long I, i, j, k, m = A.size()-1, u, v;

  I=1;
  while(!COL[I++]);
  cout<<"I= "<<I<<endl;

 start:
  for(j=I; j<A.size(); j++) {
    if ( COL[j] == 1 ){
      v = j;
      for(i=1;i<A.size();i++) {
	if (A[i][j] == 0 ) continue;
	break;
      }
      u = i;
    } else continue;
    printf("u,v,j = %d,%d,%d\n",u,v,j);
    printm(A,ROW,COL,nROW,nCOL);

    shiftCOL(I,v,A,ROW,COL,nROW,nCOL);
    exit(0);
    swapROW(I,u,A,ROW,COL,nROW,nCOL);
    printm(A,ROW,COL,nROW,nCOL);
    I++;
    goto start;
  }

  
}

void procedure2E(EnhancedMatrix& D)
{
  matrix<double> A = D.A;
  vector<long> ROW=D.ROW, COL=D.COL, nROW=D.nROW, nCOL=D.nCOL;
  long I, i, ii, j, k, m = A.size()-1, u, v;
  I=1;
  while(!COL[I++]);
  cout<<"I= "<<I<<endl;

 start:
  printf("start I= %d\n",I);
  for(j=I; j<A.size(); j++) {
    printf("for j= %d\n",j);
    if ( COL[j] == 1 ){
      v = j;
      for(i=I;i<A.size();i++) {
	if (A[i][j] == 0 ) continue;
	break;
      }
      u = i;
    } else continue;
    printf("u,v,j = %d,%d,%d\n",u,v,j);
    shiftCOL(I,v,A,ROW,COL,nROW,nCOL);
    shiftROW(I,u,A,ROW,COL,nROW,nCOL);
    I++;
    printf("I++ = %d\n",I);
    for(ii = I; ii <A.size(); ii++) if(A[I-1][ii] !=0){
	COL[ii]--;
	printf("COL = %d\n",ii);
      }
    printm(A,ROW,COL,nROW,nCOL);
    goto start;
  }
 end:
  D.A = A, D.ROW=ROW, D.COL=COL, D.nROW=nROW, D.nCOL=nCOL;
}
  

int main()
{
  EnhancedMatrix D;
  
  setvalE(D);

  procedure0E(D);

  procedure1E(D);
  printmE(D);
  procedure2E(D);
  printmE(D);


}
