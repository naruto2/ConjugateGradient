#define noGPU
#include "/home/tsukud-y/ConjugateGradient/solver.h"

void setval(matrix<double>& A)
{
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
}


void printm(matrix<double>& A)
{
  long i, j;
  printf("                     1 1 1 1 1 1 1 \n");
  printf("   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 \n");
  
  for (i = 1; i <A.size(); i++){
    printf("%2d ",i);;
    for(j = 1; j<A.size(); j++) {
      if ( A[i][j] != 0 )cout<<"* ";
      else cout<<"  ";
    }
    cout<<endl;
  }
}

void printm0(matrix<double>& A,vector<double>& ROW, vector<double>& COL)
{
  long i, j;
  printf("                     1 1 1 1 1 1 1 \n");
  printf("   1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 \n");
  
  for (i = 1; i <A.size(); i++){
    printf("%2d ",i);;
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
  cout<<endl;
}


void procedure0(matrix<double>&A, vector<double>& ROW, vector<double>& COL)
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

int p1step1(matrix<double>& A, vector<double>& ROW,long& s,long& t)
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


void p1step2(matrix<double>& A, vector<double>& ROW, long& t)
{
  long i, m= A.size()-1;
  for(i=1; i<=m; i++) if(A[i][t] !=0) ROW[i]--;
}


void procedure1(matrix<double>&A, vector<double>& ROW, vector<double>& COL)
{
  long i, j, k, m = A.size()-1;
  long s, t;

  for (k=1; k<=m; k++){
    if(p1step1(A,ROW,s,t)!=0) break;
    printf("s,t= %d, %d\n",s,t);
    p1step2(A,ROW,t);
    printm0(A,ROW,COL);
  }
}

void procedure2(matrix<double>&A, vector<double>& ROW, vector<double>& COL)
{
  long i, j, k, m = A.size()-1, u, v;

  for(u=0, k=1; k<=m; k++) {
    if (ROW[k] == 0) continue;
    if (COL[k] == 1) v=k;
  }
  if (u!=0) {
    j = v;
    for (i=1; i<m; i++){
    }
  }

}
  
int main()
{
  matrix<double> A(17);
  vector<double> ROW,COL;
  
  setval(A);
  printm(A);

  procedure0(A,ROW,COL);
  printm0(A,ROW,COL);

  procedure1(A,ROW,COL);

  procedure2(A,ROW,COL);
}
