#define noGPU
#include "/home/tsukud-y/ConjugateGradient/solver.h"

template < typename Real >
class exp_matrix : public matrix < Real > {
public:
  vector<long> ROW, COL, nROW, nCOL;
};


void setval(exp_matrix<double>& A)
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

  A.nROW.resize(A.size());
  A.nCOL.resize(A.size());

  long k;
  for ( k=0; k<A.size(); k++) {
    A.nROW[k] = k;
    A.nCOL[k] = k;
  }
}


void printm(exp_matrix<double>& A)
{
  long i, j,n;
  A.ROW.resize(A.size());
  A.COL.resize(A.size());
  
  
  cout <<"  ";
  for ( i=1; i<A.size(); i++) cout<<" "<<(A.nCOL[i]/10)%10;
  cout<<endl;
  cout <<"  ";
  for ( i=1; i<A.size(); i++) cout<<" "<<A.nCOL[i]%10;
  cout<<endl;

  for (i = 1; i <A.size(); i++){
    printf("%2d ",A.nROW[i]);;
    for(j = 1; j<A.size(); j++) {
      if ( A[i][j] != 0 )cout<<"* ";
      else cout<<"  ";
    }
    cout<<A.ROW[i]<<endl;
  }
  cout<<"  ";
  for (j= 1; j< A.size(); j++){
    cout<<" "<<A.COL[j];
  }
  cout<<endl<<endl;
}

#include "alg_single.h"
