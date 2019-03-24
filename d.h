#ifndef D_H
#define D_H

class Dclass {
  CRSdata CRS;
 public:
  Dclass(CRSdata& CRSparam) {
    CRS = CRSparam;
  }

  void matvec(device_vector<double>& bDev, CRSdata& CRS,
	      device_vector<double>& xnDev){

    y_Ax(dev(bDev), CRS, dev(xnDev));
  }
  
  void scal(long N, double* alphap, device_vector<double>& rdev, long l){
    cublasDscal(CRS.cublas, N, alphap, dev(rdev), l);
  }

  void axpy(long N, double* alphap, device_vector<double>& bDev, long l,
	    device_vector<double>& rDev, long I){
    cublasDaxpy(CRS.cublas, N, alphap, dev(bDev), l, dev(rDev), I);
  }
  

};


#endif
