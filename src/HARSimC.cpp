#include <RcppArmadillo.h>
#include "HARDataCreationC.h"
using namespace arma;


// [[Rcpp::export]]

arma::mat HARSimC(int iLength, arma::vec vLags, double dConst , arma::vec vCoef, double dSigma){
  
  int iMaxLags      = max(vLags);
  int iLags         = vLags.size();
  double uncondmean = dConst/(1-sum(vCoef));
  arma::vec vEps    = arma::randn(iLength) * dSigma;
  arma::mat mSim(iLength + 2*iMaxLags , iLags+1); mSim.fill(uncondmean);
  
  for(int i=iMaxLags+1; i<iLength; i++){
    mSim.row(i) = HARDataCreationC(mSim.submat((i-iMaxLags), 0, i, 0), vLags);
    mSim(i,0)   = dConst + sum(mSim(i, span(1,iLags)) * vCoef) + vEps(i);
  }
  
  
  return mSim;
  
}
