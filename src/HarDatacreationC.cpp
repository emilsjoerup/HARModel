#include <RcppArmadillo.h>

using namespace arma;
//[[Rcpp::export]]
arma::mat HARDataCreationC(arma::vec vRealizedMeasure , arma::vec vLags){
  int iT       = vRealizedMeasure.size();
  int iLags    = vLags.size();
  int iMaxLags = max(vLags);
  arma::mat mHARData((iT - iMaxLags) , (iLags + 1));
  int i = 0;
  
  mHARData.col(0) = vRealizedMeasure(arma::span((iMaxLags) , (iT-1)));

  if(vLags[0] == 1){
    mHARData.col(1) = vRealizedMeasure(arma::span(iMaxLags-1 , (iT - 2)));
    i = 1;
  }

  for(i = i; i<iLags; i++){
  for(int j = 0; j<(iT - iMaxLags); j++){
    mHARData(j,(i+1)) = sum(vRealizedMeasure(arma::span((iMaxLags + j - vLags[i]) , (iMaxLags + j - 1) )))/vLags[i];
  }

  }
  return(mHARData);
}

