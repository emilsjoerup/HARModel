#include <RcppArmadillo.h>

using namespace arma;
//[[Rcpp::export]]
arma::mat HARDataCreationC(arma::vec vRealizedMeasure , arma::vec vLags){
  int iT       = vRealizedMeasure.size();
  int iLags    = vLags.size();
  int iMaxLags = max(vLags);
  arma::mat mHARData((iT - iMaxLags) , (iLags + 1));
  int k = 0;
  
  mHARData.col(0) = vRealizedMeasure(arma::span((iMaxLags) , (iT-1)));
  
  if(vLags[0] == 1){
    mHARData.col(1) = vRealizedMeasure(arma::span(iMaxLags-1 , (iT - 2)));
    k = 1;
  }
  
  for(int i = k; i<iLags; i++){
    for(int j = 0; j<(iT - iMaxLags); j++){
      mHARData(j,(i+1)) = sum(vRealizedMeasure(arma::span((iMaxLags + j - vLags[i]) , (iMaxLags + j - 1) )))/vLags[i];
    }
    
  }
  return(mHARData);
}

//[[Rcpp::export]]
arma::mat HARMatCombine(arma::mat mA, arma::mat mB){
  int iRowsA = mA.n_rows;
  int iRowsB = mB.n_rows;
  arma:: vec vFoo(2); vFoo << iRowsA << iRowsB <<endr;
  arma::mat mC((arma::min(vFoo)), (mA.n_cols+mB.n_cols));
  if(iRowsA > iRowsB){
    mA.shed_rows(0, iRowsA-iRowsB-1);
  }
  if(iRowsA < iRowsB){
    mB.shed_rows(0, iRowsB-iRowsA-1);
  }
  
  mC = join_rows(mA, mB);
  
  return(mC);
  
  
  
}
