#include <RcppArmadillo.h>

using namespace arma;

//[[Rcpp::export]]
arma::mat HARDataCreationC(arma::vec vRealizedMeasure , arma::vec vPeriods, int h = 1){
  int iT       = vRealizedMeasure.size();
  int iLags    = vPeriods.size();
  int iMaxLags = max(vPeriods);
  arma::mat mHARData((iT - iMaxLags-h+1) , (iLags + 1));
  int k = 0;
  if(vPeriods[0] == 1){
    mHARData.col(1) = vRealizedMeasure(arma::span(iMaxLags-1 , (iT - h-1)));
    k = 1;
  }
  
  for(int i = k; i<iLags; i++){
    for(int j = 0; j<(iT - iMaxLags-h+1); j++){
      mHARData(j,(i+1)) = sum(vRealizedMeasure(arma::span((iMaxLags + j - vPeriods[i]) , (iMaxLags + j - 1) )))/vPeriods[i];
    }
  }
  
  if(h != 1){ //need to aggregate the dependent variable too.
    for(int j = 0; j<=(iT - iMaxLags-h); j++){
      mHARData(j,0) = sum(vRealizedMeasure(arma::span((iMaxLags + j) , (iMaxLags + j+h-1))))/h;
    }
    
  }
  else{
    mHARData.col(0) = vRealizedMeasure(arma::span((iMaxLags) , (iT-1)));
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
