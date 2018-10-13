#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::export]]

NumericMatrix HARDataCreationC(NumericVector vRealizedmeasure , NumericVector vLags){
  int iT = vRealizedmeasure.size();
  int iLags = vLags.size();
  int iMaxLags = max(vLags);
  NumericMatrix  mHarData(iT - max(vLags) , iLags+1);
  
  for(int i=0; i<iLags; i++){
    if(vLags[i] == 1){
      mHarData(_,(i+1)) = vRealizedmeasure[Range((iMaxLags-1) , (iT-1))];
      
    }
    for(int j=0; j<iT-max(vLags); j++){
      mHarData(j,(i+1)) = sum(vRealizedmeasure[Range(iMaxLags+(j) -vLags[i], iMaxLags+(j) -1)])/ vLags[i];
      
    }
    
  }
  mHarData(_,0 ) = vRealizedmeasure[Range(iMaxLags,iT)];
 
  return(mHarData);
}


