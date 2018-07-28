#include <Rcpparmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

NumericMatrix HARDataCreationC(NumericVector vRealizedmeasure , arma::vec vLags){
  int iT = vRealizedmeasure.size();
  int iLags = vLags.size();
  int iMaxLags = max(vLags);
  //double dData;
  Rcpp::NumericMatrix  mHarData(iT - max(vLags) , iLags+1);
  
  for(int i=0; i<iLags; i++){
    for(int j=0; j<iT-max(vLags); j++){
    
     
     // std::cout<<vLags[i]<<endl;
      mHarData(j,(i+1)) = sum(vRealizedmeasure[Range(iMaxLags+(j) -vLags[i], iMaxLags+(j) -1)])/ vLags[i];
      
    
    }
    
  }
  mHarData(_,0 ) = vRealizedmeasure[Range(iMaxLags,iT)];
 
  return(mHarData);
}


