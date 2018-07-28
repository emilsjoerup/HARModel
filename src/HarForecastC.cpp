#include <Rcpparmadillo.h>
using namespace Rcpp;
using namespace arma;



NumericMatrix HARDataCreationCforecast(NumericVector vRealizedmeasure , arma::vec vLags){
  int iT = vRealizedmeasure.size();
  int iLags = vLags.size();
  int iMaxLags = max(vLags);
  
  Rcpp::NumericMatrix  mHarData(iT - max(vLags) , iLags);
  
  for(int i=0; i<iLags; i++){
    for(int j=0; j<iT-max(vLags); j++){
      
      mHarData(j,i) = sum(vRealizedmeasure[Range(iMaxLags+(j) -vLags[i], iMaxLags+(j) -1)])/ vLags[i];
      
      
    }
    
  }
  
  std::cout<< mHarData<<endl;
  return(mHarData);
}




// [[Rcpp::export]]

NumericVector HarForecastLoopingC(NumericVector vRealizedmeasure , NumericVector vLags , NumericVector vCoef, int iNRoll ,int iNAhead , int j){
  int iLags = vLags.size();
  int iMaxLags = max(vLags);
  NumericVector vForecastfoo(iMaxLags + 1);
  NumericVector vHarForecastC(iNAhead);
  NumericVector vLastrowdata(iLags+1);
  NumericMatrix mData = HARDataCreationCforecast(vRealizedmeasure , vLags);
  
    vHarForecastC[0] = vCoef[0] + mData(0,0)*vCoef[1] + mData(1,0)*vCoef[2] + mData(2,0)*vCoef[3];
    for (int i = 1; i <= iMaxLags; ++i)
    {
      
      auto end = std::copy(vRealizedmeasure.end() - iMaxLags + i, vHarForecastC.end(), vHarForecastC.begin());
      std::copy(vHarForecastC.begin(), vHarForecastC.begin() + i, end);
      mData = HARDataCreationCforecast(vHarForecastC , vLags);
      vHarForecastC[i] = vCoef[0] + mData(0,0)*vCoef[1] + mData(1,0)*vCoef[2] + mData(2,0)*vCoef[3];
      //std::cout<< vForecastfoo <<endl;
    }
      
      
  return(vHarForecastC);
}
