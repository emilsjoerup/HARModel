########################################################################################################
# This package is created by Emil Sjørup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################



HARestimate = function(vRealizedMeasure , vLags = c(1,5,22),iLagSE = 5){
  ######Initialization and preparing data ######
  iLags = length(vLags)
  iT = length(vRealizedMeasure)
  vDates = index(vRealizedMeasure)
  vDates = vDates[(max(vLags)+1):iT]
  mData = HARDataCreationC(vRealizedMeasure, vLags)
  ######Initialization and data preparation end#
  
  ##### Estimate ######
  lModel = lm(mData[,1] ~ mData[,2:(iLags+1)])
  
  mVarCovar = sandwich::NeweyWest(lModel , lag = iLagSE)
  ##### Estimation end#
  
  lModel$mVarCovar = mVarCovar
  lModel$NWLagOrder = iLagSE
  names(lModel$coefficients) = paste("beta", 0:iLags , sep="")
  lModel$terms = list("RV" , "=", "c", paste("RV" , vLags[1:iLags] , sep=""))
  lModel$dates = vDates
  class(lModel) = c("HARmodel" , "lm")
  return(lModel)
  
}


#######Fast estimation routine which returns only the coefficients, it is used for forecasting and cuts down the time a bit.
#######This function is NOT user-callable######
#######iLagsPlusOne is calculated in the forecasting initialization so it can be done once isn't done every roll.
FASTHARestimate = function(vRealizedMeasure , vLags , iLagsPlusOne){
  mData = HARDataCreationC(vRealizedMeasure , vLags)
  lModel = lm(mData[,1] ~ mData[,2:(iLagsPlusOne)])
  return(lModel$coefficients)
  
}