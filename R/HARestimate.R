########################################################################################################
# This package is created by Emil Sjørup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################



HARestimate = function(vRealizedmeasure , vLags , bStandardErrors=T ,iLagSE = 5){
  mData = HARDataCreationC(vRealizedmeasure , vLags )
  lModel = lm(mData[,1] ~ mData[,2:length(mData[1,])])
  vDates = index(vRealizedmeasure)
  vDates = vDates[(max(vLags)+1):length(vRealizedmeasure)]
  if(bStandardErrors){
  mVarCovar = sandwich::NeweyWest(lModel , lag = iLagSE)
  lModel$mVarCovar = mVarCovar
  lModel$NWLagOrder = iLagSE
  names(lModel$coefficients) = paste("beta", 0:length(vLags) , sep="")
  lModel$terms = list("RV" , "=", "c", paste("RV" , vLags[1:length(vLags)] , sep=""))
  lModel$dates = vDates
  class(lModel) = c("HARmodel" , "lm")
  return(lModel)
  } # End conditional for Standard errors
  names(lModel$coefficients) = paste("beta", 0:length(vLags) , sep="")
  lModel$dates = vDates
  lModel$terms = list("RV" , "=", "c", paste("RV" , vLags[1:length(vLags)] , sep=""))
  return(lModel)
  
}