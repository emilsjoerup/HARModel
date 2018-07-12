########################################################################################################
# This package is created by Emil Sjørup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################



HARestimate = function(vRealizedmeasure , vLags , bStandardErrors=T , bPlots=F , iLagSE = 5 , sLegendPlacement = "topleft"){
  mData = HARDataCreationC(vRealizedmeasure , vLags )
  lModel = lm(mData[,1] ~ mData[,2:length(mData[1,])])
  if(bPlots){
    dates = index(vRealizedmeasure)
    dates = dates[(max(vLags)+1):length(vRealizedmeasure)]
    vFittedValues = xts(lModel$fitted.values , order.by = dates)
    vObservedRM = xts(lModel$model[,1] , order.by=dates)
    print(plot(cbind(vObservedRM , vFittedValues ) ,col = c(1,2) , main= "Observed vs. Fitted values" , lty=c(1,1) , yaxis.same =T , grid.ticks.on=1))
    print(addLegend(sLegendPlacement, on=1, legend.names= c("Observed" , "Fitted values") , col = c(1,2) , lty=c(1,1) , lwd=c(2,2)))
    
  } # End conditional for plots.
  if(bStandardErrors){
  mVarCovar = sandwich::NeweyWest(lModel , lag = iLagSE)
  lModel$mVarCovar = mVarCovar
  return(lModel)
  } # End conditional for Standard errors

  return(lModel)
  
}