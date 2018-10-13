########################################################################################################
# This package is created by Emil Sjoerrup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################



HARestimate = function(vRealizedMeasure , vLags = c(1,5,22), show=TRUE){
  start.time = Sys.time()
  ######Initialization and preparing data ######
  iLags = length(vLags)
  iT = length(vRealizedMeasure)

  mData = HARDataCreationC(vRealizedMeasure, vLags)
  ######Initialization and data preparation end#
  
  ##### Estimate ######
  Model = lm(mData[,1] ~ mData[,2:(iLags+1)])
  ##### Estimation end#
  Info = list("Lags" = vLags)
  names(Model$coefficients) = paste("beta", 0:iLags , sep="")
  #colnames(Model$mVarCovar) = paste("beta", 0:iLags , sep="")
  #rownames(Model$mVarCovar) = paste("beta", 0:iLags , sep="")
  if(is(vRealizedMeasure,"xts")){
    vDates = index(vRealizedMeasure)
    vDates = vDates[(max(vLags)+1):iT]
    Info$Dates = vDates
  }
  Info$ElapsedTime = Sys.time() - start.time
  HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , "Data" = list("Data" = vRealizedMeasure[(max(vLags)+1) : length(vRealizedMeasure)]))
  
  if(show){
    show(HARModel)
  }
  return(HARModel)
  
}


#######Fast estimation routine which returns only the coefficients, it is used for forecasting and cuts down the time a lot.
#######This function is NOT user-callable######
#######iLagsPlusOne is calculated in the forecasting initialization so it can be done once, and not every roll
FASTHARestimate = function(vRealizedMeasure , vLags , iLagsPlusOne){
  mData = HARDataCreationC(vRealizedMeasure , vLags)
  lModel = lm(mData[,1] ~ mData[,2:(iLagsPlusOne)])
  return(lModel$coefficients)
}








