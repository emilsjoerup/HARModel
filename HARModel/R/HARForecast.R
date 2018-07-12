########################################################################################################
# This package is created by Emil Sjørup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################

HARforecast = function( vRealizedmeasure , vLags , iNRoll=1 , iNAhead=10 , bPlots = F , iPlotForecastAhead=1 , sLegendPlacement = "topleft"){
  if(iNRoll > length(vRealizedmeasure)){
    stop("The amount of rolling forecasts cannot be greater than the length of the Realized measure vector.")
  }
  if(class(vRealizedmeasure)[1] == "xts" ) {
   cat("xts is unfortunately not supported yet. The realized measure is passed using the as.vector() function.","\n ")
   vRealizedmeasure = as.vector(vRealizedmeasure)
   } # end xts conditional
  if(any(vRealizedmeasure<0)){
    stop("The realized measure cannot be negative. Something is wrong.")
  } # end negative conditional
  #Initialization
  mForecast = matrix(0 , iNAhead , ncol = iNRoll+1) 
  vObservations =vRealizedmeasure[1:(length(vRealizedmeasure)-iNRoll-1)]  
  vForecastComp = vRealizedmeasure[(length(vRealizedmeasure)-iNRoll):length(vRealizedmeasure)]
  
  iLagsMax = max(vLags)
  
  if(iNAhead ==1 && iNRoll==0){
    #Produces only 1 forecast.
    mData = HARDataCreationC(vRealizedmeasure[1:(length(vRealizedmeasure)-1)] , vLags) # Initialization of data to be used for forecasting
    coef = HARestimate(vRealizedmeasure[1:length(vRealizedmeasure)] , vLags , bPlots =F , bStandardErrors = F)$coef
    mForecast[1,1] = coef[1] + sum(coef[2:length(coef)]*tail(mData,1)[2:length(mData[1,])])
    lModel = HARestimate(vRealizedmeasure[1:(length(vRealizedmeasure)-iNRoll)],vLags)
    lModel$Forecastmatrix = mForecast
    lModel$Observations = vObservations
    lModel$vForeccastComp = vForecastComp
    return(lModel)
  }
  if(iNAhead == 1){
    for (j in 1:(iNRoll+1)) {
      mData = HARDataCreationC(vRealizedmeasure[j:(length(vRealizedmeasure)-iNRoll+j-1)] , vLags) # Initialization of data to be used for forecasting.
      
      # Above makes sure the length of the model stays the same and that for each "roll" the nex period of data is used.
      coef = HARestimate(vRealizedmeasure[j:(length(vRealizedmeasure)-iNRoll+j)] , vLags , bPlots =F , bStandardErrors = F)$coef
      #Extracts the coefficients of the model
      mForecast[1,j] = coef[1] + sum(coef[2:length(coef)]*tail(mData,1)[2:length(mData[1,])])
      #Creates the j'th 1-step ahead forecast
      } # End loop
      }# End one-step ahead conditionals
  else{
    mForecast = matrix(0 , iNAhead+1 , ncol = iNRoll+1) #Initialization
    for (j in 1:(iNRoll+1)) {
      mData = HARDataCreationC(vRealizedmeasure[j:(length(vRealizedmeasure)-iNRoll+j-1)] , vLags)
      coef = HARestimate(vRealizedmeasure[j:(length(vRealizedmeasure)-iNRoll+j)] , vLags , bPlots =F , bStandardErrors = F)$coef
      
      mForecast[1,j] = coef[1] + sum(coef[2:length(coef)]*tail(mData,1)[2:length(mData[1,])])
      
      
      #Creates the j'th 1-step ahead forecast
      for (i in 2:(min(iLagsMax, iNAhead)+1)) {
        
        vForecastfoo = c(vRealizedmeasure[(length(vRealizedmeasure) - max(vLags) - (iNRoll) + (i)): (length(vRealizedmeasure) - (iNRoll))] , mForecast[1:i , j]) ## Fix so the size stays the same.
        
        vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - max(vLags)-1) : length(vForecastfoo)] , vLags)
        #Gets the data necessary to form the forecast in one vector. 
        
        mForecast[i,j] = coef[1] + sum(coef*vLastrowdata)
        
        #Creates the j'th i-step ahead forecast
        
      } #end first nested for-loop
      if(iNAhead>iLagsMax){
      for (i in (iLagsMax+1):(iNAhead+1)) {
        
        vForecastfoo = mForecast[(i-iLagsMax):i , j] ## Fix so the size stays the same.
        
        vLastrowdatafoo = HARDataCreationC(vForecastfoo , vLags)
        #Gets the data necessary to form the forecast in one vector.
        mForecast[i,j] = coef[1] + sum(coef[2:length(coef)]*vLastrowdatafoo[2:length(vLastrowdatafoo)])
        #Creates the j'th i-step ahead forecast
      } #End second nested for-loop
      } #End conditional for the forecast periods being greater than the max of the lag-vector. 
    } #end for-loop
   }# End forecasting
  if(bPlots){
    plot(vForecastComp   , type="l" , ylab = "Realized Volatility" , xlab = "Periods ahead")
    lines(mForecast[iPlotForecastAhead,] , col =2)
    legend(sLegendPlacement , legend = c("Observed (out of sample)" , "Forecasted") , lty = c(1,1)  , col = c(1,2) , box.lty = 0 , inset = c(0.005 , 0.005))
  } # end plot conditional
  lModel = HARestimate(vRealizedmeasure[1:(length(vRealizedmeasure)-iNRoll)],vLags)
  lModel$Forecastmatrix = mForecast
  lModel$Observations = vObservations
  lModel$vForeccastComp = vForecastComp
    return(lModel)
}
