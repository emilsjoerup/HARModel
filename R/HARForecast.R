########################################################################################################
# This package is created by Emil Sjørup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################

HARforecast = function(vRealizedMeasure , vLags , iNRoll=10 , iNAhead=10){
  ######Initialization section ######
  iT = length(vRealizedMeasure)
  iLags = length(vLags)
  iLagsPlusOne = iLags+1
  iLagsMax = max(vLags)
  vObservations =vRealizedMeasure[iLagsMax:(iT-iNRoll-1)]  
  vForecastComp = vRealizedMeasure[(iT-iNRoll):iT]
  ######Initialization end    #######
  ######Checks section #######
  if(iNRoll > length(vRealizedMeasure)){
    stop("The amount of rolling forecasts cannot be greater than the length of the Realized measure vector.")
  }
  if(class(vRealizedMeasure)[1] == "xts" ) {
    vRealizedMeasure = as.vector(vRealizedMeasure)
  } # end xts conditional
  if(any(vRealizedMeasure<0)){
    stop("The realized measure cannot be negative. Something is wrong.")
  } # end negative conditional
  ######  Checks end   #######
  
  ######Forecasting routine #######
  
  if(iNAhead ==1 && iNRoll==1){
    #Produces only 1 forecast.
    mForecast = matrix(0 , iNAhead , ncol = iNRoll) 
    mData = HARDataCreationC(vRealizedMeasure[1:iT] , vLags) # Initialization of data to be used for forecasting
    vCoef = FASTHARestimate(vRealizedMeasure[1:iT] , vLags , iLagsPlusOne)
    mForecast[1,1] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*tail(mData,1)[2:iLagsPlusOne])
    lModel = HARestimate(vRealizedMeasure[1:(iT-iNRoll)],vLags, iLagSE = 5)
    lModel$Forecastmatrix = mForecast
    lModel$Observations = vObservations
    lModel$vForeccastComp = vForecastComp
    return(lModel)
  }
  if(iNAhead == 1){
    mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
    for (j in 1:(iNRoll)) {
      mData = HARDataCreationC(vRealizedMeasure[j:(iT-iNRoll+j-1)] , vLags) # Initialization of data to be used for forecasting.
      
      # Above makes sure the length of the model stays the same and that for each "roll" the nex period of data is used.
      vCoef = FASTHARestimate(vRealizedMeasure[j:(iT-iNRoll+j)] , vLags , iLagsPlusOne)
      #Extracts the coefficients of the model
      mForecast[1,j] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*tail(mData,1)[2:iLagsPlusOne])
      #Creates the j'th 1-step ahead forecast
    } # End loop
  }# End one-step ahead conditionals
  else{
    mForecast = matrix(1, nrow = iNAhead, ncol = iNRoll) #Initialization. Fills the matrix with ones which is used later 
    
    for (j in 1:(iNRoll)) {
      mData = HARDataCreationC(vRealizedMeasure[j:(iT-iNRoll+j-1)] , vLags)
      vCoef = FASTHARestimate(vRealizedMeasure[j:(iT-iNRoll+j)] , vLags , iLagsPlusOne)
      #Creates the j'th 1-step ahead forecast
      mForecast[1,j] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*tail(mData,1)[2:iLagsPlusOne])
      for (i in 2:(min(iLagsMax, iNAhead)+1)) {
        
        vForecastfoo = c(vRealizedMeasure[(iT - max(vLags) - (iNRoll) + (i)): (iT - (iNRoll))] , mForecast[1:i , j]) 
        
        #Gets the data necessary to form the forecast in one vector. 
        vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax-1) : length(vForecastfoo)] , vLags)
        
        #Creates the j'th i-step ahead forecast
        mForecast[i,j] = sum(vCoef*vLastrowdata)
      } #end first nested for-loop
      if(iNAhead>iLagsMax){
        for (i in (iLagsMax+1):(iNAhead)) {
          
          vForecastfoo = mForecast[(i-iLagsMax):i , j] ## Fix so the size stays the same.
          
          vLastrowdatafoo = HARDataCreationC(vForecastfoo , vLags)
          #Gets the data necessary to form the forecast in one vector.
          mForecast[i,j] = sum(vCoef*vLastrowdatafoo)
          #Creates the j'th i-step ahead forecast
        } #End second nested for-loop
      } #End conditional for the forecast periods being greater than the max of the lag-vector. 
    } #end for-loop
  }
  ######Forecasting end #######
  
  lModel = HARestimate(vRealizedMeasure[1:(iT-iNRoll)],vLags, iLagSE = 5)
  lModel$ForecastMatrix = as.data.frame(mForecast)
  names(lModel$ForecastMatrix) = paste("roll" , 1:iNRoll)
  lModel$ForecastDates = index(vForecastComp)
  lModel$Observations = vObservations
  lModel$vForecastComp = vForecastComp
  class(lModel) = "HARforecast"
  return(lModel)
}