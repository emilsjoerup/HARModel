########################################################################################################
# This package is created by Emil Sjoerup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################

HARforecast = function(vRealizedMeasure , vLags , iNRoll=10 , iNAhead=10){
  FCstart.time = Sys.time()
  ######Initialization section ######
  iT = length(vRealizedMeasure)
  iLags = length(vLags)
  iLagsPlusOne = iLags+1
  iLagsMax = max(vLags)
  #vAllDates = index(vRealizedMeasure) # Will be used in some sort of fix of the
  vObservations = vRealizedMeasure[(iLagsMax+1):(iT-iNRoll)]
  vForecastComp = vRealizedMeasure[(iT-iNRoll+1):iT] #Up to here should be fine
  iTForecast = length(vForecastComp)
  #browser()
  length(c(vObservations , vForecastComp))
  length(vRealizedMeasure)
  # #####Initialization end #######
  ######Checks section #######
  if(iNRoll > length(vRealizedMeasure)){
    stop("The amount of rolling forecasts cannot be greater than the length of the Realized measure vector.")
  }
  if(class(vRealizedMeasure)[1] == "xts" ) {
    vRealizedMeasure = as.vector(vRealizedMeasure)
  } # end xts conditional
  ######  Checks end   #######
  
  ######Forecasting routine #######
  
  if(iNAhead ==1 && iNRoll==1){
    #Produces only 1 forecast.
    mForecast = matrix(0 , iNAhead , ncol = iNRoll) 
    mData = HARDataCreationC(vRealizedMeasure[1:iT] , vLags) # Initialization of data to be used for forecasting
    vCoef = FASTHARestimate(vRealizedMeasure[1:iT] , vLags , iLagsPlusOne)
    mForecast[1,1] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*tail(mData,1)[2:iLagsPlusOne])
  }
  else if(iNAhead == 1){
    mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
    for (j in 1:(iNRoll)) {
      mData = HARDataCreationC(vRealizedMeasure[j:(iT-iNRoll+j-1)] , vLags) # Initialization of data to be used for forecasting.
      
      # Above makes sure the length of the model stays the same and that for each "roll" the nex period of data is used.
      vCoef = FASTHARestimate(vRealizedMeasure[j:(iT-iNRoll+j-1)] , vLags , iLagsPlusOne)
      #Extracts the coefficients of the model
      mForecast[1,j] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*tail(mData,1)[2:iLagsPlusOne])
      #Creates the j'th 1-step ahead forecast
    } # End loop
  }# End one-step ahead conditionals
  else{
    mForecast = matrix(NA, nrow = iNAhead, ncol = iNRoll) 
    
    for (j in 1:(iNRoll)) {
      vThisRoll = vRealizedMeasure[(iT-iNRoll+j-1 - iNRoll):(iT-iNRoll+j-1)]
      
      mData = HARDataCreationC(vThisRoll , vLags)
      #print(dim(mData))
      vCoef = FASTHARestimate(vThisRoll , vLags , iLagsPlusOne)
      #Creates the j'th 1-step ahead forecast
      mForecast[1,j] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*tail(mData,1)[2:iLagsPlusOne])
      for (i in 2:(min((iLagsMax), iNAhead))) {
        
        vForecastfoo = c(vThisRoll[i:iTForecast] , mForecast[1:(i-1) , j]) #TODO: FIX
        #browser()
        
        #Gets the data necessary to form the forecast in one vector. 
        vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax) : length(vForecastfoo)] , vLags)
        
        #Creates the j'th i-step ahead forecast
        mForecast[i,j] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*vLastrowdata[2:iLagsPlusOne])
      } #end first nested for-loop
      if(iNAhead>iLagsMax){
        for (i in (iLagsMax+1):(iNAhead)) {
          
          vForecastfoo = mForecast[(i-iLagsMax):i , j] ## Fix so the size stays the same.
          
          vLastrowdatafoo = HARDataCreationC(vForecastfoo , vLags)
          
          #Gets the data necessary to form the forecast in one vector.
          mForecast[i,j] = vCoef[1] + sum(vCoef[2:iLagsPlusOne]*vLastrowdatafoo[2:iLagsPlusOne])
          #Creates the j'th i-step ahead forecast
        } #End second nested for-loop
      } #End conditional for the forecast periods being greater than the max of the lag-vector. 
    } #end for-loop
  }
  ######Forecasting end #######
  FCElapsedTime = Sys.time() - FCstart.time
  
  HARForecast = new("HARForecast" , "Model" = HARestimate(vRealizedMeasure[1:(iT-iNRoll)],vLags,  show=F) , "Forecast" =as.data.frame(mForecast), "Info" = list("ElapsedTime" =FCElapsedTime , "Rolls" = iNRoll ,"Horizon" = iNAhead),"Data" = list( "ForecastDates" = index(vForecastComp), "Observations" = vObservations , "ForecastComparison" = vForecastComp ) )
  names(HARForecast@Forecast) = paste("roll" , 1:iNRoll)
  show(HARForecast)
  return(HARForecast)
}