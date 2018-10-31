HARforecast = function(vRealizedMeasure, vAuxData = NULL , vLags = c(1,5,22) , vJumpLags = NULL , iNRoll=10 , iNAhead=10 , type = "HAR" , HARQargs = list(demean = T)){
  FCstart.time = Sys.time()
  ######Initialization section ######
  iT            = length(vRealizedMeasure)
  iLags         = length(vLags)
  iLagsPlusOne  = iLags+1
  iLagsMax      = max(vLags , vJumpLags)
  vObservations = vRealizedMeasure[(iLagsMax+1):(iT-iNRoll)]
  vForecastComp = vRealizedMeasure[(iT-iNRoll+1):iT] 
  iTForecast    = length(vForecastComp)
  ######Initialization end #######
  ######Checks section #######
  if(iNRoll > (length(vRealizedMeasure)/2)){
    stop("The amount of rolling forecasts cannot be greater than half of the length of the Realized measure vector.")
  }
  if(class(vRealizedMeasure)[1] == "xts" ) {
    vRealizedMeasure = as.vector(vRealizedMeasure)
  } # end xts conditional
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ")
  
  if(!any(grepl(type, vImplementedTypes))){
    
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  ######  Checks end   #######
  
  ######Forecasting #######
  ### type: "HAR"
  if(type == "HAR"){
  if(iNAhead ==1 && iNRoll==1){
    #Produces only 1 forecast.
    mForecast      = matrix(0 , iNAhead , ncol = iNRoll) 
    lFit           = FASTHARestimate(vRealizedMeasure[-iT] , vLags , iLagsPlusOne)
    mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne)]"]]
    vCoef          = lFit$coefficients
    mForecast[1,1] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1],])
  }
  else if(iNAhead == 1){
    mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
    for (j in 1:(iNRoll)) {
    
      lFit = FASTHARestimate(vRealizedMeasure[j:(iT-iNRoll+j-1)] , vLags , iLagsPlusOne)
      mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne)]"]]
      #Extracts the coefficients of the model
      vCoef          = lFit$coefficients
      mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
      #Creates the j'th 1-step ahead forecast
    } # End loop
  }# End one-step ahead conditionals
  else{
    mForecast = matrix(NA, nrow = iNAhead, ncol = iNRoll) 
    
    for (j in 1:(iNRoll)) {
      vThisRoll = vRealizedMeasure[j:(iT-iNRoll+j-1)]

      lFit           = FASTHARestimate(vThisRoll , vLags , iLagsPlusOne)
      mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne)]"]]
      vCoef          = lFit$coefficients
      #Creates the j'th 1-step ahead forecast
      mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
      
      for (i in 2:(min((iLagsMax), iNAhead))) {
        
        vForecastfoo = c(vThisRoll[(iT-iNRoll-iLagsMax) : (iT-iNRoll)] , mForecast[1:(i-1) , j]) 
        
        
        #Gets the data necessary to form the forecast in one vector. 
        vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax) : length(vForecastfoo)] , vLags)
        
        #Creates the j'th i-step ahead forecast
        mForecast[i,j] = vCoef[1] + sum(vCoef[-1]*vLastrowdata[2:iLagsPlusOne])
      } #end first nested for-loop
      if(iNAhead>iLagsMax){
        for (i in (iLagsMax+1):(iNAhead)) {
          vForecastfoo    = mForecast[(i-iLagsMax):i , j] 
          vLastrowdatafoo = HARDataCreationC(vForecastfoo , vLags)
          
          #Gets the data necessary to form the forecast in one vector.
          mForecast[i,j] = vCoef[1] + sum(vCoef[-1]*vLastrowdatafoo[2:iLagsPlusOne])
          #Creates the j'th i-step ahead forecast
        } #End second nested for-loop
      } #End conditional for the forecast periods being greater than the max of the lag-vector. 
    } #end for-loop
  }
    FCElapsedTime = Sys.time() - FCstart.time
    
    HARForecast = new("HARForecast" , "Model" = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], vAuxData ,vLags, vJumpLags,  show=F , type = type) ,
                      "Forecast" =as.data.frame(mForecast), "Info" = list("ElapsedTime" =FCElapsedTime , "Rolls" = iNRoll ,"Horizon" = iNAhead , "type" = type),
                      "Data" = list( "ForecastDates" = index(vForecastComp), "Observations" = vObservations , "ForecastComparison" = vForecastComp ))
    
    names(HARForecast@Forecast) = paste("roll" , 1:iNRoll)
    show(HARForecast)
    return(HARForecast)
  }# type: "HAR" end
  
  # type: "HARJ" end
  if(type == "HARJ"){
    if(is.null(vJumpLags)){ # check if vJumpLags is provided
      warning("Jump Lags not provided, changed to match vLags")
      vJumpLags = vLags
    }
    iJumpLags        = length(vJumpLags)
    iJumpLagsPlusOne = iJumpLags+1
    
    if(iNAhead ==1 && iNRoll==1){
      #Produces only 1 forecast.
      mForecast      = matrix(0 , iNAhead , ncol = iNRoll) 
      lFit           = FASTHARJestimate(vRealizedMeasure[-iT], vAuxData, vLags , vJumpLags , iLagsPlusOne, iJumpLags, iJumpLagsPlusOne)
      vCoef          = lFit$coefficients
      mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne + iJumpLags)]"]]
      mForecast[1,1] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
      
    }
    else if(iNAhead == 1){
      mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
      for (j in 1:(iNRoll)) {

        lFit = FASTHARJestimate(vRealizedMeasure[j:(iT-iNRoll+j-1)], vAuxData[j:(iT-iNRoll+j-1)] , vLags , vJumpLags , iLagsPlusOne, iJumpLags , iJumpLagsPlusOne)
        #Extracts the coefficients of the model
        vCoef          = lFit$coefficients
        mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne + iJumpLags)]"]]
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      mForecast = matrix(NA, nrow = iNAhead, ncol = iNRoll) 
      vJumpForecast = rep(NA , iNAhead)
      for (j in 1:iNRoll) {
        vThisRoll     = vRealizedMeasure[j:(iT-iNRoll+j-1)]
        vThisRollJump = vAuxData[j:(iT-iNRoll+j-1)]
        lFit          = FASTHARJestimate(vThisRoll, vThisRollJump , vLags , vJumpLags , iLagsPlusOne, iJumpLags , iJumpLagsPlusOne)
        vCoef         = lFit$coefficients
        mData         = lFit[["model"]][["mData[, 2:(iLagsPlusOne + iJumpLags)]"]]
        
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
        
        vJumpForecast[1] = predict(ar.ols(vThisRollJump, FALSE ,order.max = 1))$pred
        for (i in 2:(min((iLagsMax), iNAhead))) {
          vForecastfoo     = c(vThisRoll[(iT-iNRoll-iLagsMax) : (iT-iNRoll)] , mForecast[1:(i-1) , j]) 
          vJumpForecastfoo = c(vThisRollJump[(iT-iNRoll-iLagsMax) : (iT-iNRoll)] , vJumpForecast[1:(i-1)]) 
          vJumpForecast[i] = predict(ar.ols(vThisRollJump, FALSE ,order.max = 1))$pred
          vLastrowdata     = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax) : length(vForecastfoo)] , vLags)
          vJumpLastrowdata = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) - iLagsMax) : length(vJumpForecastfoo)], vJumpLags)
          mForecast[i,j]   = vCoef[1] + sum(vCoef[-1] * c(vLastrowdata[2:iLagsPlusOne], vJumpLastrowdata[2:iJumpLagsPlusOne]))
          #Creates the j'th i-step ahead forecast
          
        } #end first nested for-loop
        if(iNAhead>iLagsMax){
          for (i in (iLagsMax+1):(iNAhead)) {
            vForecastfoo     = mForecast[(i-iLagsMax):i , j] ## Fix so the size stays the same.
            vJumpForecastfoo = vJumpForecast[(i-iLagsMax):i]
            vJumpForecast[i] = predict(ar.ols(vThisRollJump, FALSE ,order.max = 1))$pred
            vLastrowdatafoo  = HARDataCreationC(vForecastfoo , vLags)
            vJumpLastrowdata = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) - iLagsMax) : length(vJumpForecastfoo)], vJumpLags)
            mForecast[i,j]   = vCoef[1] + sum(vCoef[-1]*c(vLastrowdata[2:iLagsPlusOne], vJumpLastrowdata[2:iJumpLagsPlusOne]))
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } #End conditional for the forecast periods being greater than the max of the lag-vector. 
      } #end for-loop
    }
    FCElapsedTime = Sys.time() - FCstart.time
    
    HARForecast = new("HARForecast" , "Model" = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], vAuxData ,vLags, vJumpLags,  show=F , type = type) ,
                      "Forecast" =as.data.frame(mForecast), "Info" = list("ElapsedTime" =FCElapsedTime , "Rolls" = iNRoll ,"Horizon" = iNAhead , "type" = type),
                      "Data" = list( "ForecastDates" = index(vForecastComp), "Observations" = vObservations , "ForecastComparison" = vForecastComp ))
    
    names(HARForecast@Forecast) = paste("roll" , 1:iNRoll)
    show(HARForecast)
    return(HARForecast)
  }#type: "HARJ" end
  
  if(type == "HARQ"){
    if(iNAhead ==1 && iNRoll==1){
      #Produces only 1 forecast.
      mForecast      = matrix(0 , iNAhead , ncol = iNRoll) 
      lFit           = FASTHARQestimate(vRealizedMeasure[-iT] , vAuxData[-iT] , vLags , iLagsPlusOne ,HARQargs)
      mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne)]"]]
      vCoef          = lFit$coefficients
      mForecast[1,1] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1],])
    }
    else if(iNAhead == 1){
      mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
      for (j in 1:(iNRoll)) {
        
        lFit = FASTHARQestimate(vRealizedMeasure[j:(iT-iNRoll+j-1)], vAuxData[j:(iT-iNRoll+j-1)] , vLags , iLagsPlusOne,HARQargs)
       
        mData          = lFit[["model"]][["mData[, 2:(iLagsPlusOne)]"]]
        #Extracts the coefficients of the model
        vCoef          = lFit$coefficients
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      mForecast    = matrix(NA, nrow = iNAhead, ncol = iNRoll) 
      vAuxForecast = rep(NA, iNAhead)
      for (j in 1:(iNRoll)) {
        vThisRoll        = vRealizedMeasure[j:(iT-iNRoll+j-1)]
        vThisRollAuxData = vAuxData[j:(iT-iNRoll+j-1)]
        lFit             = FASTHARQestimate(vThisRoll ,vThisRollAuxData , vLags , iLagsPlusOne , HARQargs)
       
        mData            = lFit[["model"]][["mData[, 2:dim(mData)[2]]"]]
        vCoef            = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        vAuxForecast[1]  = predict(ar.ols(vThisRollAuxData, FALSE ,order.max = 1))$pred
        for (i in 2:(min((iLagsMax), iNAhead))) {
          
          vForecastfoo    = c(vThisRoll[(iT-iNRoll-iLagsMax) : (iT-iNRoll)] , mForecast[1:(i-1) , j]) 
          vForecastAuxFoo = c(vThisRollAuxData[(iT-iNRoll-iLagsMax) : (iT - iNRoll)] , vAuxForecast[1:(i-1)] )
          vAuxForecast[i] = predict(ar.ols(vThisRollAuxData, FALSE ,order.max = 1))$pred
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax) : length(vForecastfoo)] , vLags)

          #Creates the j'th i-step ahead forecast
          mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[2:iLagsPlusOne] , vAuxForecast[i]))
        } #end first nested for-loop
        if(iNAhead>iLagsMax){
          for (i in (iLagsMax+1):(iNAhead)) {
            vForecastfoo    = c(vThisRoll[(iT-iNRoll-iLagsMax) : (iT-iNRoll)] , mForecast[1:(i-1) , j]) 
            vForecastAuxFoo = c(vThisRollAuxData[(iT-iNRoll-iLagsMax) : (iT - iNRoll)] , vAuxForecast[1:(i-1)] )
            vAuxForecast[i] = predict(ar.ols(vThisRollAuxData, FALSE ,order.max = 1))$pred
            #Gets the data necessary to form the forecast in one vector. 
            vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax) : length(vForecastfoo)] , vLags)
            #Gets the data necessary to form the forecast in one vector.
            mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[2:iLagsPlusOne] , vAuxForecast[i]))
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } #End conditional for the forecast periods being greater than the max of the lag-vector. 
      } #end for-loop
    }
    FCElapsedTime = Sys.time() - FCstart.time
    
    HARForecast = new("HARForecast" , "Model" = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], vAuxData[1:(iT-iNRoll)] ,vLags, vJumpLags,  show=F , type = type) ,
                      "Forecast" =as.data.frame(mForecast), "Info" = list("ElapsedTime" =FCElapsedTime , "Rolls" = iNRoll ,"Horizon" = iNAhead , "type" = type),
                      "Data" = list( "ForecastDates" = index(vForecastComp), "Observations" = vObservations , "ForecastComparison" = vForecastComp ))
    
    names(HARForecast@Forecast) = paste("roll" , 1:iNRoll)
    show(HARForecast)
    return(HARForecast)
  }#type: "HARQ" end
  
  
  
  ######Forecasting end #######
  sError = "something unexpected happened, please report bug with code to replicate."
  show(sError)
  return(sError)
}