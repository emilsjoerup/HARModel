HARforecast = function(vRealizedMeasure, vJumpComponent= NULL, vAuxData = NULL , vLags = c(1,5,22) , 
                       vJumpLags = NULL, vAuxLags = NULL, iNRoll=10 , iNAhead=10 , type = "HAR",
                       InsanityFilter = TRUE ,  HARQargs = list(demean = T)){
  FCstart.time = Sys.time()
  ######Initialization section ######
  iT            = length(vRealizedMeasure)
  iLags         = length(vLags)
  iLagsPlusOne  = iLags+1
  iLagsMax      = max(vLags , vJumpLags)
  vObservations = vRealizedMeasure[(iLagsMax+1):(iT-iNRoll)]
  vForecastComp = vRealizedMeasure[(iT-iNRoll+1):iT] 
  vDates = as.Date((iLagsMax+1):(iT-iNRoll), origin = "1970/01/01")
  vForecastDates = as.Date((iT-iNRoll+1):iT, origin = "1970/01/01")
  if(is(vRealizedMeasure,"xts")){
    vDates = index(vRealizedMeasure)
    vDates = vDates[(max(vLags)+1):(iT-iNRoll)]
    vForecastDates = index(vForecastComp)
  }
  
  
  
  ######Initialization end #######
  ######Checks section #######
  if(iNRoll > (length(vRealizedMeasure) + iLagsMax)) {
    stop("The amount of rolling forecasts cannot be greater the length of the Realized measure vector 
         plus the maximum lag.")
  }
  
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "HARQ-J")
  
  if(!any(grepl(type, vImplementedTypes))){
    cat("type argument not correctly specifiec or is not implemented, available types are:", 
        paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  if(length(vAuxLags) > length(vLags)){
    stop("vAuxLags cannot be longer than vLags")
  }
  ######  Checks end   #######
  
  ######Forecasting #######
  ### type: "HAR"
  if(type == "HAR"){
  mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
  if(iNAhead == 1){
    
    for (j in 1:iNRoll) {
     lFit           = FASTHARestimate(vRealizedMeasure[j:(iT-iNRoll+j)] , vLags)
     mData          = lFit$mData[,-1]
     vCoef          = lFit$coefficients
     mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[nrow(mData),])
      #Creates the j'th 1-step ahead forecast
    } # End loop
  }# End one-step ahead conditionals
  else{
    for (j in 1:(iNRoll)){
      vThisRoll      = as.numeric(vRealizedMeasure[j:(iT-iNRoll+j)]) #pass as.numeric
      lFit           = FASTHARestimate(vThisRoll , vLags)
      mData          = lFit$mData[,-1]
      vCoef          = lFit$coefficients
      #Creates the j'th 1-step ahead forecast
      mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[nrow(mData),])
      for (i in 2:(min((iLagsMax), iNAhead))) {
        vForecastfoo = c(vThisRoll[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)] , mForecast[1:(i-1) , j]) 
        #Gets the data necessary to form the forecast in one vector. 
        vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax):length(vForecastfoo)],
                                        vLags)
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
      } 

    } #end for-loop
  }
    FCElapsedTime = Sys.time() - FCstart.time
    
    HARForecast = new("HARForecast" , "Model" = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], 
                                                            vLags = vLags, show=FALSE , type = type),
                      "Forecast" = mForecast, 
                      "Info" = list("ElapsedTime" =FCElapsedTime,
                                    "Rolls" = iNRoll ,"Horizon" = iNAhead , "type" = type),
                      "Data" = list("ForecastDates" = vForecastDates,
                                    "Observations" = xts(vObservations, vDates),
                                    "ForecastComparison" = xts(vForecastComp, vForecastDates)))
    
    colnames(HARForecast@Forecast) = paste0("roll", 1:iNRoll)
    rownames(HARForecast@Forecast) = paste0("step", 1:iNAhead)
    show(HARForecast)
    return(HARForecast)
  }# type: "HAR" end
  
  # type: "HARJ" end
  if(type == "HARJ"){
    if(is.null(vJumpComponent)){
      stop("Jump component must be provided as the vJumpComponent input")
    }
    if(is.null(vJumpLags)){ # check if vJumpLags is provided
      cat("\nJump Lags not provided, changed to match vLags")
      vJumpLags = vLags
    }
    if(max(vLags) < max(vJumpLags)){
      cat("\nHigher maximum value of jump lag vector than RV lag vector is unfortunately not 
          implemented in forecasting yet. vJumpLags is set to vLags")
      vJumpLags = vLags
    }
    iJumpLags        = length(vJumpLags)
    iJumpLagsPlusOne = iJumpLags+1
    mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
   if(iNAhead == 1){
      for (j in 1:(iNRoll)) {
        lFit = FASTHARJestimate(vRealizedMeasure[j:(iT-iNRoll+j)], vJumpComponent[j:(iT-iNRoll+j)],
                                vLags , vJumpLags)
        #Extracts the coefficients of the model
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
        if(InsanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                            max(vRealizedMeasure[j:(iT-iNRoll+j)]),
                                            mean(vRealizedMeasure[j:(iT-iNRoll+j)]))
        } # end insanityfilter conditional
        #Creates the j'th 1-step ahead forecast
        } # End loop
      }
    else{
      mJumpForecast = matrix(0, nrow = iNAhead, ncol = iNRoll)
      for (j in 1:iNRoll) {
        #Creates the j'th 1-step ahead forecast for both the realized measure and the jump component.
        vThisRoll        = as.numeric(vRealizedMeasure[j:(iT-iNRoll+j)])
        vThisRollJump    = as.numeric(vJumpComponent[j:(iT-iNRoll+j)])
        lFit             = FASTHARJestimate(vThisRoll, vThisRollJump , vLags , vJumpLags)
        mJumpFoo         = HARDataCreationC(vThisRollJump , c(1))
        dJumpY           = vJumpComponent[(iT-iNRoll+j-1)]
        lJumpFit         = fastLMcoef(cbind(1,mJumpFoo[,2]) , mJumpFoo[,1])
        dAlpha           = lJumpFit[1]
        dBeta            = lJumpFit[2]
        dMu              = dAlpha/(1-dBeta)
        mJumpForecast[1,j] = dMu + dBeta*(dJumpY - dMu)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
        
        for (i in 2:(min((iLagsMax), iNAhead))) {
          vForecastfoo = c(vThisRoll[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)] , mForecast[1:(i-1) , j]) 
          vJumpForecastfoo = c(vThisRollJump[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)] , mJumpForecast[1:(i-1),j]) 
          mJumpForecast[i,j]  = dMu + (dBeta^i)*(dJumpY - dMu)
          vLastrowdata     = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - 
                                                              iLagsMax):length(vForecastfoo)] 
                                              , vLags)
          vLastrowJumpdata = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) -
                                                                  iLagsMax):length(vJumpForecastfoo)],
                                              vJumpLags)
          mForecast[i,j]   = vCoef[1] + sum(vCoef[-1] * 
                                            c(vLastrowdata[2:iLagsPlusOne], 
                                              vLastrowJumpdata[2:iJumpLagsPlusOne]))
          #Creates the j'th i-step ahead forecast
          
        } #end first nested for-loop
        if(iNAhead>iLagsMax){
          for (i in (iLagsMax+1):(iNAhead)) {
            vForecastfoo       = mForecast[(i-iLagsMax):i , j] ## Fix so the size stays the same.
            vJumpForecastfoo   = mJumpForecast[(i-iLagsMax):i,j]
            mJumpForecast[i,j] = dMu + (dBeta^i)*(dJumpY - dMu)
            vLastrowdatafoo    = HARDataCreationC(vForecastfoo , vLags)
            vLastrowJumpdata   = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) -
                                                                    iLagsMax):length(vJumpForecastfoo)], 
                                                vJumpLags)
            mForecast[i,j]     = vCoef[1] + sum(vCoef[-1]*c(vLastrowdata[2:iLagsPlusOne], 
                                                            vLastrowJumpdata[2:iJumpLagsPlusOne]))
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(InsanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , mean(vThisRoll))
        } # end insanityfilter conditional
      } #end for-loop
    }
    FCElapsedTime = Sys.time() - FCstart.time
    
    HARForecast = new("HARForecast" , 
                      "Model" = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], 
                                            vJumpComponent = vJumpComponent,
                                            vLags = vLags, vJumpLags = vJumpLags,  
                                            show=F , type = type) ,
                      "Forecast" = mForecast, 
                      "Info" = list("ElapsedTime" =FCElapsedTime, 
                                    "Rolls" = iNRoll ,"Horizon" = iNAhead , "type" = type),
                      "Data" = list("ForecastDates" = vForecastDates,
                                    "Observations" = xts(vObservations, vDates),
                                    "ForecastComparison" = xts(vForecastComp, vForecastDates)))
    
    colnames(HARForecast@Forecast) = paste0("roll", 1:iNRoll)
    rownames(HARForecast@Forecast) = paste0("step", 1:iNAhead)
    show(HARForecast)
    return(HARForecast)
  }#type: "HARJ" end
  
  if(type == "HARQ"){
    if(is.null(vAuxData)){
      stop("Auxiliary data must be proveded as the vAuxData input")
    }
    if(is.null(vAuxLags)){
      vAuxLags = vLags 
    }
    if(max(vLags) < max(vAuxLags)){
      cat("\nHigher maximum value of auxiliary lag vector than RV lag vector is unfortunately not 
          implemented in forecasting yet. vAuxLags is set to vLags")
      vAuxLags = vLags
    }
    iMaxAuxLags      = max(vAuxLags)
    mForecast = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
    mAuxForecast = matrix(0, nrow = iNAhead, ncol = iNRoll)
    if(iNAhead == 1){
      for (j in 1:(iNRoll)) {
        
        lFit = FASTHARQestimate(vRealizedMeasure[j:(iT-iNRoll+j)], vAuxData[j:(iT-iNRoll+j)], 
                                vLags, vAuxLags, HARQargs)
       
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients

        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        if(InsanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                            max(vRealizedMeasure[j:(iT-iNRoll+j)]),
                                            mean(vRealizedMeasure[j:(iT-iNRoll+j)]))
        } # end insanityfilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(iNRoll)) {
        vThisRoll        = as.numeric(vRealizedMeasure[j:(iT-iNRoll+j)])
        vThisRollAuxData = as.numeric(vAuxData[j:(iT-iNRoll+j)])
        lFit             = FASTHARQestimate(vThisRoll, vThisRollAuxData, vLags, vAuxLags, HARQargs)
        mAuxFoo          = HARDataCreationC(vThisRollAuxData , c(1)) 
        dAuxY            = sqrt(vAuxData[(iT-iNRoll+j-1)]) #remember square root
        lAuxFit          = fastLMcoef(cbind(1,sqrt(mAuxFoo[,2])),sqrt(mAuxFoo[,1])) 
        dAlpha           = lAuxFit[1]
        dBeta            = lAuxFit[2]
        dMu              = dAlpha/(1 - dBeta)

        mAuxForecast[1,j]= dMu + dBeta * (dAuxY - dMu)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        
        for (i in 2:(min((iLagsMax), iNAhead))) {
         
          vForecastfoo      = c(vThisRoll[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)] , mForecast[1:(i-1) , j]) 
          vForecastAuxFoo   = c(sqrt(vThisRollAuxData[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)]),
                                mAuxForecast[1:(i-1) , j])
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) -
                                                               iLagsMax):length(vForecastfoo)]
                                               , vLags)
          vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                  iMaxAuxLags) : length(vForecastAuxFoo)] 
                                               , vAuxLags)
          #Creates the j'th i-step ahead forecast
          mAuxForecast[i,j] = dMu + (dBeta^i)*(dAuxY - dMu)
          mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1] , vLastrowAuxdata[-1]))
          
        } #end first nested for-loop
        if(iNAhead>iLagsMax){
          for (i in (iLagsMax+1):(iNAhead)) {
            vForecastfoo      = mForecast[(i-iLagsMax):i , j] 
            vForecastAuxFoo   =  mAuxForecast[(i-iLagsMax):i, j]
            #Gets the data necessary to form the forecast in one vector. 
            vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) -
                                                                 iLagsMax) : length(vForecastfoo)],
                                                 vLags)
            vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                    iMaxAuxLags) : length(vForecastAuxFoo)],
                                                 vAuxLags)
            mAuxForecast[i,j] = dMu + (dBeta^i)*(dAuxY - dMu)
            
            #Gets the data necessary to form the forecast in one vector.
            
            mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1] , vLastrowAuxdata[,-1]))
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(InsanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , mean(vThisRoll))
        } # end insanityfilter conditional
      } #end for-loop
    }
    FCElapsedTime = Sys.time() - FCstart.time

    HARForecast = new("HARForecast" , 
                      "Model"    = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], 
                                               vAuxData = vAuxData[1:(iT-iNRoll)],
                                               vLags = vLags, vAuxLags = vAuxLags, show=F, type = type),
                      "Forecast" = mForecast, 
                      "Info"     = list("ElapsedTime" = FCElapsedTime , "Rolls" = iNRoll,
                                        "Horizon" = iNAhead , "type" = type),
                      "Data" = list("ForecastDates" = vForecastDates,
                                    "Observations" = xts(vObservations, vDates),
                                    "ForecastComparison" = xts(vForecastComp, vForecastDates)))
  if(iNAhead>1){
    HARForecast@Data$"AuxiliaryForecast" = mAuxForecast 
  } 
    
    colnames(HARForecast@Forecast) = paste0("roll", 1:iNRoll)
    rownames(HARForecast@Forecast) = paste0("step", 1:iNAhead)
    show(HARForecast)
    return(HARForecast)
  }#type: "HARQ" end
  
  if(type == "HARQ-J"){
    if(is.null(vAuxData)){
      stop("Auxiliary data must be proveded as the vAuxData input")
    }
    if(is.null(vJumpComponent)){
      stop("Jump component must be provided as the vJumpComponent input")
    }
    if(is.null(vJumpLags)){ # check if vJumpLags is provided
      cat("\nJump Lags not provided, changed to match vLags")
      vJumpLags = vLags
    }
    if(is.null(vAuxLags)){
      vAuxLags = vLags 
    }
    if(max(vLags) < max(vJumpLags)){
      cat("\nHigher maximum value of jump lag vector than RV lag vector is 
          unfortunately not implemented in forecasting yet. vJumpLags is set to vLags")
      vJumpLags = vLags
    }
    if(max(vLags) < max(vAuxLags)){
      cat("\nHigher maximum value of auxiliary lag vector than RV lag vector is 
          unfortunately not implemented in forecasting yet. vAuxLags is set to vLags")
      vJumpLags = vLags
    }
    iMaxAuxLags      = max(vAuxLags)
    iJumpLags        = length(vJumpLags)
    iJumpLagsPlusOne = iJumpLags+1
    mJumpForecast    = matrix(0 , nrow = iNAhead, ncol = iNRoll )
    mForecast        = matrix(0 , nrow = iNAhead , ncol = iNRoll) 
    mAuxForecast     = matrix(0, nrow = iNAhead, ncol = iNRoll)
    if(iNAhead == 1){
      for (j in 1:(iNRoll)) {
        
        lFit = FASTHARQJestimate(vRealizedMeasure[j:(iT-iNRoll+j)], vAuxData[j:(iT-iNRoll+j)], 
                                 vJumpComponent[j:(iT-iNRoll+j)], vLags, 
                                 vAuxLags,vJumpLags, HARQargs)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        if(InsanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                            max(vRealizedMeasure[j:(iT-iNRoll+j)]),
                                            mean(vRealizedMeasure[j:(iT-iNRoll+j)]))
        } # end insanityfilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(iNRoll)) {
        vThisRoll          = as.numeric(vRealizedMeasure[j:(iT-iNRoll+j)])
        vThisRollAuxData   = as.numeric(vAuxData[j:(iT-iNRoll+j)])
        vThisRollJump      = as.numeric(vJumpComponent[j:(iT-iNRoll+j)])
        mJumpFoo           = HARDataCreationC(vThisRollJump , c(1))
        dJumpY             = vJumpComponent[(iT-iNRoll+j-1)]
        lJumpFit           = fastLMcoef(cbind(1,mJumpFoo[,2]) ,mJumpFoo[,1])
        dJumpAlpha         = lJumpFit[1]
        dJumpBeta          = lJumpFit[2]
        dJumpMu            = dJumpAlpha/(1-dJumpBeta)
        lFit               = FASTHARQJestimate(vThisRoll, vThisRollAuxData, vThisRollJump,
                                               vLags, vAuxLags, vJumpLags , HARQargs)
        mAuxFoo            = HARDataCreationC(vThisRollAuxData , c(1)) 
        dAuxY              = sqrt(vAuxData[(iT-iNRoll+j-1)])
        lAuxFit            = fastLMcoef(cbind(1,sqrt(mAuxFoo[,2])), sqrt(mAuxFoo[,1])) 
        dAuxAlpha          = lAuxFit[1]
        dAuxBeta           = lAuxFit[2]
        dAuxMu             = dAuxAlpha/(1 - dAuxBeta)
        mJumpForecast[1,j] = dJumpMu + dJumpBeta*(dJumpY - dJumpMu)
        mAuxForecast[1,j]  = dAuxMu + dAuxBeta * (dAuxY - dAuxMu)
        mData              = lFit$mData[,-1]
        vCoef              = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        
        for (i in 2:(min((iLagsMax), iNAhead))) {
          
          vForecastfoo      = c(vThisRoll[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)] , mForecast[1:(i-1) , j]) 
          vForecastAuxFoo   = c(sqrt(vThisRollAuxData[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)]),
                                mAuxForecast[1:(i-1) , j])
          vJumpForecastfoo  = c(vThisRollJump[(iT-iNRoll-iLagsMax+i) : (iT-iNRoll+1)],
                                mJumpForecast[1:(i-1)],j) 
          vLastrowJumpdata  = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) -
                                                                   iLagsMax) : length(vJumpForecastfoo)],
                                               vJumpLags)
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - 
                                                               iLagsMax):length(vForecastfoo)] , vLags)
          vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                  iMaxAuxLags):length(vForecastAuxFoo)], 
                                               vAuxLags)
          #Creates the j'th i-step ahead forecast
          
          mAuxForecast[i,j] = dAuxMu + (dAuxBeta^i)*(dAuxY - dAuxMu)
          mJumpForecast[i,j]= dJumpMu + (dJumpBeta^i)*(dJumpY - dJumpMu)
          mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1] , 
                                                          vLastrowAuxdata[-1] , 
                                                          vLastrowJumpdata[-1]))
          
        } #end first nested for-loop
        if(iNAhead>iLagsMax){
          for (i in (iLagsMax+1):(iNAhead)) {
            vForecastfoo      = mForecast[(i-iLagsMax):i , j] 
            vForecastAuxFoo   = mAuxForecast[(i-iLagsMax):i, j] 
            vJumpForecastfoo  = mJumpForecast[(i-iLagsMax):i, j]
            mJumpForecast[i,j]= dJumpMu + (dJumpBeta^i)*(dJumpY - dJumpMu)
            vLastrowdatafoo   = HARDataCreationC(vForecastfoo , vLags)
            vLastrowJumpdata  = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) 
                                                                   - iLagsMax) : length(vJumpForecastfoo)], 
                                                 vJumpLags)
            #Gets the data necessary to form the forecast in one vector. 
            vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) -
                                                                 iLagsMax) : length(vForecastfoo)] 
                                                 , vLags)
            vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                    iMaxAuxLags) : length(vForecastAuxFoo)] ,
                                                 vAuxLags)
            mAuxForecast[i,j] = dAuxMu + (dAuxBeta^i)*(dAuxY - dAuxMu)
            #Gets the data necessary to form the forecast in one vector.
            
            mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1],
                                                            vLastrowAuxdata[,-1],
                                                            vLastrowJumpdata[-1]))
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(InsanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , mean(vThisRoll))
        } # end insanityfilter conditional
      } #end for-loop
    }
    FCElapsedTime = Sys.time() - FCstart.time
    
    HARForecast = new("HARForecast" , 
                      "Model"    = HARestimate(vRealizedMeasure[1:(iT-iNRoll)], vAuxData = vAuxData[1:(iT-iNRoll)]
                                               , vJumpComponent = vJumpComponent[1:(iT-iNRoll)],
                                               vLags = vLags ,vJumpLags = vJumpLags, vAuxLags=vAuxLags , 
                                               show=F, type = type) ,
                      "Forecast" = mForecast, 
                      "Info"     = list("ElapsedTime" = FCElapsedTime , "Rolls" = iNRoll,
                                        "Horizon" = iNAhead , "type" = type),
                      "Data"     = list("ForecastDates" = vForecastDates,
                                        "Observations" = xts(vObservations, vDates),
                                        "ForecastComparison" = xts(vForecastComp, vForecastDates)))
    
    if(iNAhead>1){
      HARForecast@Data$"JumpFOrecast" = mJumpForecast
      HARForecast@Data$"AuxiliaryForecast" = mAuxForecast 
    } 
    
    colnames(HARForecast@Forecast) = paste0("roll", 1:iNRoll)
    rownames(HARForecast@Forecast) = paste0("step", 1:iNAhead)
    show(HARForecast)
    return(HARForecast)
  }#type: "HARQ-J" end
  
  
  ######Forecasting end #######
  sError = "something unexpected happened, please report bug."
  show(sError)
  return(sError)
}