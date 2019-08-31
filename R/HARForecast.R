HARForecast = function(RM, BPV = NULL, RQ = NULL , periods = c(1,5,22) , 
                       periodsJ = NULL, periodsRQ = NULL, nRoll=10 , nAhead=1 , type = "HAR", windowType = "rolling",
                       insanityFilter = TRUE, h = 1){
  FCstart.time = Sys.time()
  ######Initialization section ######
  iT            = length(RM)
  iLags         = length(periods)
  iLagsPlusOne  = iLags+1
  iLagsMax      = max(periods)
  vObservations = RM[(iLagsMax+1):(iT-nRoll)]
  vForecastComp = RM[(iT-nRoll+1):iT] 
  vDates = as.Date((iLagsMax+1):(iT-nRoll), origin = "1970/01/01")
  vForecastDates = as.Date((iT-nRoll+1):iT, origin = "1970/01/01")
  vRollingSynonyms = c("rolling","Rolling", "fixed", "Fixed", "RW", "rw")
  vIncreasingSynonyms = c("increasing","Increasing", "IW", "iw", "expanding", "ew", "EW")
  lFit = list()
  if(!(windowType %in% c(vRollingSynonyms, vIncreasingSynonyms))){
    warning(paste("windowType", windowType, "not in synonyms of windowtypes, thus it is set to rolling."))
    windowType = "rolling"
  }
  
  if(is(RM,"xts")){
    vDates = index(RM)
    vDates = vDates[(max(periods)+1):(iT-nRoll)]
    vForecastDates = index(vForecastComp)
  }
  if(h != 1){
    vForecastComp = HARDataCreationC(RM, periods, h = h)[-(1:(iT - nRoll - iLagsMax-h+1)),1]
    
    
  }
  
  ######Initialization end #######
  ######Checks section #######
  if(nRoll > (length(RM) + iLagsMax)) {
    stop("The amount of rolling forecasts cannot be greater the length of the Realized measure vector 
         plus the maximum lag.")
  }
  
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "HARQ-J", "CHAR", "CHARQ", "TV-HAR")
  
  if(!any(grepl(type, vImplementedTypes))){
    cat("type argument not correctly specifiec or is not implemented, available types are:", 
        paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  if(length(periodsRQ) > length(periods)){
    stop("periodsRQ cannot be longer than periods")
  }
  if(h %% 1 != 0 | h < 1){
    stop("h must be a positive integer")
  }
  if(type %in% c("HARQ", "HARQ-J", "CHARQ")){
    if(h !=1 & length(periodsRQ) == 1){
      if(!h %in% periods){
        stop("h is not in periods, that is not allowed as of now when periodsRQ is of length 1")
      }
    }
  }
  if(nAhead > 1 && h!=1){
    stop("More than 1 step ahead forecasts when h>1 have not been implemented, and currently I don't plan on doing it.\n
         Feel free to implement it yourself and send the code to me.")
  }
  ######  Checks end   #######
  
  ######Forecasting #######
  ### type: "HAR"
  if(type == "HAR"){
    mForecast = matrix(0 , nrow = nAhead , ncol = nRoll) 
    if(nAhead == 1){
      
      for (j in 1:nRoll) {
        if(windowType %in% vRollingSynonyms){
          lFit           = FASTHARestimate(RM[j:(iT-nRoll+j)], periods, h = h)
        }else{
          lFit           = FASTHARestimate(RM[1:(iT-nRoll+j)], periods, h = h)  
        }
        mData          = lFit$mData[,-1]
        vCoef          = lFit$coefficients
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[nrow(mData),])
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(nRoll)){
        if(windowType %in% vRollingSynonyms){
          vThisRoll      = as.numeric(RM[j:(iT-nRoll+j)]) #pass as.numeric
        }
        if(windowType %in% vIncreasingSynonyms){
          vThisRoll      = as.numeric(RM[1:(iT-nRoll+j)]) #pass as.numeric
        }
        lFit           = FASTHARestimate(vThisRoll , periods)
        mData          = lFit$mData[,-1]
        vCoef          = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[nrow(mData),])
        for (i in 2:(min((iLagsMax), nAhead))) {
          vForecastfoo = c(vThisRoll[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)] , mForecast[1:(i-1) , j]) 
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax):length(vForecastfoo)],
                                          periods)
          #Creates the j'th i-step ahead forecast
          mForecast[i,j] = vCoef[1] + sum(vCoef[-1]*vLastrowdata[-1])
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastfoo    = mForecast[(i-iLagsMax):i , j] 
            vLastrowdatafoo = HARDataCreationC(vForecastfoo , periods)
            #Gets the data necessary to form the forecast in one vector.
            mForecast[i,j] = vCoef[1] + sum(vCoef[-1]*vLastrowdatafoo[-1])
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        
      } #end for-loop
    }
    
  }# type: "HAR" end
  
  # type: "HARJ" end
  if(type == "HARJ"){
    JumpComponent = pmax(RM - BPV , 0)
    if(is.null(periodsJ)){ # check if periodsJ is provided
      cat("\nJump Lags not provided, changed to match periods")
      periodsJ = periods
    }
    if(max(periods) < max(periodsJ)){
      cat("\nHigher maximum value of jump lag vector than RV lag vector is unfortunately not 
          implemented in forecasting yet. periodsJ is set to periods")
      periodsJ = periods
    }
    iMaxJumpLags     = max(periodsJ)
    iJumpLags        = length(periodsJ)
    iJumpLagsPlusOne = iJumpLags+1
    mForecast = matrix(0 , nrow = nAhead , ncol = nRoll) 
    if(nAhead == 1){
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          lFit = FASTHARJestimate(RM[j:(iT-nRoll+j)], JumpComponent[j:(iT-nRoll+j)],
                                  periods , periodsJ, h = h)
        } 
        if(windowType %in% vIncreasingSynonyms){
          lFit = FASTHARJestimate(RM[1:(iT-nRoll+j)], JumpComponent[1:(iT-nRoll+j)],
                                  periods , periodsJ, h = h)
        }
        #Extracts the coefficients of the model
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
        if(insanityFilter){
          if(windowType %in% vRollingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[, 1]),
                                              mean(lFit$mData[,1]))
          }
          if(windowType %in% vIncreasingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          }
        } # end insanityFilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }
    else{
      mJumpForecast = matrix(0, nrow = nAhead, ncol = nRoll)
      for (j in 1:nRoll) {
        #Creates the j'th 1-step ahead forecast for both the realized measure and the jump component.
        if(windowType %in% vRollingSynonyms){
          vThisRoll        = as.numeric(RM[j:(iT-nRoll+j)]) #pass as.numeric
          vThisRollJump    = as.numeric(JumpComponent[j:(iT-nRoll+j)])
        }
        if(windowType %in% vIncreasingSynonyms){
          vThisRoll        = as.numeric(RM[1:(iT-nRoll+j)])
          vThisRollJump    = as.numeric(JumpComponent[1:(iT-nRoll+j)])
        }
        lFit             = FASTHARJestimate(vThisRoll, vThisRollJump , periods , periodsJ)
        mJumpFoo         = HARDataCreationC(vThisRollJump , c(1))
        dJumpY           = JumpComponent[(iT-nRoll+j-1)]
        lJumpFit         = fastLMcoef(cbind(1,mJumpFoo[,2]) , mJumpFoo[,1])
        dAlpha           = lJumpFit[1]
        dBeta            = lJumpFit[2]
        dMu              = dAlpha/(1-dBeta)
        mJumpForecast[1,j] = dMu + dBeta*(dJumpY - dMu)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1], ])
        
        for (i in 2:(min((iLagsMax), nAhead))) {
          vForecastfoo = c(vThisRoll[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)] , mForecast[1:(i-1) , j]) 
          vJumpForecastfoo = c(vThisRollJump[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)] , mJumpForecast[1:(i-1),j]) 
          mJumpForecast[i,j]  = dMu + (dBeta^i)*(dJumpY - dMu)
          vLastrowdata     = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - 
                                                              iLagsMax):length(vForecastfoo)] 
                                              , periods)
          vLastrowJumpdata = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) -
                                                                  iMaxJumpLags):length(vJumpForecastfoo)],
                                              periodsJ)
          mForecast[i,j]   = vCoef[1] + sum(vCoef[-1] * 
                                              c(vLastrowdata[-1], 
                                                vLastrowJumpdata[-1]))
          #Creates the j'th i-step ahead forecast
          
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastfoo       = mForecast[(i-iLagsMax):i , j] ## Fix so the size stays the same.
            vJumpForecastfoo   = mJumpForecast[(i-iLagsMax):i,j]
            mJumpForecast[i,j] = dMu + (dBeta^i)*(dJumpY - dMu)
            vLastrowdatafoo    = HARDataCreationC(vForecastfoo , periods)
            vLastrowJumpdata   = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) -
                                                                      iMaxJumpLags):length(vJumpForecastfoo)], 
                                                  periodsJ)
            mForecast[i,j]     = vCoef[1] + sum(vCoef[-1]*c(vLastrowdata[-1], 
                                                            vLastrowJumpdata[-1]))
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(insanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , sum(vThisRoll)/length(vThisRoll))
        } # end insanityFilter conditional
      } #end for-loop
    }
    
  }#type: "HARJ" end
  
  if(type == "HARQ"){
    if(is.null(RQ)){
      stop("Auxiliary data must be proveded as the RQ input")
    }
    if(is.null(periodsRQ)){
      periodsRQ = periods 
    }
    if(max(periods) < max(periodsRQ)){
      cat("\nHigher maximum value of auxiliary lag vector than RV lag vector is unfortunately not 
          implemented in forecasting yet. periodsRQ is set to periods")
      periodsRQ = periods
    }
    iMaxAuxLags      = max(periodsRQ)
    mForecast = matrix(0 , nrow = nAhead , ncol = nRoll) 
    mAuxForecast = matrix(0, nrow = nAhead, ncol = nRoll)
    if(nAhead == 1){
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          lFit = FASTHARQestimate(RM[j:(iT-nRoll+j)], RQ[j:(iT-nRoll+j)], 
                                  periods, periodsRQ, h = h)
        } 
        if(windowType %in% vIncreasingSynonyms){
          lFit = FASTHARQestimate(RM[1:(iT-nRoll+j)], RQ[1:(iT-nRoll+j)], 
                                  periods, periodsRQ, h = h)
         
        }
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        if(insanityFilter){
          if(windowType %in% vRollingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          } 
          if(windowType %in% vIncreasingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          }
        } # end insanityFilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          vThisRoll        = as.numeric(RM[j:(iT-nRoll+j)])
          vThisRollRQ      = as.numeric(RQ[j:(iT-nRoll+j)])
        } 
        if(windowType %in% vIncreasingSynonyms){
          vThisRoll        = as.numeric(RM[1:(iT-nRoll+j)])
          vThisRollRQ      = as.numeric(RQ[1:(iT-nRoll+j)])
        }
        lFit             = FASTHARQestimate(vThisRoll, vThisRollRQ, periods, periodsRQ)
        mAuxFoo          = HARDataCreationC(vThisRollRQ , c(1)) 
        dAuxY            = sqrt(RQ[(iT-nRoll+j-1)]) #remember square root
        lAuxFit          = fastLMcoef(cbind(1,sqrt(mAuxFoo[,2])),sqrt(mAuxFoo[,1])) 
        dAlpha           = lAuxFit[1]
        dBeta            = lAuxFit[2]
        dMu              = dAlpha/(1 - dBeta)
        
        mAuxForecast[1,j]= dMu + dBeta * (dAuxY - dMu)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        
        for (i in 2:(min((iLagsMax), nAhead))) {
          
          vForecastfoo      = c(vThisRoll[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)] , mForecast[1:(i-1) , j]) 
          vForecastAuxFoo   = c(sqrt(vThisRollRQ[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)]),
                                mAuxForecast[1:(i-1) , j])
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) -
                                                               iLagsMax):length(vForecastfoo)]
                                               , periods)
          vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                  iMaxAuxLags) : length(vForecastAuxFoo)] 
                                               , periodsRQ)
          #Creates the j'th i-step ahead forecast
          mAuxForecast[i,j] = dMu + (dBeta^i)*(dAuxY - dMu)
          mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1] , vLastrowAuxdata[-1]))
          
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastfoo      = mForecast[(i-iLagsMax):i , j] 
            vForecastAuxFoo   =  mAuxForecast[(i-iLagsMax):i, j]
            #Gets the data necessary to form the forecast in one vector. 
            vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) -
                                                                 iLagsMax) : length(vForecastfoo)],
                                                 periods)
            vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                    iMaxAuxLags) : length(vForecastAuxFoo)],
                                                 periodsRQ)
            mAuxForecast[i,j] = dMu + (dBeta^i)*(dAuxY - dMu)
            
            #Gets the data necessary to form the forecast in one vector.
            
            mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1] , vLastrowAuxdata[,-1]))
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(insanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , sum(vThisRoll)/length(vThisRoll))
        } # end insanityFilter conditional
      } #end for-loop
    }
    
  }#type: "HARQ" end
  
  if(type == "HARQ-J"){
    JumpComponent = pmax(RM - BPV , 0)
    if(is.null(RQ)){
      stop("Auxiliary data must be proveded as the RQ input")
    }
    if(is.null(periodsJ)){ # check if periodsJ is provided
      cat("\nJump Lags not provided, changed to match periods")
      periodsJ = periods
    }
    if(is.null(periodsRQ)){
      periodsRQ = periods 
    }
    if(max(periods) < max(periodsJ)){
      cat("\nHigher maximum value of jump lag vector than RV lag vector is 
          unfortunately not implemented in forecasting yet. periodsJ is set to periods")
      periodsJ = periods
    }
    if(max(periods) < max(periodsRQ)){
      cat("\nHigher maximum value of auxiliary lag vector than RV lag vector is 
          unfortunately not implemented in forecasting yet. periodsRQ is set to periods")
      periodsJ = periods
    }
    iMaxAuxLags      = max(periodsRQ)
    iMaxJumpLags     = max(periodsJ)
    mJumpForecast    = matrix(0 , nrow = nAhead, ncol = nRoll )
    mForecast        = matrix(0 , nrow = nAhead , ncol = nRoll) 
    mAuxForecast     = matrix(0, nrow = nAhead, ncol = nRoll)
    if(nAhead == 1){
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          lFit = FASTHARQJestimate(RM[j:(iT-nRoll+j)], RQ[j:(iT-nRoll+j)], 
                                   JumpComponent[j:(iT-nRoll+j)], periods, 
                                   periodsRQ,periodsJ, h = h)
        } 
        if(windowType %in% vIncreasingSynonyms){
          lFit = FASTHARQJestimate(RM[1:(iT-nRoll+j)], RQ[1:(iT-nRoll+j)], 
                                   JumpComponent[1:(iT-nRoll+j)], periods, 
                                   periodsRQ,periodsJ, h = h)
        }
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        if(insanityFilter){
          if(windowType %in% vRollingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          } 
          if(windowType %in% vIncreasingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          }
        } # end insanityFilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          vThisRoll          = as.numeric(RM[j:(iT-nRoll+j)])
          vThisRollRQ        = as.numeric(RQ[j:(iT-nRoll+j)])
          vThisRollJump      = as.numeric(JumpComponent[j:(iT-nRoll+j)])
        } else{
          vThisRoll          = as.numeric(RM[1:(iT-nRoll+j)])
          vThisRollRQ        = as.numeric(RQ[1:(iT-nRoll+j)])
          vThisRollJump      = as.numeric(JumpComponent[1:(iT-nRoll+j)])
        }
        mJumpFoo           = HARDataCreationC(vThisRollJump , c(1))
        dJumpY             = JumpComponent[(iT-nRoll+j-1)]
        lJumpFit           = fastLMcoef(cbind(1,mJumpFoo[,2]) ,mJumpFoo[,1])
        dJumpAlpha         = lJumpFit[1]
        dJumpBeta          = lJumpFit[2]
        dJumpMu            = dJumpAlpha/(1-dJumpBeta)
        lFit               = FASTHARQJestimate(vThisRoll, vThisRollRQ, vThisRollJump,
                                               periods, periodsRQ, periodsJ)
        mAuxFoo            = HARDataCreationC(vThisRollRQ , c(1)) 
        dAuxY              = sqrt(RQ[(iT-nRoll+j-1)])
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
        
        for (i in 2:(min((iLagsMax), nAhead))) {
          
          vForecastfoo      = c(vThisRoll[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)] , mForecast[1:(i-1) , j]) 
          vForecastAuxFoo   = c(sqrt(vThisRollRQ[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)]),
                                mAuxForecast[1:(i-1) , j])
          vJumpForecastfoo  = c(vThisRollJump[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)],
                                mJumpForecast[1:(i-1)],j) 
          vLastrowJumpdata  = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) -
                                                                   iMaxJumpLags) : length(vJumpForecastfoo)],
                                               periodsJ)
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - 
                                                               iLagsMax):length(vForecastfoo)] , periods)
          vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                  iMaxAuxLags):length(vForecastAuxFoo)], 
                                               periodsRQ)
          #Creates the j'th i-step ahead forecast
          
          mAuxForecast[i,j] = dAuxMu + (dAuxBeta^i)*(dAuxY - dAuxMu)
          mJumpForecast[i,j]= dJumpMu + (dJumpBeta^i)*(dJumpY - dJumpMu)
          mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1] , 
                                                          vLastrowAuxdata[-1] , 
                                                          vLastrowJumpdata[-1]))
          
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastfoo      = mForecast[(i-iLagsMax):i , j] 
            vForecastAuxFoo   = mAuxForecast[(i-iLagsMax):i, j] 
            vJumpForecastfoo  = mJumpForecast[(i-iLagsMax):i, j]
            vLastrowJumpdata  = HARDataCreationC(vJumpForecastfoo[(length(vJumpForecastfoo) 
                                                                   - iMaxJumpLags) : length(vJumpForecastfoo)], 
                                                 periodsJ)
            #data preparation
            vLastrowdata      = HARDataCreationC(vForecastfoo[(length(vForecastfoo) -
                                                                 iLagsMax) : length(vForecastfoo)] 
                                                 , periods)
            vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                    iMaxAuxLags) : length(vForecastAuxFoo)] ,
                                                 periodsRQ)
            
            
            #Creates the j'th i-step ahead forecast
            mJumpForecast[i,j]= dJumpMu + (dJumpBeta^i)*(dJumpY - dJumpMu)
            mAuxForecast[i,j] = dAuxMu + (dAuxBeta^i)*(dAuxY - dAuxMu)
            mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[-1],
                                                            vLastrowAuxdata[,-1],
                                                            vLastrowJumpdata[-1]))
            
          } #End second nested for-loop
        } 
        if(insanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , sum(vThisRoll)/length(vThisRoll))
        } # end insanityFilter conditional
      } #end for-loop
    }
    
    
    
  }#type: "HARQ-J" end
  
  if(type == "CHAR"){
    mForecast = matrix(0 , nrow = nAhead , ncol = nRoll) 
    mBPVForecast = matrix(0, nrow = nAhead, ncol = nRoll)
    iN = iT-nRoll+1 #Used for fast CHAR estimation, need better name for it tho
    if(nAhead == 1){
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          lFit = FASTCHARestimate(RM[j:(iT-nRoll+j)], BPV[j:(iT-nRoll+j)], 
                                  periods, iN, h = h)
        }
        if(windowType %in% vIncreasingSynonyms){
          lFit = FASTCHARestimate(RM[1:(iT-nRoll+j)], BPV[1:(iT-nRoll+j)], 
                                  periods, iN + j -1, h = h)
        }
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        if(insanityFilter){
          if(windowType %in% vRollingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          } 
          if(windowType %in% vIncreasingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          }
        } # end insanityFilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          vThisRoll        = as.numeric(RM[j:(iT-nRoll+j)])
          vThisRollBPV     = as.numeric(BPV[j:(iT-nRoll+j)])
        } 
        if(windowType %in% vIncreasingSynonyms){
          vThisRoll        = as.numeric(RM[1:(iT-nRoll+j)])
          vThisRollBPV     = as.numeric(BPV[1:(iT-nRoll+j)])
          iN = iT-nRoll+j
        }
        lFit             = FASTCHARestimate(vThisRoll, vThisRollBPV, periods, iN)
        mBPVFoo          = HARDataCreationC(vThisRollBPV , c(1)) 
        dBPVY            = BPV[(iT-nRoll+j-1)]
        lBPVFit          = fastLMcoef(cbind(1,mBPVFoo[,2]),mBPVFoo[,1]) 
        dAlpha           = lBPVFit[1]
        dBeta            = lBPVFit[2]
        dMu              = dAlpha/(1 - dBeta)
        mBPVForecast[1,j]= dMu + dBeta * (dBPVY - dMu)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        
        for (i in 2:(min((iLagsMax), nAhead))) {
          
          vForecastBPVFoo   = c(vThisRollBPV[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)],
                                mBPVForecast[1:(i-1) , j])
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowBPVdata   = HARDataCreationC(vForecastBPVFoo[(length(vForecastBPVFoo) -
                                                                  iLagsMax) : length(vForecastBPVFoo)] 
                                               , periods)
          #Creates the j'th i-step ahead forecast
          mBPVForecast[i,j] = dMu + (dBeta^i)*(dBPVY - dMu)
          mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowBPVdata[-1]))
          
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastBPVFoo   = mBPVForecast[(i-iLagsMax):i, j]
            #Gets the data necessary to form the forecast in one vector. 
            vLastrowBPVdata   = HARDataCreationC(vForecastBPVFoo[(length(vForecastBPVFoo) -
                                                                    iLagsMax) : length(vForecastBPVFoo)],
                                                 periods)
            mBPVForecast[i,j] = dMu + (dBeta^i)*(dBPVY - dMu)
            
            #Gets the data necessary to form the forecast in one vector.
            
            mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* vLastrowBPVdata[,-1])
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(insanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , mean(vThisRoll))
        } # end insanityFilter conditional
      } #end for-loop
    }
    
    
  }#type: "CHAR" end
  
  
  if(type == "CHARQ"){
    mForecast = matrix(0 , nrow = nAhead , ncol = nRoll) 
    mBPVForecast = matrix(0, nrow = nAhead, ncol = nRoll)
    mAuxForecast     = matrix(0, nrow = nAhead, ncol = nRoll)
    iMaxAuxLags      = max(periodsRQ)
    iN = iT-nRoll+1 #Used for fast CHAR estimation, need better name for it tho
    if(nAhead == 1){
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          lFit = FASTCHARQestimate(RM[j:(iT-nRoll+j)], BPV[j:(iT-nRoll+j)], RQ[j:(iT-nRoll+j)],
                                   periods, periodsRQ, iN, h = h)
        } 
        if(windowType %in% vIncreasingSynonyms){
          lFit = FASTCHARQestimate(RM[1:(iT-nRoll+j)], BPV[1:(iT-nRoll+j)], RQ[1:(iT-nRoll+j)], 
                                   periods, periodsRQ, iN + j -1, h = h)
        }
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        if(insanityFilter){
          if(windowType %in% vRollingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          } else {
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          }
        } # end insanityFilter conditional
        #Creates the j'th 1-step ahead forecast
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(nRoll)) {
        if(windowType %in% vRollingSynonyms){
          vThisRollRQ      = as.numeric(RQ[j:(iT-nRoll+j)])
          vThisRoll        = as.numeric(RM[j:(iT-nRoll+j)])
          vThisRollBPV     = as.numeric(BPV[j:(iT-nRoll+j)])
        } 
        if(windowType %in% vIncreasingSynonyms){
          vThisRollRQ      = as.numeric(RQ[1:(iT-nRoll+j)])
          vThisRoll        = as.numeric(RM[1:(iT-nRoll+j)])
          vThisRollBPV     = as.numeric(BPV[1:(iT-nRoll+j)])
          iN = iT-nRoll+j
        }
        lFit             = FASTCHARQestimate(vThisRoll, vThisRollBPV, vThisRollRQ, periods, periodsRQ,iN)
        ##BPV modeling
        mBPVFoo          = HARDataCreationC(vThisRollRQ , c(1)) 
        dBPVY            = BPV[(iT-nRoll+j-1)]
        lBPVFit          = fastLMcoef(cbind(1,mBPVFoo[,2]),mBPVFoo[,1]) 
        dAlpha           = lBPVFit[1]
        dBeta            = lBPVFit[2]
        dMu              = dAlpha/(1 - dBeta)
        mBPVForecast[1,j]= dMu + dBeta * (dBPVY - dMu)
        ##RQ modeling
        mAuxFoo          = HARDataCreationC(vThisRollRQ , c(1)) 
        dAuxY            = sqrt(RQ[(iT-nRoll+j-1)]) #remember square root
        lAuxFit          = fastLMcoef(cbind(1,sqrt(mAuxFoo[,2])),sqrt(mAuxFoo[,1])) 
        dAlphaRQ         = lAuxFit[1]
        dBetaRQ          = lAuxFit[2]
        dMuRQ            = dAlphaRQ/(1 - dBetaRQ)
        mAuxForecast[1,j]= dMuRQ + dBetaRQ * (dAuxY - dMuRQ)
        mData            = lFit$mData[,-1]
        vCoef            = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j]   = vCoef[1] + sum(vCoef[-1]*mData[dim(mData)[1] ,])
        
        for (i in 2:(min((iLagsMax), nAhead))) {
          
          
          vForecastBPVFoo   = c(vThisRollRQ[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)],
                                mBPVForecast[1:(i-1) , j])
          vForecastAuxFoo   = c(sqrt(vThisRollRQ[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)]),
                                mAuxForecast[1:(i-1) , j])
          #data preparation
          vLastrowBPVdata   = HARDataCreationC(vForecastBPVFoo[(length(vForecastBPVFoo) -
                                                                  iLagsMax) : length(vForecastBPVFoo)] 
                                               , periods)
          vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                  iMaxAuxLags):length(vForecastAuxFoo)], 
                                               periodsRQ)
          #Creates the j'th i-step ahead forecasts
          mAuxForecast[i,j] = dMuRQ + (dBetaRQ^i)*(dAuxY - dMuRQ)
          mBPVForecast[i,j] = dMu + (dBeta^i)*(dBPVY - dMu)
          mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowBPVdata[-1], vLastrowAuxdata[,-1]))
          
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastAuxFoo   = mAuxForecast[(i-iLagsMax):i, j] 
            vForecastBPVFoo   =  mBPVForecast[(i-iLagsMax):i, j]
            #Gets the data necessary to form the forecast in one vector. 
            vLastrowBPVdata   = HARDataCreationC(vForecastBPVFoo[(length(vForecastBPVFoo) -
                                                                    iLagsMax) : length(vForecastBPVFoo)],
                                                 periods)
            vLastrowAuxdata   = HARDataCreationC(vForecastAuxFoo[(length(vForecastAuxFoo) -
                                                                    iMaxAuxLags) : length(vForecastAuxFoo)] ,
                                                 periodsRQ)
            mBPVForecast[i,j] = dMu + (dBeta^i)*(dBPVY - dMu)
            mAuxForecast[i,j] = dMuRQ + (dBetaRQ^i)*(dAuxY - dMuRQ)
            mForecast[i,j]    = vCoef[1] + sum(vCoef[-1]* c(vLastrowBPVdata[-1], vLastrowAuxdata[,-1]))
            
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(insanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , mean(vThisRoll))
        } # end insanityFilter conditional
      } #end for-loop
    }
    
    
    
  }#type: "CHARQ" end
  
  ### type: "TV-HAR"
  if(type == "TV-HAR"){
    mForecast = matrix(0 , nrow = nAhead , ncol = nRoll) 
    if(nAhead == 1){
      
      for (j in 1:nRoll) {
        if(windowType %in% vRollingSynonyms){
          lFit           = FASTTVHARestimate(RM[j:(iT-nRoll+j)] , periods, h = h)
        }else{
          lFit           = FASTTVHARestimate(RM[1:(iT-nRoll+j)] , periods, h = h)  
        }
        mData          = lFit$mData[,-1]
        vCoef          = lFit$coefficients
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1]*mData[nrow(mData),])
        #Creates the j'th 1-step ahead forecast
        if(insanityFilter){
          if(windowType %in% vRollingSynonyms){
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          } else {
            mForecast[,j] = HARinsanityFilter(mForecast[,j], 0,
                                              max(lFit$mData[,1]),
                                              mean(lFit$mData[,1]))
          }
        } # end insanityFilter conditional
      } # End loop
    }# End one-step ahead conditionals
    else{
      for (j in 1:(nRoll)){
        if(windowType %in% vRollingSynonyms){
          vThisRoll      = as.numeric(RM[j:(iT-nRoll+j)]) #pass as.numeric
        }
        if(windowType %in% vIncreasingSynonyms){
          vThisRoll      = as.numeric(RM[1:(iT-nRoll+j)]) #pass as.numeric
        }
        lFit           = FASTTVHARestimate(vThisRoll , periods)
        mData          = lFit$mData[,-1]
        vCoef          = lFit$coefficients
        #Creates the j'th 1-step ahead forecast
        mForecast[1,j] = vCoef[1] + sum(vCoef[-1] * mData[nrow(mData),])
        
        
        for (i in 2:(min((iLagsMax), nAhead))) {
          vForecastfoo = c(vThisRoll[(iT-nRoll-iLagsMax+i) : (iT-nRoll+1)] , mForecast[1:(i-1) , j]) 
          #Gets the data necessary to form the forecast in one vector. 
          vLastrowdata = HARDataCreationC(vForecastfoo[(length(vForecastfoo) - iLagsMax):length(vForecastfoo)],
                                          periods)
          #Creates the j'th i-step ahead forecast
          mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[2], 
                                                       abs(vLastrowdata[2]-vLastrowdata[iLagsPlusOne]), #difference
                                                       vLastrowdata[3:iLagsPlusOne]))
          
        } #end first nested for-loop
        if(nAhead>iLagsMax){
          for (i in (iLagsMax+1):(nAhead)) {
            vForecastfoo    = mForecast[(i-iLagsMax):i , j] 
            vLastrowdatafoo = HARDataCreationC(vForecastfoo , periods)
            #Gets the data necessary to form the forecast in one vector.
            mForecast[i,j] = vCoef[1] + sum(vCoef[-1]* c(vLastrowdata[2], 
                                                         abs(vLastrowdata[2]-vLastrowdata[iLagsPlusOne]), #difference
                                                         vLastrowdata[3:iLagsPlusOne]))
            #Creates the j'th i-step ahead forecast
          } #End second nested for-loop
        } 
        if(insanityFilter){
          mForecast[,j] = HARinsanityFilter(mForecast[,j] , 0 , max(vThisRoll) , sum(vThisRoll)/length(vThisRoll))
        }
      } #end for-loop
    }
    
  }# type: "TV-HAR" end
  FCElapsedTime = Sys.time() - FCstart.time
  HARForecast = new("HARForecast" , 
                    "model"    = HAREstimate(RM[1:(iT-nRoll)], 
                                             BPV = BPV[1:(iT-nRoll)],
                                             RQ  = RQ[1:(iT-nRoll)],
                                             periods = periods, periodsRQ = periodsRQ, type = type, h = h),
                    "forecast" = mForecast, 
                    "info"     = list("elapsedTime" = FCElapsedTime , "rolls" = nRoll,
                                      "horizon" = nAhead , "type" = type, "windowType" = windowType),
                    "data" = list("forecastDates" = vForecastDates,
                                  "observations" = xts(vObservations, vDates),
                                  "forecastComparison" = xts(vForecastComp, vForecastDates)))
  colnames(HARForecast@forecast) = paste0("roll", 1:nRoll)
  rownames(HARForecast@forecast) = paste0("step", 1:nAhead)
  
  return(HARForecast)
  
}