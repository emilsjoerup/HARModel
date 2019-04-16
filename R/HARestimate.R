HARestimate = function(RM, BPV = NULL, RQ = NULL, periods = c(1,5,22), periodsJ = NULL, periodsRQ = NULL,
                       type = "HAR" , InsanityFilter = TRUE){
  ###### Initialization and preparing data ######
  iT = length(RM)
  
  ###### Initialization and data preparation end ######
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "HARQ-J", "CHAR", "CHARQ")
  if(!(type %in% vImplementedTypes)){
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  if(length(periodsRQ) > length(periods)){
    stop("periodsRQ cannot be longer than periods")
  }
  
  ##### Type: "HAR"  - vanilla
  if(type == "HAR"){
  mData = HARDataCreationC(RM, periods)
  Model = lm(mData[,1] ~ mData[,-1])
  
  names(Model$coefficients) = paste("beta", c(0,periods) , sep="")
  Info = list("Lags" = periods , "type" = type)
  vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
  if(is(RM,"xts")){
    vDates = index(RM)
    vDates = vDates[(max(periods)+1):(iT)]
    Info$Dates = vDates
  }
  
  HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                 "Data" = list("RealizedMeasure" = xts(mData[,1], vDates)) )
  
  return(HARModel)
  }
  ##### Type: "HAR"  - vanilla end
  
  
  ##### Type: "HARJ" 
  if(type == "HARJ"){
    mData = HARDataCreationC(RM, periods)
    JumpComponent = pmax(RM - BPV,0)
    if(is.null(periodsJ)){
      periodsJ = periods
    }
    
    mJumpData = HARDataCreationC(JumpComponent , periodsJ)
    
    # Dimension check
    if(dim(mData)[1] != dim(mJumpData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mJumpData)[1]
      
      if(iDimCheck<0){
        
        mJumpData = mJumpData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mJumpData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    } # Dimension check end
    
    mData = cbind(mData , mJumpData[ ,-1])
    
    Model = lm(mData[,1] ~ mData[,-1])
    if(InsanityFilter){
    Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(RM) , mean(RM) )  
    }
    
    names(Model$coefficients) = c(paste("beta", c(0,periods), sep="") , paste("j" , periodsJ, sep =""))
    
    Info = list("Lags" = periods , "JumpLags" = periodsJ , "type" = type)
    vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
    if(is(RM,"xts")){
      vDates = index(RM)
      vDates = vDates[(max(c(periods , periodsJ))+1):iT]
      Info$Dates = vDates
    }
    
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(RM[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "JumpComponent" = xts(JumpComponent[(iT-nrow(mData)+1):iT], order.by = vDates)))
    
    
    
    return(HARModel)
  }
  ##### Type: "HARJ" end

  ##### Type: "HARQ" - BPQ
  if(type == "HARQ"){
    mData = HARDataCreationC(RM, periods)
    if(is.null(periodsRQ)){
      periodsRQ = periods 
    }
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux
    
    if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
      
      if(iDimCheck<0){
        
        mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    }

    if(length(periods) != length(periodsRQ)){
      mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(periodsRQ)+1)])
    }
    else{
      mData = cbind(mData,  mAuxData[,-1] * mData[,-1])
    }
    
    Model = lm(mData[,1] ~ mData[,-1])
    if(InsanityFilter){
      Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(RM) , mean(RM) )
    }
    names(Model$coefficients) = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    
    
    Info = list("Lags" = periods , "RQLags" = periodsRQ , "type" = type)
    vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
    if(is(RM,"xts")){
      vDates = index(RM)
      vDates = vDates[(max(periods)+1):iT]
      Info$Dates = vDates
    }
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(RM[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "RealizedQuarticity" = xts(RQ[(iT-nrow(mData)+1):iT], order.by = vDates)))
    
    return(HARModel)
    
  }
  ##### Type: "HARQ" end
  
  ##### Type: "HARQ-J" - BPQ
  if(type == "HARQ-J"){
    mData = HARDataCreationC(RM, periods)
    JumpComponent = pmax(RM - BPV,0)
    if(is.null(periodsRQ)){
      periodsRQ = periods 
    }
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux
    if(is.null(periodsJ)){
      periodsJ = periods
    }
    
    #First bind RQ to RV, then Jump to the combination
    if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
      
      if(iDimCheck<0){
        
        mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    }
    if(length(periods) != length(periodsRQ)){
      mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(periodsRQ)+1)])
    }
    else{
      mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
    }
    mJumpData = HARDataCreationC(JumpComponent , periodsJ)
    # Dimension check on Jump component
    if(dim(mData)[1] != dim(mJumpData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mJumpData)[1]
      
      if(iDimCheck<0){
        
        mJumpData = mJumpData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mJumpData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    } # Dimension check end

    mData = cbind(mData , mJumpData[ ,-1])

    
    Model = lm(mData[,1] ~ mData[,-1])
    if(InsanityFilter){
      Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(RM) , mean(RM) )
    }
    names(Model$coefficients) = c(paste("beta", c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""),
                                  paste("beta_j" , periodsJ , sep = ""))
    
    
    Info = list("Lags" = periods , "RQLags" = periodsRQ, "JumpLags" = periodsJ, "type" = type)
    vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
    if(is(RM,"xts")){
      vDates = index(RM)
      vDates = vDates[(max(periods)+1):iT]
      Info$Dates = vDates
    }
    
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(RM[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "RealizedQuarticity" = xts(RQ[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "JumpComponent" = xts(JumpComponent[(iT-nrow(mData)+1):iT], order.by = vDates)))
    
    return(HARModel)
  }
  #####Type: "HARQ-J" end
  
  if(type == "CHAR"){
    mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
    Model = lm(mData[,1] ~ mData[,-1])
    
    names(Model$coefficients) = paste("beta", c(0,periods) , sep="")
    Info = list("Lags" = periods , "type" = type)
    vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
    if(is(RM,"xts")){
      vDates = index(RM)
      vDates = vDates[(max(periods)+1):(iT)]
      Info$Dates = vDates
    }
    
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(mData[,1], vDates)) )
    
    return(HARModel)
  }#####Type: "CHAR" end
  
  ##### Type: "CHARQ" - BPQ
  if(type == "CHARQ"){
    mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux

    if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
      
      if(iDimCheck<0){
        
        mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    }
    
    if(length(periods) != length(periodsRQ)){
      mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(periodsRQ)+1)])
    }
    else{
      mData = cbind(mData,  mAuxData[,-1] * mData[,-1])
    }
    
    Model = lm(mData[,1] ~ mData[,-1])
    if(InsanityFilter){
      Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(RM) , mean(RM) )
    }
    names(Model$coefficients) = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    
    
    Info = list("Lags" = periods , "RQLags" = periodsRQ , "type" = type)
    vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
    if(is(RM,"xts")){
      vDates = index(RM)
      vDates = vDates[(max(periods)+1):iT]
      Info$Dates = vDates
    }
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(RM[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "RealizedQuarticity" = xts(RQ[(iT-nrow(mData)+1):iT], order.by = vDates)))
    
    return(HARModel)
    
  }
  
  
  sError = "something unexpected happened, please report bug."
  show(sError)
  return(sError)
  
}


#######Fast estimation routine which returns only the coefficients, it is used for forecasting and cuts down the time a lot.
#######These functions are undocumented and not exported

FASTHARestimate = function(RM , periods){
  
  mData = HARDataCreationC(RM , periods)
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
}

FASTHARJestimate = function(RM, JumpComponent, periods, periodsJ ){
  mData = HARDataCreationC(RM , periods)
  
  mJumpData = HARDataCreationC(JumpComponent , periodsJ)
    # Dimension check
    if(dim(mData)[1] != dim(mJumpData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mJumpData)[1]
      if(iDimCheck<0){
        
        mJumpData = mJumpData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mJumpData will have row(s) removed
      }
      if(iDimCheck>0){
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
    } # Dimension check end
    mData = cbind(mData , mJumpData[ ,-1])
    
    vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1]) 
    #leave the last row out of the estimation and pass it back to form the forecast
    return(list("coefficients" = vCoef, "mData" = mData))
  
}

FASTHARQestimate =function(RM, RQ, periods, periodsRQ = NULL){
  if(is.null(periodsRQ)){
    periodsRQ = periods
  }
  mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
  meanAux = mean(sqrt(RQ))
  
  mAuxData = sqrt(mAuxData) - meanAux
  
  mData = HARDataCreationC(RM , periods)

  if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
    iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
    
    if(iDimCheck<0){
      
      mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
    }
    
    if(iDimCheck>0){
      
      mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
    }
    
  }
  if(length(periods) != length(periodsRQ)){
    mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(periodsRQ)+1)])
  }
  else{
    mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
  }
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


FASTHARQJestimate = function(RM , RQ , JumpComponent , periods , periodsRQ , periodsJ){
  mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
  meanAux = mean(sqrt(RQ))
  
  mAuxData = sqrt(mAuxData) - meanAux
  mData = HARDataCreationC(RM , periods)
  
  #First bind RQ to RV, then Jump to the combination
  if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
    if(iDimCheck<0){
      mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
    }
    if(iDimCheck>0){
      mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
    }
  }#dimension check end
  if(length(periods) != length(periodsRQ)){
    mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(periodsRQ)+1)])
  }
  else{
    mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
  }
  
  mJumpData = HARDataCreationC(JumpComponent , periodsJ)
  # Dimension check on Jump component
  if(dim(mData)[1] != dim(mJumpData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
    iDimCheck = dim(mData)[1] - dim(mJumpData)[1]
    if(iDimCheck<0){
      mJumpData = mJumpData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mJumpData will have row(s) removed
    }
    if(iDimCheck>0){
      mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
    }
  } # Dimension check end

  mData = cbind(mData , mJumpData[ ,-1])


  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


FASTCHARestimate = function(RM, BPV, periods, iT){
  mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
}



FASTCHARQestimate =function(RM, BPV, RQ, periods, periodsRQ, iT){
  
  mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
  
  mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
  meanAux = mean(sqrt(RQ))
  
  mAuxData = sqrt(mAuxData) - meanAux
  
  if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
    iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
    
    if(iDimCheck<0){
      
      mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
    }
    
    if(iDimCheck>0){
      
      mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
    }
    
  }
  if(length(periods) != length(periodsRQ)){
    mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(periodsRQ)+1)])
  }
  else{
    mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
  }
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


HARinsanityFilter = function(vY , dL, dU , Replacement){
  vIndex = (vY<dL | vY>dU)
  vY[vIndex] = Replacement
  return(vY)
}

