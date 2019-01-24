HARestimate = function(vRealizedMeasure, vJumpComponent = NULL, vAuxData = NULL ,  vLags = c(1,5,22), vJumpLags = NULL , vAuxLags = NULL,
                       type = "HAR" , InsanityFilter = TRUE, HARQargs = list(demean = TRUE ) ,show=TRUE ){
  start.time = Sys.time()
  ######Initialization and preparing data ######
  iT = length(vRealizedMeasure)
  mData = HARDataCreationC(vRealizedMeasure, vLags)
  ######Initialization and data preparation end#
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "HARQ-J")
  if(!any(grepl(type, vImplementedTypes))){
  
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  if(length(vAuxLags) > length(vLags)){
    stop("vAuxLags cannot be longer than vLags")
  }
  
  ##### Type: "HAR"  - vanilla
  if(type == "HAR"){
  Model = lm(mData[,1] ~ mData[,-1])
  
  names(Model$coefficients) = paste("beta", c(0,vLags) , sep="")
  Info = list("Lags" = vLags , "type" = type)
  vDates = as.Date((max(vLags)+1):iT, origin = "1970/01/01")
  if(is(vRealizedMeasure,"xts")){
    vDates = index(vRealizedMeasure)
    vDates = vDates[(max(vLags)+1):(iT)]
    Info$Dates = vDates
  }
  
  Info$"ElapsedTime" = Sys.time() - start.time
  
  HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                 "Data" = list("RealizedMeasure" = xts(mData[,1], vDates)) )
  
  if(show){
    show(HARModel)
  }
  return(HARModel)
  }
  ##### Type: "HAR"  - vanilla end
  
  
  ##### Type: "HARJ" 
  if(type == "HARJ"){
    if(is.null(vJumpComponent)){
      stop("Jump component must be provided as the vJumpComponent input")
    }
    if(is.null(vJumpLags)){
      vJumpLags = vLags
    }
    
    mJumpData = HARDataCreationC(vJumpComponent , vJumpLags)
    
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
    Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(vRealizedMeasure) , mean(vRealizedMeasure) )  
    }
    
    names(Model$coefficients) = c(paste("beta", c(0,vLags), sep="") , paste("j" , vJumpLags, sep =""))
    
    Info = list("Lags" = vLags , "JumpLags" = vJumpLags , "type" = type)
    vDates = as.Date((max(vLags)+1):iT, origin = "1970/01/01")
    if(is(vRealizedMeasure,"xts")){
      vDates = index(vRealizedMeasure)
      vDates = vDates[(max(c(vLags , vJumpLags))+1):iT]
      Info$Dates = vDates
    }
    
    Info$"ElapsedTime" = Sys.time() - start.time
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(vRealizedMeasure[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "JumpComponent" = xts(vJumpComponent[(iT-nrow(mData)+1):iT], order.by = vDates)))
    
    
    if(show){
      show(HARModel)
    }
    return(HARModel)
  }
  ##### Type: "HARJ" end

  ##### Type: "HARQ" - BPQ
  if(type == "HARQ"){
    if(is.null(vAuxLags)){
    vAuxLags = vLags 
    }
    mAuxData = as.matrix(HARDataCreationC(vAuxData , vAuxLags))
    if(HARQargs$demean){
      vMeanAux = colMeans(mAuxData)
      
      mAuxData = sqrt(mAuxData) - sqrt(vMeanAux)
    }
    
    if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
      
      if(iDimCheck<0){
        
        mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    }

    if(length(vLags) != length(vAuxLags)){
      mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(vAuxLags)+1)])
    }
    else{
      mData = cbind(mData,  mAuxData[,-1] * mData[,-1])
    }
    
    Model = lm(mData[,1] ~ mData[,-1])
    if(InsanityFilter){
      Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(vRealizedMeasure) , mean(vRealizedMeasure) )
    }
    names(Model$coefficients) = c(paste("beta",c(0,vLags), sep="") , paste("beta_q" , vAuxLags, sep =""))
    
    
    Info = list("Lags" = vLags , "RQLags" = vAuxLags , "type" = type)
    vDates = as.Date((max(vLags)+1):iT, origin = "1970/01/01")
    if(is(vRealizedMeasure,"xts")){
      vDates = index(vRealizedMeasure)
      vDates = vDates[(max(vLags)+1):iT]
      Info$Dates = vDates
    }

    Info$"ElapsedTime" = Sys.time() - start.time
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(vRealizedMeasure[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "RealizedQuarticity" = xts(vAuxData[(iT-nrow(mData)+1):iT], order.by = vDates)))
    if(show){
      show(HARModel)
    }
    return(HARModel)
    
  }
  ##### Type: "HARQ" end
  ##### Type: "HARQ-J" - BPQ
  if(type == "HARQ-J"){
    if(is.null(vAuxLags)){
      vAuxLags = vLags 
    }
    mAuxData = as.matrix(HARDataCreationC(vAuxData , vAuxLags))
    if(HARQargs$demean){
      vMeanAux = colMeans(mAuxData)
      
      mAuxData = sqrt(mAuxData) - sqrt(vMeanAux)
    }
    else{
      mAuxData = sqrt(mAuxData)
    }
    if(is.null(vJumpComponent)){
      stop("Jump component must be provided as the vJumpComponent input")
    }
    if(is.null(vJumpLags)){
      vJumpLags = vLags
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
    if(length(vLags) != length(vAuxLags)){
      mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(vAuxLags)+1)])
    }
    else{
      mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
    }
    mJumpData = HARDataCreationC(vJumpComponent , vJumpLags)
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
      Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(vRealizedMeasure) , mean(vRealizedMeasure) )
    }
    names(Model$coefficients) = c(paste("beta", c(0,vLags), sep="") , paste("beta_q" , vAuxLags, sep =""),
                                  paste("beta_j" , vJumpLags , sep = ""))
    
    
    Info = list("Lags" = vLags , "RQLags" = vAuxLags, "JumpLags" = vJumpLags, "type" = type)
    vDates = as.Date((max(vLags)+1):iT, origin = "1970/01/01")
    if(is(vRealizedMeasure,"xts")){
      vDates = index(vRealizedMeasure)
      vDates = vDates[(max(vLags)+1):iT]
      Info$Dates = vDates
    }
    
    Info$"ElapsedTime" = Sys.time() - start.time
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("RealizedMeasure" = xts(vRealizedMeasure[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "RealizedQuarticity" = xts(vAuxData[(iT-nrow(mData)+1):iT], order.by = vDates),
                                 "JumpComponent" = xts(vJumpComponent[(iT-nrow(mData)+1):iT], order.by = vDates)))
    if(show){
      show(HARModel)
    }
    return(HARModel)
  }
  #####Type: "HARQ-J" end
  
  sError = "something unexpected happened, please report bug with code to replicate."
  show(sError)
  return(sError)
  
}


#######Fast estimation routine which returns only the coefficients, it is used for forecasting and cuts down the time a lot.
#######These functions are undocumented and not directly user-callable

FASTHARestimate = function(vRealizedMeasure , vLags){
  
  mData = HARDataCreationC(vRealizedMeasure , vLags)
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
}

FASTHARJestimate = function(vRealizedMeasure, vAuxData, vLags, vJumpLags ){
  mData = HARDataCreationC(vRealizedMeasure , vLags)
  
  mJumpData = HARDataCreationC(vAuxData , vJumpLags)
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

FASTHARQestimate =function(vRealizedMeasure, vAuxData, vLags, vAuxLags = NULL, HARQargs = list(demean = T)){
  if(is.null(vAuxLags)){
    vAuxLags = vLags
  }
  mAuxData = as.matrix(HARDataCreationC(vAuxData , vAuxLags))
  # if(HARQargs$demean){
  #   vMeanAux = colMeans(mAuxData)
  #   mAuxData = sqrt(mAuxData) - sqrt(vMeanAux)
  # }
  # else{
  #   
  # }
  mAuxData = sqrt(mAuxData)

  mData = HARDataCreationC(vRealizedMeasure , vLags)

  if(dim(mData)[1] != dim(mAuxData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
    iDimCheck = dim(mData)[1] - dim(mAuxData)[1]
    
    if(iDimCheck<0){
      
      mAuxData = mAuxData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mAuxData will have row(s) removed
    }
    
    if(iDimCheck>0){
      
      mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
    }
    
  }
  if(length(vLags) != length(vAuxLags)){
    mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(vAuxLags)+1)])
  }
  else{
    mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
  }
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


FASTHARQJestimate = function(vRealizedMeasure , vAuxData , vJumpComponent , vLags , vAuxLags , vJumpLags , HARQargs = list(demean = T)){
  mAuxData = as.matrix(HARDataCreationC(vAuxData , vAuxLags))
  # if(HARQargs$demean){
  #   vMeanAux = colMeans(mAuxData)
  #   mAuxData = sqrt(mAuxData) - sqrt(vMeanAux)
  # }
  # else{
  #   
  # }
  mAuxData = sqrt(mAuxData)
  mData = HARDataCreationC(vRealizedMeasure , vLags)
  

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
  if(length(vLags) != length(vAuxLags)){
    mData = cbind(mData , mAuxData[,-1] * mData[,2:(length(vAuxLags)+1)])
  }
  else{
    mData = cbind(mData ,  mAuxData[,-1]*mData[,-1])
  }
  
  mJumpData = HARDataCreationC(vJumpComponent , vJumpLags)
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

HARinsanityFilter = function(vY , dL, dU , Replacement){
  vIndex = (vY<dL | vY>dU)
  vY[vIndex] = Replacement
  return(vY)
}

