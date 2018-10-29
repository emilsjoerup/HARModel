HARestimate = function(vRealizedMeasure, vAuxData = NULL ,  vLags = c(1,5,22), vJumpLags = NULL ,show=TRUE , type = "HAR" , HARQargs = list(demean = T)){
  start.time = Sys.time()
  ######Initialization and preparing data ######
  iLags = length(vLags)
  iT = length(vRealizedMeasure)
  mData = HARDataCreationC(vRealizedMeasure, vLags)
  ######Initialization and data preparation end#
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "Full-HARQ")
  
  if(!any(grepl(type, vImplementedTypes))){
  
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  
  ##### Type: "HAR"  - vanilla
  if(type == "HAR"){
  Model = lm(mData[,1] ~ mData[,2:(iLags+1)])
  names(Model$coefficients) = paste("beta", 0:iLags , sep="")
  Info = list("Lags" = vLags , "Type" = type)
  if(is(vRealizedMeasure,"xts")){
    vDates = index(vRealizedMeasure)
    vDates = vDates[(max(vLags)+1):iT]
    Info$Dates = vDates
  }
  Info$ElapsedTime = Sys.time() - start.time
  HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                 "Data" = list("Realized Measure" = mData[,1]))
  if(show){
    show(HARModel)
  }
  return(HARModel)
  }
  ##### Type: "HAR"  - vanilla end
  
  
  ##### Type: "HARJ" 
  if(type == "HARJ"){
    if(is.null(vAuxData)){
      stop("Jump component must be provided as the vAuxData input")
    }
    if(is.null(vJumpLags)){
      warning("Jump Lags not provided, changed to match vLags")
      vJumpLags = vLags
    }
    iJumpLags = length(vJumpLags)
    mJumpData = HARDataCreationC(vAuxData , vJumpLags)
    
    # Dimension check
    if(dim(mData)[1] != dim(mJumpData)[1]){ # Checking whether the length of the realized measure and jump matrix match up
      iDimCheck = dim(mData)[1] - dim(mJumpData)[1]
      
      if(iDimCheck<0){
        #browser()
        mJumpData = mJumpData[-(1:abs(iDimCheck)),] #If iDimCheck is negative, mJumpData will have row(s) removed
      }
      
      if(iDimCheck>0){
        
        mData = mData[-(1:iDimCheck) , ] #If iDimCheck is positive, mData will have row(s) removed
      }
      
    } # Dimension check end
    
    mData = cbind(mData , mJumpData[ ,2:(iJumpLags+1)])
    
    Model = lm(mData[,1] ~ mData[,2:(iLags+1 + iJumpLags)])
    
    names(Model$coefficients) = c(paste("beta", 0:iLags, sep="") , paste("j" , 1:iJumpLags, sep =""))
    
    Info = list("Lags" = vLags , "JumpLags" = vJumpLags , "Type" = type)
    
    if(is(vRealizedMeasure,"xts")){
      vDates = index(vRealizedMeasure)
      vDates = vDates[(max(c(vLags , vJumpLags))+1):iT]
      Info$Dates = vDates
    }
    
    Info$ElapsedTime = Sys.time() - start.time
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("Realized Measure" = mData[,1],
                                 "Jump Component" = mJumpData[,1]))
    if(show){
      show(HARModel)
    }
    return(HARModel)
  }
  ##### Type: "HARJ" end
  
  ##### Type: "HARQ" - BPQ
  if(type == "HARQ"){
    
    
    
    
    vAuxDataQ = sqrt(vAuxData[max(vLags) : (length(vAuxData)-1)])
    if(HARQargs$demean){
      vAuxDataQ =  vAuxDataQ - mean(vAuxDataQ)
    }
    
    
    mData = mData
    mData = cbind(mData[,1], mData[,2:(iLags+1)] ,  vAuxDataQ*mData[,2])
    
    Model = lm(mData[,1] ~ mData[,2:dim(mData)[2]])
    names(Model$coefficients) = c(paste("beta", 0:iLags, sep="") , paste("beta_q" , 1, sep =""))
    
    
    Info = list("Lags" = vLags , "Type" = type)
    if(is(vRealizedMeasure,"xts")){
      vDates = index(vRealizedMeasure)
      vDates = vDates[(max(vLags)+1):iT]
      Info$Dates = vDates
    }
    Info$ElapsedTime = Sys.time() - start.time
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("Realized Measure" = mData[,1],
                                 "Realized Quarticity" = vAuxData[max(vLags) : (length(vAuxData)-1)]))
    if(show){
      show(HARModel)
    }
    return(HARModel)
    
  }
  ##### Type: "HARQ" end
  ##### Type: "Full-HARQ" - BPQ
  if(type == "Full-HARQ"){
    mAuxData = HARDataCreationC(vAuxData , vLags)[,-1]
    if(HARQargs$demean){
      vMeanAux = apply(mAuxData , 2, mean)
      mAuxData = mAuxData - vMeanAux
    }
    
    
    
    mData = cbind(mData[,1], mData[,2:(iLags+1)] ,  mAuxData * mData[,-1])
    
    Model = lm(mData[,1] ~ mData[,2:dim(mData)[2]])
    names(Model$coefficients) = c(paste("beta", 0:iLags, sep="") , paste("beta_q" , 1:iLags, sep =""))
    
    
    Info = list("Lags" = vLags , "Type" = type)
    if(is(vRealizedMeasure,"xts")){
      vDates = index(vRealizedMeasure)
      vDates = vDates[(max(vLags)+1):iT]
      Info$Dates = vDates
    }
    Info$ElapsedTime = Sys.time() - start.time
    HARModel = new("HARModel" , "Model" = Model ,"Info" = Info , 
                   "Data" = list("Realized Measure" = mData[,1],
                                 "Realized Quarticity" = vAuxData[max(vLags) : (length(vAuxData)-1)]))
    if(show){
      show(HARModel)
    }
    return(HARModel)
    
  }
  ##### Type: "Full-HARQ" end

  
  sError = "something unexpected happened, please report bug with code to replicate."
  return(sError)
  
}


#######Fast estimation routine which returns only the coefficients, it is used for forecasting and cuts down the time a lot.
#######These functions are undocumented and not directly user-callable######
#######iLagsPlusOne is calculated in the forecasting initialization so it can be done once, and not every roll
FASTHARestimate = function(vRealizedMeasure , vLags , iLagsPlusOne){
  mData = HARDataCreationC(vRealizedMeasure , vLags)
  lModel = lm(mData[,1] ~ mData[,2:(iLagsPlusOne)])
  return(lModel$coefficients)
}

FASTHARJestimate = function(vRealizedMeasure , vAuxData, vLags , vJumpLags , iLagsPlusOne, iJumpLags , iJumpLagsPlusOne){
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
  
    mData = cbind(mData , mJumpData[ ,2:(iJumpLagsPlusOne)])
    
    lModel = lm(mData[,1] ~ mData[,2:(iLagsPlusOne + iJumpLags)])
    
   
    return(lModel)
  
}

FASTHARQestimate =function(vRealizedMeasure , vLags , iLagsPlusOne){
  
}




