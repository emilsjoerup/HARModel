HARestimate = function(RM, BPV = NULL, RQ = NULL, periods = c(1,5,22), periodsJ = NULL, periodsRQ = NULL,
                       type = "HAR" , InsanityFilter = TRUE){
  ###### Initialization and preparing data ######
  iT = length(RM)
  Info = list()
  ###### Initialization and data preparation end ######
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "HARQ-J", "CHAR", "CHARQ", "TV-HAR")
  if(!(type %in% vImplementedTypes)){
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  if(length(periodsRQ) > length(periods)){
    stop("periodsRQ cannot be longer than periods")
  }
  vDates = as.Date((max(periods)+1):iT, origin = "1970/01/01")
  if(is(RM,"xts")){
    vDates = index(RM)
    vDates = vDates[(max(periods)+1):(iT)]
  }
  ##### Type: "HAR"  - vanilla
  if(type == "HAR"){
  mData = HARDataCreationC(RM, periods)
  coefNames = paste("beta", c(0,periods) , sep="")
  Info = list("Lags" = periods , "type" = type, "dates" = vDates)
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
    mData = HARMatCombine(mA = mData, mB = mJumpData[,-1, drop = FALSE])
    coefNames = c(paste("beta", c(0,periods), sep="") , paste("j" , periodsJ, sep =""))
    
    Info = list("Lags" = periods , "JumpLags" = periodsJ , "type" = type, "dates" = vDates)


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
    
    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])

    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term

    
    coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    Info = list("Lags" = periods , "RQLags" = periodsRQ , "type" = type, "dates" = vDates)
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
    mJumpData = HARDataCreationC(JumpComponent , periodsJ)

    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
    
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term
    mData = HARMatCombine(mData, mJumpData[,-1,drop = FALSE])
    
    coefNames = c(paste("beta", c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""),
                                  paste("beta_j" , periodsJ , sep = ""))
    Info = list("Lags" = periods , "RQLags" = periodsRQ, "JumpLags" = periodsJ, "type" = type, "dates" = vDates)

  }
  #####Type: "HARQ-J" end
  
  if(type == "CHAR"){
    mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
    coefNames = paste("beta", c(0,periods) , sep="")
    Info = list("Lags" = periods , "type" = type, "dates" = vDates)

  }#####Type: "CHAR" end
  
  ##### Type: "CHARQ" - BPQ
  if(type == "CHARQ"){
    mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux

    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)]
    
    
    coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    Info = list("Lags" = periods , "RQLags" = periodsRQ , "type" = type, "dates" = vDates)

    
  }####Type: "CHARQ" - end
  
  if(type == "TV-HAR"){####Type: "TV-HAR" 
    mData = HARDataCreationC(RM, periods)
    mData = cbind(mData[,1], mData[,2] , abs(mData[,2] - mData[,ncol(mData)]) * mData[,2], mData[,3:ncol(mData)])
    coefNames = c("beta", "gamma", "alpha", paste("beta", periods[-1] , sep=""))
    Info = list("Lags" = periods , "type" = type, "dates" = vDates)
  }####Type: "TV-HAR" -end
  
  Model = lm(mData[,1] ~ mData[,-1])
  if(InsanityFilter){
    Model$fitted.values = HARinsanityFilter(Model$fitted.values , 0 , max(RM) , mean(RM) )  
  }
  names(Model$coefficients) = coefNames 
  
  HARModel = new("HARModel" , "Model" = Model ,"Info" = Info)

  return(HARModel)
  
  
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
  
  mData = HARMatCombine(mA = mData, mB = mJumpData[,-1, drop = FALSE])
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
  mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
  
  mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)]
  
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


FASTHARQJestimate = function(RM , RQ , JumpComponent , periods , periodsRQ , periodsJ){
  mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
  meanAux = mean(sqrt(RQ))
  
  mAuxData = sqrt(mAuxData) - meanAux
  mData = HARDataCreationC(RM , periods)
  
  mJumpData = HARDataCreationC(JumpComponent , periodsJ)
  
  mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
  
  mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)]  ##interaction term
  
  mData = HARMatCombine(mData, mJumpData[,-1,drop = FALSE])
  

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



FASTCHARQestimate = function(RM, BPV, RQ, periods, periodsRQ, iT){
  mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods)[,-1])
  
  mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ))
  meanAux = mean(sqrt(RQ))
  
  mAuxData = sqrt(mAuxData) - meanAux
  mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
  mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)]
  
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}



FASTTVHARestimate = function(RM, periods){
  
  mData = HARDataCreationC(RM , periods)
  mData = cbind(mData[,1], mData[,2] , abs(mData[,2] - mData[,ncol(mData)]) * mData[,2], mData[,3:ncol(mData)])
  vCoef = fastLMcoef(cbind(1, mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}

HARinsanityFilter = function(vY , dL, dU , Replacement){
  vIndex = (vY<dL | vY>dU)
  vY[vIndex] = Replacement
  return(vY)
}

