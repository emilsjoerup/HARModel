HAREstimate = function(RM, BPV = NULL, RQ = NULL, periods = c(1,5,22), periodsJ = NULL, periodsRQ = NULL,
                       type = "HAR" , insanityFilter = TRUE, h = 1){
  
  iT = length(RM)
  Info = list()
  vImplementedTypes = c("HAR" , "HARJ" , "HARQ" , "HARQ-J", "CHAR", "CHARQ", "TV-HAR")
  if(!(type %in% vImplementedTypes)){
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  if(length(periodsRQ) > length(periods)){
    stop("periodsRQ cannot be longer than periods")
  }
  vDates = as.Date((max(periods)+h):iT, origin = "1970/01/01")
  if(is(RM,"xts")){
    vDates = index(RM)
    vDates = vDates[(max(periods)+h):(iT)]
  }
  if(h %% 1 != 0 | h < 1){
    stop("h must be a positive integer")
  }
  
  ##### Type: "HAR"  - vanilla
  if(type == "HAR"){
    mData = HARDataCreationC(RM, periods, h)
    coefNames = paste("beta", c(0,periods) , sep="")
    Info = list("periods" = periods , "type" = type, "dates" = vDates)
  }
  ##### Type: "HAR"  - vanilla end
  
  
  ##### Type: "HARJ" 
  if(type == "HARJ"){
    mData = HARDataCreationC(RM, periods, h)
    JumpComponent = pmax(RM - BPV,0)
    if(is.null(periodsJ)){
      periodsJ = periods
    }
    mJumpData = HARDataCreationC(JumpComponent , periodsJ,  h = h)
    mData = HARMatCombine(mA = mData, mB = mJumpData[,-1, drop = FALSE])
    coefNames = c(paste("beta", c(0,periods), sep="") , paste("j" , periodsJ, sep =""))
    
    Info = list("periods" = periods , "preiodsJ" = periodsJ , "type" = type, "dates" = vDates)
    
    
  }
  ##### Type: "HARJ" end
  
  ##### Type: "HARQ"
  if(type == "HARQ"){
    mData = HARDataCreationC(RM, periods, h)
    if(is.null(periodsRQ)){
      periodsRQ = periods 
    }
    if(h!=1 & length(periodsRQ) == 1){
      if(!h %in% periods){
        stop("h is not in periods, that is not allowed as of now when periodsRQ is of length 1")
      }
      mAuxData = matrix(1, nrow = length(RQ) - 2*h +1, ncol = length(periods))
      mAuxData[,h == periods] = as.matrix(HARDataCreationC(RQ , c(h), h))[,-1]
      meanAux = mean(sqrt(RQ))
      mAuxData[, h == periods] = sqrt(mAuxData[,h == periods]) - meanAux
      
      mData = HARMatCombine(mData, mAuxData[,h == periods, drop = FALSE])
      mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * 
        mData[,c(FALSE, h==periods, rep(FALSE, length(periodsRQ)))] ##interaction term
      coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , h, sep =""))
    }else{
      mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ, h))
      meanAux = mean(sqrt(RQ))
      mAuxData = sqrt(mAuxData) - meanAux
      
      mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
      
      mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term
      coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    }
    
    Info = list("periods" = periods , "periodsRQ" = periodsRQ , "type" = type, "dates" = vDates)
    
  }
  ##### Type: "HARQ" end
  
  ##### Type: "HARQ-J"
  if(type == "HARQ-J"){
    mData = HARDataCreationC(RM, periods, h)
    JumpComponent = pmax(RM - BPV,0)
    if(is.null(periodsRQ)){
      periodsRQ = periods 
    }
    if(is.null(periodsJ)){
      periodsJ = periods
    }
    
    
    if(h!=1 & length(periodsRQ) == 1){
      if(!h %in% periods){
        stop("h is not in periods, that is not allowed as of now when periodsRQ is of length 1")
      }
      mAuxData = matrix(1, nrow = length(RQ) - 2*h +1, ncol = length(periods))
      mAuxData[,h == periods] = as.matrix(HARDataCreationC(RQ , c(h), h))[,-1]
      meanAux = mean(sqrt(RQ))
      mAuxData[, h == periods] = sqrt(mAuxData[,h == periods]) - meanAux
      
      mData = HARMatCombine(mData, mAuxData[,h == periods, drop = FALSE])
      mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * 
        mData[,c(FALSE, h==periods, rep(FALSE, length(periodsRQ)))] ##interaction term
      coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , h, sep =""))
    }else{
      mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ, h))
      meanAux = mean(sqrt(RQ))
      mAuxData = sqrt(mAuxData) - meanAux
      
      mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
      
      mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term
      coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    }
    

    mJumpData = HARDataCreationC(JumpComponent , periodsJ, h)
    mData = HARMatCombine(mData, mJumpData[,-1,drop = FALSE])
    
    coefNames = c(coefNames, paste("beta_j" , periodsJ , sep = ""))
    Info = list("periods" = periods , "periodsRQ" = periodsRQ, "periodsJ" = periodsJ, "type" = type, "dates" = vDates)
    
  }
  #####Type: "HARQ-J" end
  
  if(type == "CHAR"){
    mData = cbind(RM[(max(periods)+h):(iT)], HARDataCreationC(BPV, periods, h)[,-1]) #Aggregate the explanatory variables
    if(h != 1){
      mData[,1] = HARDataCreationC(RM, periods, h)[,1] #Aggregate the dependent variable too
    }
    coefNames = paste("beta", c(0,periods) , sep="")
    Info = list("periods" = periods , "type" = type, "dates" = vDates)
    
  }#####Type: "CHAR" end
  
  ##### Type: "CHARQ" - BPQ
  if(type == "CHARQ"){
    mData = cbind(RM[(max(periods)+h):(iT)], HARDataCreationC(BPV, periods, h)[,-1])
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ, h))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)]
    
    
    coefNames = c(paste("beta",c(0,periods), sep="") , paste("beta_q" , periodsRQ, sep =""))
    Info = list("Lags" = periods , "RQLags" = periodsRQ , "type" = type, "dates" = vDates)
    
    
  }####Type: "CHARQ" - end
  
  ####Type: "TV-HAR" 
  if(type == "TV-HAR"){
    mData = HARDataCreationC(RM, periods, h)
    mData = cbind(mData[,1], mData[,2] , abs(mData[,2] - mData[,ncol(mData)]) * mData[,2], mData[,3:ncol(mData)])
    coefNames = c("beta", "gamma", "alpha", paste("beta", periods[-1] , sep=""))
    Info = list("Lags" = periods , "type" = type, "dates" = vDates)
  }####Type: "TV-HAR" -end
  
  model = lm(mData[,1] ~ mData[,-1])
  if(insanityFilter){
    model$fitted.values = HARinsanityFilter(model$fitted.values , 0 , max(mData[,1]) , mean(mData[,1]) )  
  }
  names(model$coefficients) = coefNames 
  Info$h = h
  HARModel = new("HARModel" , "model" = model ,"info" = Info)
  
  return(HARModel)
  
}


#######Fast estimation routine which returns only the coefficients, it is used for forecasting and cuts down the time a lot.
#######These functions are undocumented and not exported

FASTHARestimate = function(RM , periods, h = 1){
  
  mData = HARDataCreationC(RM , periods, h)
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
}

FASTHARJestimate = function(RM, JumpComponent, periods, periodsJ, h = 1 ){
  mData = HARDataCreationC(RM , periods, h)
  
  mJumpData = HARDataCreationC(JumpComponent , periodsJ, h)
  
  mData = HARMatCombine(mA = mData, mB = mJumpData[,-1, drop = FALSE])
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1]) 
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}

FASTHARQestimate =function(RM, RQ, periods, periodsRQ = NULL, h = 1){
  if(is.null(periodsRQ)){
    periodsRQ = periods
  }
  
  mData = HARDataCreationC(RM , periods , h)
  
  if(h!=1 & length(periodsRQ) == 1){
    mAuxData = matrix(1, nrow = length(RQ) - 2*h +1, ncol = length(periods))
    mAuxData[,h == periods] = as.matrix(HARDataCreationC(RQ , c(h), h))[,-1]
    meanAux = mean(sqrt(RQ))
    mAuxData[, h == periods] = sqrt(mAuxData[,h == periods]) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,h == periods, drop = FALSE])
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * 
      mData[,c(FALSE, h==periods, rep(FALSE, length(periodsRQ)))] ##interaction term
    
  }else{
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ, h))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
    
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term
    
  }

  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


FASTHARQJestimate = function(RM , RQ , JumpComponent , periods , periodsRQ , periodsJ, h = 1){
  mData = HARDataCreationC(RM , periods, h)
  mJumpData = HARDataCreationC(JumpComponent , periodsJ, h)
  if(h!=1 & length(periodsRQ) == 1){
    mAuxData = matrix(1, nrow = length(RQ) - 2*h +1, ncol = length(periods))
    mAuxData[,h == periods] = as.matrix(HARDataCreationC(RQ , c(h), h))[,-1]
    meanAux = mean(sqrt(RQ))
    mAuxData[, h == periods] = sqrt(mAuxData[,h == periods]) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,h == periods, drop = FALSE])
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * 
      mData[,c(FALSE, h==periods, rep(FALSE, length(periodsRQ)))] ##interaction term
    
  }else{
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ, h))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
    
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term
    
  }
  
  mData = HARMatCombine(mData, mJumpData[,-1,drop = FALSE])
  
  
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}


FASTCHARestimate = function(RM, BPV, periods, iT, h = 1){
  mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods, h)[,-1])
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
}



FASTCHARQestimate = function(RM, BPV, RQ, periods, periodsRQ, iT, h = 1){
  mData = cbind(RM[(max(periods)+1):(iT)], HARDataCreationC(BPV, periods, h)[,-1])
  if(h!=1 & length(periodsRQ) == 1){
    mAuxData = matrix(1, nrow = length(RQ) - 2*h +1, ncol = length(periods))
    mAuxData[,h == periods] = as.matrix(HARDataCreationC(RQ , c(h), h))[,-1]
    meanAux = mean(sqrt(RQ))
    mAuxData[, h == periods] = sqrt(mAuxData[,h == periods]) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,h == periods, drop = FALSE])
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * 
      mData[,c(FALSE, h==periods, rep(FALSE, length(periodsRQ)))] ##interaction term
    
  }else{
    mAuxData = as.matrix(HARDataCreationC(RQ , periodsRQ, h))
    meanAux = mean(sqrt(RQ))
    mAuxData = sqrt(mAuxData) - meanAux
    
    mData = HARMatCombine(mData, mAuxData[,-1, drop = FALSE])
    
    mData[,-(1:(length(periods) +1))] = mData[,-(1:(length(periods) +1))] * mData[,2:(length(periodsRQ)+1)] ##interaction term
    
  }
  vCoef = fastLMcoef(cbind(1,mData[-nrow(mData),-1]) , mData[-nrow(mData),1])
  #leave the last row out of the estimation and pass it back to form the forecast
  return(list("coefficients" = vCoef, "mData" = mData))
  
}



FASTTVHARestimate = function(RM, periods, h = 1){
  
  mData = HARDataCreationC(RM , periods, h)
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

