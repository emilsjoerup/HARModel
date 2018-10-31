HARsimulate = function(iLength=1500, vLags = c(1, 5, 22) , vCoef = c(0.002, 0.36 ,0.28 , 0.28), dSigma = 0.001 , type = "HAR"){
  ######Initialization section ######
  start.time = Sys.time()
  iLags = length(vLags)
  iLagsPlusOne = iLags+1
  iLagsMax = max(vLags)
  mSim = matrix(nrow = iLength + iLagsMax , ncol = iLagsPlusOne)
  vErrorTermSim = rnorm(iLength+iLagsMax , 0 , sd = dSigma)
  ######Initialization end #########
  vImplementedTypes = c("HAR")
  
  if(!any(grepl(type, vImplementedTypes))){
    
    cat("type argument not correctly specifiec or is not implemented, available types are:", paste(dQuote(vImplementedTypes)))
    return(NULL)
  }
  mSim[1:iLagsMax,]  = vCoef[1]/(1-sum(vCoef[-1]))
  for (i in (iLagsMax+1):(iLength + iLagsMax)){
    mSim[(i), ] = HARDataCreationC(mSim[(i-iLagsMax):i,1], vLags)
    mSim[(i),1] =  vCoef[1] + sum(mSim[i , 2:iLagsPlusOne] * vCoef[-1]) + vErrorTermSim[i]
  }
  
  ElapsedTime = Sys.time() - start.time
  Info = list("Length" = iLength , "Lags" = vLags, "Coefficients" = vCoef, "ErrorTermSD" = dSigma ,"ElapsedTime" = ElapsedTime)
  names(Info$Coefficients) = paste("beta", 0:iLags , sep="")
  HARSim = new("HARSim" , "Simulation" = mSim[(iLagsMax+1):(iLagsMax+iLength),1] , "Info" = Info)
  show(HARSim)
  return(HARSim)
}


