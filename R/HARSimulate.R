HARsimulate = function(iLength=10000, vLags = c(1, 5, 22) , vCoef = c(0.002, 0.36 ,0.28 , 0.28), dSigma = 0.001 , type = "HAR" , warning = T){
  ######Initialization section ######
  start.time = Sys.time()
  iLags = length(vLags)
  iLagsPlusOne = iLags+1
  iLagsMax = max(vLags)
  mSim = matrix(nrow = iLength + iLagsMax , ncol = iLagsPlusOne)
  vErrorTermSim = rnorm(iLength+iLagsMax , 0 , sd = dSigma)
  ######Initialization end #########
  if(sum(vCoef[2:iLagsPlusOne])>1 && warning){
    print("Sum of coefficients are above 1 - Watch out for stationarity - proceeding as normal")
  }
  mSim[1:iLagsMax,]  = vCoef[1]/(1-sum(vCoef[2:iLagsPlusOne]))
  for (i in (iLagsMax+1):(iLength + iLagsMax)){
    mSim[(i), ] = HARDataCreationC(mSim[(i-iLagsMax):i,1], vLags)
    mSim[(i),1] =  vCoef[1] + sum(mSim[i , 2:iLagsPlusOne] * vCoef[2:iLagsPlusOne]) + vErrorTermSim[i]
  }
  
  ElapsedTime = Sys.time() - start.time
  Info = list("Length" = iLength , "Lags" = vLags, "Coefficients" = vCoef, "ErrorTermSD" = dSigma ,"ElapsedTime" = ElapsedTime)
  names(Info$Coefficients) = paste("beta", 0:iLags , sep="")
  HARSim = new("HARSim" , "Simulation" = mSim[(iLagsMax+1):(iLagsMax+iLength),1] , "Info" = Info)
  show(HARSim)
  return(HARSim)
}



HARMonteCarlo = function(iLength=1000, vLags = c(1, 5, 22) , vCoef = c(1, 0.36 ,0.28 , 0.28), iBurnin=100 , dSigma = 1 , iLagsPlusOne = length(vLags)+1){
  lSim = HARsimulate(iLength, vLags, vCoef , dSigma)
  vCoef = FASTHARestimate(lSim@Simulation , vLags , iLagsPlusOne)
  return(vCoef)
}


