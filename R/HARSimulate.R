HARsimulate = function(iLength=1500, periods = c(1, 5, 22), 
                       coef = c(0.01, 0.36 ,0.28 , 0.28), 
                       Sigma = 0.001){
  ######Initialization section ######
  
  start.time = Sys.time()
  iLags = length(periods)
  iLagsMax = max(periods)
  ######Initialization end #########

  mSim = HARSimC(iLength+2*iLagsMax, periods, dConst = coef[1], coef = coef[-1], dSigma = Sigma)

  ElapsedTime = Sys.time() - start.time
  Info = list("Length" = iLength , "Lags" = periods, "Coefficients" = coef, 
              "ErrorTermSD" = Sigma, "ElapsedTime" = ElapsedTime)
  names(Info$Coefficients) = paste("beta", 0:iLags, sep="")
  HARSim = new("HARSim",
               "Simulation" = mSim[(iLagsMax*2+1):(iLagsMax*2 + iLength),1], 
               "Info" = Info)
  
  return(HARSim)
  
  
}


