HARSimulate = function(len=1500, periods = c(1, 5, 22), 
                       coef = c(0.01, 0.36 ,0.28 , 0.28), 
                       errorTermSD = 0.001){
  ######Initialization section ######
  
  start.time = Sys.time()
  iLags = length(periods)
  iLagsMax = max(periods)
  ######Initialization end #########

  mSim = HARSimC(len+2*iLagsMax, periods, dConst = coef[1], coef = coef[-1], dSigma = errorTermSD)
  ElapsedTime = Sys.time() - start.time
  info = list("len" = len , "periods" = periods, "coefficients" = coef, 
              "errorTermSD" = errorTermSD, "elapsedTime" = ElapsedTime)
  names(info$coefficients) = paste("beta", 0:iLags, sep="")
  HARSim = new("HARSim",
               "simulation" = mSim[(iLagsMax*2+1):(iLagsMax*2 + len),1], 
               "info" = info)
  
  return(HARSim)
  
  
}


