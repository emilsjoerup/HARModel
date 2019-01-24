HARsimulate = function(iLength=1500, vLags = c(1, 5, 22), 
                       vCoef = c(0.01, 0.36 ,0.28 , 0.28), 
                       dSigma = 0.001, show = TRUE){
  ######Initialization section ######
  
  start.time = Sys.time()
  iLags = length(vLags)
  iLagsMax = max(vLags)
  ######Initialization end #########

  mSim = HARSimC(iLength+2*iLagsMax, vLags, dConst = vCoef[1], vCoef = vCoef[-1], dSigma = dSigma)

  ElapsedTime = Sys.time() - start.time
  Info = list("Length" = iLength , "Lags" = vLags, "Coefficients" = vCoef, 
              "ErrorTermSD" = dSigma, "ElapsedTime" = ElapsedTime)
  names(Info$Coefficients) = paste("beta", 0:iLags, sep="")
  HARSim = new("HARSim",
               "Simulation" = mSim[(iLagsMax*2+1):(iLagsMax*2 + iLength),1], 
               "Info" = Info)
  if(show) show(HARSim)
  return(HARSim)
  
  
}


