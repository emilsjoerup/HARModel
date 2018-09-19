########################################################################################################
# This package is created by Emil Sjørup, at the time of beginning a bachelor of economics student
# at the university of Aarhus in Denmark.
# Any bugs should be reported to Emilsjoerup@live.dk  
########################################################################################################

HARsimulate = function(iLength=1000, vLags = c(1, 5, 22) , vCoef = c(1, 0.36 ,0.28 , 0.28), iBurnin=100 , dSigma = 1){
  ######Initialization section ######
  
  start.time = Sys.time()
  iLags = length(vLags)
  iLagsPlusOne = iLags+1
  iLagsMax = max(vLags)
  vBurnin = rnorm(iLagsMax+1)^2
  mBurnin = matrix(1 ,nrow =iBurnin , ncol = iLagsPlusOne )
  mSim = matrix(NA , nrow = iLength , ncol= iLagsPlusOne)
  if(sum(vCoef[2:iLagsPlusOne])>0.95){
    print("Sum of coefficients are close to or above 1 - Watch out for stationarity.")
  }
  ######Initialization end    #######
  #browser()
  mBurnin[1,] = HARDataCreationC(vBurnin , vLags)
  #Burnin loop#
  for (i in 2:iBurnin) {
    mBurnin[(i),] = tail(HARDataCreationC(c(vBurnin, mBurnin[1:(i),1]) , vLags),1)
    mBurnin[(i),1] = max(0, vCoef[1] + sum(mBurnin[(i),2:iLagsPlusOne] * vCoef[2:iLagsPlusOne]) + rnorm(1 , sd = dSigma))
  }# end burnin loop
  #Simulation loop
  mSim[1,]  = tail(HARDataCreationC(mBurnin[,1] , vLags),1)
  for (i in 1:(iLagsMax)) {
    mSim[(i), ] = tail(HARDataCreationC(c(mBurnin[(iBurnin-iLagsMax + i):iBurnin, 1] , mSim[1:i,1]) , vLags),1 )
    mSim[(i),1] = max(0 , vCoef[1] + sum(mSim[i , 2:iLagsPlusOne] * vCoef[2:iLagsPlusOne]) + rnorm(1,sd =dSigma)) 
  }
  
  
  for (i in (iLagsMax+1):iLength){
    #browser()
    mSim[(i), ] = tail(HARDataCreationC(mSim[(i-iLagsMax):i,1], vLags),1 )
    #print(length(mSim[(i-iLagsMax) : i,1]))
    mSim[(i),1] = max(0 , vCoef[1] + sum(mSim[i , 2:iLagsPlusOne] * vCoef[2:iLagsPlusOne]) + rnorm(1, sd = dSigma))
  }
  
  ElapsedTime = Sys.time() - start.time
  Info = list("Length" = iLength , "Lags" = vLags, "Coefficients" = vCoef, "ErrorTermSD" = dSigma , "BurninLength" = iBurnin , "ElapsedTime" = ElapsedTime)
  
  names(Info$Coefficients) = paste("beta", 0:iLags , sep="")
  
  
  
  HARSim = new("HARSim" , "Simulation" = mSim[,1] , "Info" = Info)
  
  return(HARSim)
}


HARMonteCarlo = function(iLength=1000, vLags = c(1, 5, 22) , vCoef = c(1, 0.36 ,0.28 , 0.28), iBurnin=100 , dSigma = 1 , iLagsPlusOne = length(vLags)+1){
  lSim = HARsimulate(iLength, vLags, vCoef , iBurnin , dSigma)
  vCoef = FASTHARestimate(lSim@Simulation , vLags , iLagsPlusOne)
  return(vCoef)
}