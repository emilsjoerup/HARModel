#Using HAR model Package
rm(list = ls())
#########
library(HARModel)
library(rugarch)
##loading data from the rugarch package.
data("spyreal")
spyreal = spyreal[,2]*100
#########

#Estimate a standard HAR-RV model with the lag-vector of (1,5,22), plots the observations and fitted values, and returns the Newey-West standard errors.

lFit = HARestimate(spyreal, c(1,5,22) , bPlots = T , bStandardErrors = T , iLagSE = 5  , sLegendPlacement = "topright")
summary(lFit)
print(lFit$mVarCovar)
#HARestimate function also returns Newey-West standard-Errors in a manner consistent with Corsi(2009),
#where the auto-correlation is argued to be set to 5. Use iLagSE to change from the standard of 5.
#sLegendPlacement to change the location from the standard of "topright".

#Produce a rolling forecast of the realized measure.
lForecast = HARforecast(vRealizedmeasure = spyreal , vLags = c(1,5,22) ,iNRoll = 501, iNAhead = 247 ,
                        bPlots = T , iPlotForecastAhead = 1 , sLegendPlacement = "topleft")
#HARforecast function produces a forecast consisting of iNRoll+1 rolling windows with a forecast-horizon of iNAhead+1.
#The forecasts are done out of sample, thus the realized measure provided will be estimated with length(vRealizedmeasure)-iNRoll observations.
#A plot of the forecasted values are presented, the standard is to plot the 1-step ahead forecasted values.
#iPlotForecastAhead can be changed to plot h<iNAhead-step ahead forecasts.
#Plots can be disabled by setting bPlots=F, and as above the placement of the legend can be changed by changing
#sLegendPlacement.



