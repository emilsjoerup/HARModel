print.HARmodel = function(HARmodel){
  with(HARmodel , 
       modelterms = unlist(terms) ,
    cat("\nFormula" , "\n" , terms[[1]] , terms[[2]] , terms[[3]] , "+" , paste(terms[[4]][1:(length(terms[[4]])-1)] ,"", sep=" +") ,terms[[4]][length(terms[[4]])], "\n",
   "\nCoefficients","\n",formatC(names(coefficients) , width = 6),"\n", formatC(coefficients,width=6 , format="f"),"\n" ,
   "\nNewey-West standard errors" , "\n" , formatC(names(coefficients) , width = 6) , "\n" , formatC(diag(mVarCovar) , width=6 , format = "f") , "\n",
   "\nT-stats", "\n" ,formatC(names(coefficients) , width=6),"\n" , formatC((coefficients/diag(mVarCovar)) , width=6), "\n",
   "\nP-values", "\n" , formatC(names(coefficients) , width=6),"\n" , formatC(dt(coefficients/diag(mVarCovar) , df=df.residual) , width=6 , format="f"), "\n",
   sep="  ")
)
  invisible(HARmodel)
}

print.HARforecast = function(HARforecast){
  with(HARforecast , 
       modelterms = unlist(terms) ,
       cat("\nHARmodel for the first" , length(Observations) , "observations",
      "\nFormula" , "\n" , terms[[1]] , terms[[2]] , terms[[3]] , "+" , paste(terms[[4]][1:(length(terms[[4]])-1)] ,"", sep=" +") ,terms[[4]][length(terms[[4]])], "\n",
      "\nCoefficients","\n",formatC(names(coefficients) , width = 6),"\n", formatC(coefficients,width=6 , format="f"),"\n" ,
      "\nNewey-West standard errors" , "\n" , formatC(names(coefficients) , width = 6) , "\n" , formatC(diag(mVarCovar) , width=6 , format = "f") , "\n",
      "\nT-stats", "\n" ,formatC(names(coefficients) , width=6),"\n" , formatC((coefficients/diag(mVarCovar)) , width=6), "\n",
      "\nP-values", "\n" , formatC(names(coefficients) , width=6),"\n" , formatC(dt(coefficients/diag(mVarCovar) , df=df.residual) , width=6 , format="f"), "\n",
      "\n",
      "\n----------------------------------------------",
      "\nRolling forecasts",dim(ForecastMatrix)[2],
      "\nHorizon per roll" , dim(ForecastMatrix)[1],
      "\n----------------------------------------------", "\n",
      sep="  ")
  )
  invisible(HARforecast)
} 
  
plot.HARmodel = function(HARmodel ,which = c(1:5), caption = list("Residuals vs Fitted", 
"Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage", 
expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))), sLegendPlacement = "topleft", 
panel = if (add.smooth) panel.smooth else points, sub.caption = NULL, 
main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(), 
..., id.n = 3, labels.id = names(residuals(HARmodel)), cex.id = 0.75, 
qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"), 
label.pos = c(4, 2), cex.caption = 1){
      ##For now, the plot function is largely the same as the plot.harModel function from the highfrequency package.
      #browser()
      vObservations = HARmodel$model$`mData[, 1]`
      vFittedValues = HARmodel$fitted.values
      vDates    = HARmodel$dates
      vDates    = as.POSIXct(vDates)
      vObservations = xts(vObservations, order.by=vDates)
      vFittedValues   = xts(vFittedValues, order.by=vDates)
      g_range = range(vFittedValues,vObservations)
      g_range[1] = 0.95*g_range[1] 
      g_range[2]= 1.05 * g_range[2] 
      title = paste("Observed and forecasted RV based on HAR Model:")
      plot(cbind(vObservations, vFittedValues), col=c(1:2), main=title,ylab="Realized Volatility", lty=c(2,1))
      addLegend(sLegendPlacement, on=1, 
                legend.names = c("Observed RV","Forecasted RV"), 
                lty=c(1, 1), lwd=c(2, 2),
                col=c(1:2))
      
}

plot.HARforecast = function(HARforecast, iPlotForecastAhead = 1 , sLegendPlacement = "topleft"){
  iRolls = dim(HARforecast$ForecastMatrix)[2]
  iNAhead = dim(HARforecast$ForecastMatrix)[1]
  vDatesRoll = tail(HARforecast$ForecastDates , iRolls)
  vDatesiNAhead = vDatesRoll[1:iNAhead]
  vRollingForecastplot = t(as.vector(HARforecast$ForecastMatrix[iPlotForecastAhead,]))
  vRollingForecastplot = xts(vRollingForecastplot , order.by = vDatesRoll)
  
  vForecastResiduals = vRollingForecastplot - HARforecast$vForecastComp
  
  vINAheadForecastplot = as.vector(HARforecast$ForecastMatrix[,1])
  vINAheadForecastplot = xts(vINAheadForecastplot , order.by = vDatesiNAhead)
  
  
  print(plot(cbind(HARforecast$vForecastComp, vRollingForecastplot), main="Observed vs. forecasted", ylab = "Realized Volatility" , xlab = "Periods ahead"))
  print(addLegend(sLegendPlacement, on=1, legend.names= c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(1,1) , lwd=c(2,2)))
  
  invisible(readline(prompt="Press [enter] to continue"))
  
  print(plot(vForecastResiduals , main = "forecast-residuals"))
  
  invisible(readline(prompt="Press [enter] to continue"))
  print(plot(vINAheadForecastplot , main = "h-step ahead forecast"))
}
  
  
