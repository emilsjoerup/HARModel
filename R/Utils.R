print.HARmodel = function(object){
  UseMethod("HARmodel" , object)
}


setMethod("print.HARmodel" , signature(object="HARmodel" ) , function(object){
  with(object , 
       modelterms = unlist(terms) ,
       cat("\nFormula" , "\n" , terms[[1]] , terms[[2]] , terms[[3]] , "+" , paste(terms[[4]][1:(length(terms[[4]])-1)] ,"", sep=" +") ,terms[[4]][length(terms[[4]])], "\n",
           "\nCoefficients","\n",formatC(names(coefficients) , width = 6),"\n", formatC(coefficients,width=6 , format="f"),"\n" ,
           "\nNewey-West standard errors" , "\n" , formatC(names(coefficients) , width = 6) , "\n" , formatC(diag(mVarCovar) , width=6 , format = "f") , "\n",
           "\nT-stats", "\n" ,formatC(names(coefficients) , width=6),"\n" , formatC((coefficients/diag(mVarCovar)) , width=6), "\n",
           "\nP-values", "\n" , formatC(names(coefficients) , width=6),"\n" , formatC(dt(coefficients/diag(mVarCovar) , df=df.residual) , width=6 , format="f"), "\n",
           sep="  ")
  )
  invisible(object)
})

      

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
plot.HARmodel = function(object){
  UseMethod("HARmodel" , object)
  NextMethod(object = object)
  
}
setMethod("plot.HARmodel" , signature(object = "HARmodel"), function(object){
  ##For now, the plot function is largely the same as the plot.harModel function from the highfrequency package.
  vObservations = object$model$`mData[, 1]`
  vFittedValues = object$fitted.values
  vDates    =    object$dates
  vDates    = as.POSIXct(vDates)
  vObservations = xts(vObservations, order.by=vDates)
  vFittedValues   = xts(vFittedValues, order.by=vDates)
  g_range = range(vFittedValues,vObservations)
  g_range[1] = 0.95*g_range[1] 
  g_range[2]= 1.05 * g_range[2] 
  title = paste("Observed and forecasted RV based on HAR Model:")
  plot(cbind(vObservations, vFittedValues), col=c(1:2), main=title,ylab="Realized Volatility", lty=c(2,1))
  addLegend("topleft", on=1, 
            legend.names = c("Observed RV","Fitted values"), 
            lty=c(2, 1), lwd=c(2, 2),
            col=c(1:2))
})







plot.HARforecast = function(object){
  UseMethod("HARforecast" , object)
}
  
  
setMethod("plot.HARforecast" , signature(object = "HARforecast")  , function(object){
  iRolls = dim(object$ForecastMatrix)[2]
  iNAhead = dim(object$ForecastMatrix)[1]
  vDatesRoll = tail(object$ForecastDates , iRolls)
  vDatesiNAhead = tail(object$ForecastDates, iNAhead)
  vRollingForecastplot = t(as.vector(object$ForecastMatrix[1,]))
  vRollingForecastplot = xts(vRollingForecastplot , order.by = vDatesRoll)
  
  vForecastResiduals = vRollingForecastplot - object$vForecastComp
  vINAheadForecastplot = as.vector(object$ForecastMatrix[,1])
  
  if(length(vDatesiNAhead)==iNAhead){
    vINAheadForecastplot = xts(vINAheadForecastplot , order.by = vDatesiNAhead)
    print(plot(cbind(object$vForecastComp, vRollingForecastplot), main="Observed vs. forecasted", ylab = "Realized Volatility", col = c(1,2),  lty=c(2,1)))
    print(addLegend("topleft", on=1, legend.names= c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(2,1) , lwd=c(2,2)))
    
    invisible(readline(prompt="Press [enter] to continue"))
    
    print(plot(vForecastResiduals , main = "forecasting residuals"))
    
    invisible(readline(prompt="Press [enter] to continue"))
    print(plot(vINAheadForecastplot , main = "First roll forecast"))
    
  }else{
    print(plot(cbind(object$vForecastComp, vRollingForecastplot), main="Observed vs. forecasted", ylab = "Realized Volatility" , xlab = "Periods ahead", col = c(1,2),  lty=c(2,1)))
    print(addLegend("topleft", on=1, legend.names= c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(2,1) , lwd=c(2,2)))
    
    invisible(readline(prompt="Press [enter] to continue"))
    
    print(plot(vForecastResiduals , main = "forecasting residuals"))
    
    invisible(readline(prompt="Press [enter] to continue"))
    print(plot(vINAheadForecastplot , main = "First roll forecast" , type="l"))
  }
  invisible(object)
})
  
  
