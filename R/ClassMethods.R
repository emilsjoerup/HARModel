setClass("HARModel", representation(Model = "lm", Info = "list", Data = "list"))
setClass("HARForecast" , representation(Model = "HARModel" , Forecast = "data.frame", Info = "list" , Data = "list"))
setClass("HARSim" , representation(Simulation = "numeric" , Info = "list"))

setMethod("show", "HARModel" , function(object) {
  coefficients = object@Model$coefficients
 
  cat("\n----------------------------HARModel-------------------------\n")
  cat("\n Observations used:", length(object@Model$res),"\n", "Maximum lags",max(object@Info$Lags , object@Info$JumpLags), "\n")
  if(object@Info$Type == "HARJ"){
    cat("\n Jump observations used", max(object@Info$JumpLags), "\n")
  }
  cat(paste("\n Specification:"))
  cat("\n Type:" , object@Info$Type)
  cat("\n Lags:", object@Info$Lags)
  if(object@Info$Type == "HARJ"){
    cat("\n JumpLags:", object@Info$JumpLags)  
  }
  cat("\n")
  cat(paste("\n Estimates:\n"))
  cat("\n")
  cat("\n")
  print(round(coefficients , 6) , digits = 5)  
  cat("\n")
  cat("\n Elapsed Time:" , object@Info$ElapsedTime , "seconds\n")
  cat("\n-------------------------------------------------------------\n")
})

setMethod("show" , "HARForecast" , function(object){
  
  cat("\n First model estimated:")
  show(object@Model)
  cat("\n-------------------------HARForecast-------------------------\n")
  cat("\n Forecast specification:\n")
  cat("\n Rolls performed:" , dim(object@Forecast)[1],"\n")
  cat("\n Length of rolls:" , dim(object@Forecast)[2],"\n")
  cat("\n Elapsed time:" , object@Info$ElapsedTime,"seconds\n")
  cat("\n-------------------------------------------------------------\n")
})

setMethod("show" , "HARSim" , function(object){
  coefficients = object@Info$Coefficients
  cat("\n-------------------------HARSim-------------------------\n")
  cat("\n Simulation length:" , object@Info$Length,"\n")
  cat("\n Standard deviation of the error term:", object@Info$ErrorTermSD,"\n")
  cat("\n Lags used:", object@Info$Lags, "\n")
  cat("\n Coefficients:\n")
  print(round(coefficients , 6) , digits = 5)  
  cat("\n")
  cat("\n Elapsed Time:" , object@Info$ElapsedTime , "seconds\n")
  cat("\n--------------------------------------------------------\n")
})

setMethod("plot" , signature(x= "HARModel", y = "missing"), function(x, which=NULL){
  vY = x@Data$`Realized Measure`
  vRes = x@Model$residuals
  vFitted.Val = x@Model$fitted.values
  if(is(vY, "xts")){
    vDates = x@Info$Dates
    print(plot(cbind(vY , vFitted.Val), main = "Realized measure vs. Fitted values" , lty = c(2,1)))
    print(addLegend("topright" , col = c(1,2), on=1 , legend.names = c("Observed RM" , "Fitted Values") , lty = c(2,1) , lwd = c(2,2)))
    invisible(readline(prompt="Press [enter] to continue"))
    print(plot(xts(vRes , order.by = vDates) , main = "Residuals"))
  }
  else{
    vDates = 1:length(vY)
    plot(y= vY ,x = vDates , type="l" , lty = 2 , main = "Realized measure vs. Fitted values" , ylab = "Realized measure" , xlab = "dates")
    lines(vFitted.Val , col = 2)
    legend(x = "topright" , legend = c("Realized Measure" , "Fitted Values") , lty = c(2,1) , col = c(1,2))
    invisible(readline(prompt="Press [enter] to continue"))
    plot(y = vRes , x = vDates , type="l", main = "Residual plot" , ylab = "Residuals" , xlab = "dates")
    }
})

setMethod("plot" , signature(x = "HARForecast" , y = "missing") , function(x , which = NULL){
  if(is(x@Data$Observations , "xts")){
  iRolls = dim(x@Forecast)[2]
  iNAhead = dim(x@Forecast)[1]
  vDatesRoll = tail(x@Data$ForecastDates , iRolls)
  vDatesiNAhead = tail(x@Data$ForecastDates, iNAhead)
  vRollingForecastplot = t(as.vector(x@Forecast[1,]))
  vRollingForecastplot = xts(vRollingForecastplot , order.by = vDatesRoll)
  vForecastResiduals = vRollingForecastplot - x@Data$ForecastComparison[,1]
  vINAheadForecastplot = as.vector(x@Forecast[,1])
  vForecastComp = x@Data$ForecastComparison[,1]
  if(length(vDatesiNAhead)==iNAhead){
    vINAheadForecastplot = xts(vINAheadForecastplot , order.by = vDatesiNAhead)
    print(plot(cbind(vForecastComp, vRollingForecastplot), main="Observed vs. forecasted", ylab = "Realized Volatility", col = c(1,2),  lty=c(2,1)))
    print(addLegend("topleft", on=1, legend.names= c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(2,1) , lwd=c(2,2)))
    invisible(readline(prompt="Press [enter] to continue"))
    print(plot(vForecastResiduals , main = "forecasting residuals"))
  }else{
    print(plot(cbind(vForecastComp, vRollingForecastplot), main="Observed vs. forecasted", ylab = "Realized Volatility" , xlab = "Periods ahead", col = c(1,2),  lty=c(2,1)))
    print(addLegend("topleft", on=1, legend.names= c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(2,1) , lwd=c(2,2)))
    invisible(readline(prompt="Press [enter] to continue"))
    print(plot(vForecastResiduals , main = "forecasting residuals"))
  }
  }
  else{
    iRolls = dim(x@Forecast)[2]
    iNAhead = dim(x@Forecast)[1]
    vDatesRoll = tail(x@Data$ForecastDates , iRolls)
    vDatesiNAhead = tail(x@Data$ForecastDates, iNAhead)
    vRollingForecastplot = t(as.vector(x@Forecast[1,]))
    vForecastResiduals = vRollingForecastplot - x@Data$ForecastComparison
    vINAheadForecastplot = as.vector(x@Forecast[,1])
    vForecastComp = x@Data$ForecastComparison
    if(length(vDatesiNAhead)==iNAhead){
      plot(vForecastComp , main="Observed vs. forecasted", ylab = "Realized Volatility",type = "l" ,col =1,  lty=2)
      lines(vRollingForecastplot , col = 2 )
      legend("topleft", legend =  c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(2,1) , lwd=c(2,2))
      invisible(readline(prompt="Press [enter] to continue"))
      plot(vForecastResiduals , main = "forecasting residuals" , type = "l")
    }else{
      print(plot(vForecastComp , main="Observed vs. forecasted", ylab = "Realized Volatility",type = "l" ,col =1,  lty=2))
      print(lines(vRollingForecastplot , col = 2 ))      
      print(legend("topleft", legend =  c("Observed (out of sample) " , "Forecasted") , col = c(1,2) , lty=c(2,1) , lwd=c(2,2)))
      invisible(readline(prompt="Press [enter] to continue"))
      print(plot(vForecastResiduals , main = "forecasting residuals" , type = "l"))
    }
  }
})

setMethod("plot" , signature(x = "HARSim" , y = "missing") , function(x , length = "ALL" , ctrl = "start"){
  if(length == "ALL"){
  plot(x@Simulation , x = (1:length(x@Simulation)) ,  type = "l" , main = "Simulated HARmodel" , xlab = "Time" , ylab = "Realized Measure")}
  else if (ctrl == "start" & is(length,"numeric")){
    plot(x@Simulation[1:length] , x = (1:length) ,  type = "l" , main = "Simulated HARmodel" , xlab = "Time" , ylab = "Realized Measure")
  }else if(ctrl == "end" & is(length,"numeric")){
    plot(x@Simulation[(length(x@Simulation) - length):length(x@Simulation)] , x = (length(x@Simulation) - length):length(x@Simulation) ,  type = "l" , main = "Simulated HARmodel" , xlab = "Time" , ylab = "Realized Measure")
  }else{
    cat("invalid plot type - if you believe your plot shoud be possible - send me an email")
  }
})

setMethod("coef" , "HARModel" , function(object){
  vCoef = object@Model$coefficients
  return(vCoef)
})

setMethod("coef" , "HARForecast" , function(object){
  vCoef = object@Model@Model$coefficients
  return(vCoef)
})

setMethod("coef" , "HARSim" , function(object){
  vCoef = object@Info$Coefficients
  return(vCoef)
})


setGeneric("uncondmean", function(object)
standardGeneric("uncondmean")
)

setMethod("uncondmean" , "HARModel" , function(object){
  if(object@Info$Type!="HAR"){
    print("Unconditional mean is only implemented for HAR type")
    return(NULL)
  }
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[2:length(vCoef)]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setMethod("uncondmean" , "HARForecast" , function(object){
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[2:length(vCoef)]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setMethod("uncondmean" , "HARSim" , function(object){
  
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[2:length(vCoef)]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setGeneric("SandwichNeweyWest", function(object , lags)
  standardGeneric("SandwichNeweyWest")
)

setMethod("SandwichNeweyWest" , "HARModel" , function(object , lags = 5){
  mVarCovar = sandwich::NeweyWest(object@Model)
  coefficients = object@Model$coefficients
  cat(paste("\n Newey-West Standard errors using a lag order of ", lags , ":\n", sep=""))
  cat(paste(formatC(names(coefficients), width = 8 , format="s")))
  cat("\n")
  cat(paste(formatC(diag(mVarCovar) , width = 8 , format = "fg")))
  cat("\n")
  cat("\n T-Statistics: \n")
  cat(paste(formatC(names(coefficients), width = 8 , format="s")))
  cat("\n", formatC((coefficients/diag(mVarCovar)) , width=8 , format = "fg"))
  cat("\n")
  cat("\n P-Values: \n")
  cat(formatC(names(coefficients) , width=6),"\n" , formatC(dt(coefficients/diag(mVarCovar) , df=object@Model$df.residual) , width=6 , format="f"))
  cat("\n")
  return("HACmatrix" = mVarCovar)
})