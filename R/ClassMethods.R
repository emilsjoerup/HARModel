setClass("HARModel", slots = c("Model", "Info"))
setClass("HARForecast", slots = c("Model", "Forecast", "Info", "Data" ))
setClass("HARSim", slots = c("Simulation", "Info"))

setMethod("show", signature(object = "HARModel") , function(object) {
  coefficients = object@Model$coefficients
 
  cat("\n----------------------------------HARModel-------------------------------\n")
  cat("\n Observations used:", length(object@Model$res),"\n", "Maximum lags",max(object@Info$Lags , object@Info$JumpLags), "\n")
  if(object@Info$type == "HARJ"){
    cat("\n Jump observations used", max(object@Info$JumpLags), "\n")
  }
  cat(paste("\n Specification:"))
  cat("\n Type:" , object@Info$type)
  cat("\n Lags:", object@Info$Lags)
  if(object@Info$type == "HARJ"){
    cat("\n JumpLags:", object@Info$JumpLags)  
  }
  if(object@Info$type == "HARQ"){
    cat("\n Realized Quarticity lags:", object@Info$RQLags)  
  }
  cat("\n")
  cat(paste("\n Estimates:\n"))
  cat("\n")
  cat("\n")
  print(round(coefficients , 6) , digits = 4)  
  cat("\n")
  cat("\n")
  cat("R-squared: ", round(summary(object)$r.squared, 4))
  cat("\n")
  cat("Adjusted R-squared: ", round(summary(object)$adj.r.squared, 4))
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("show" , signature(object = "HARForecast") , function(object){
  cat("\n First model estimated:")
  show(object@Model)
  cat("\n-------------------------------HARForecast-------------------------------\n")
  cat("\n Forecast specification:" , object@Info[["windowType"]], "\n")
  cat("\n Length of rolls:" , dim(object@Forecast)[1],"\n")
  cat("\n Rolls performed:" , dim(object@Forecast)[2],"\n")
  cat("\n Elapsed Time:" , round(as.double(object@Info[["ElapsedTime"]] , units = "secs" ) , digits = 3)  , "seconds\n")
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("show" , signature(object = "HARSim") , function(object){
  coefficients = object@Info$Coefficients
  
  cat("\n----------------------------------HARSim---------------------------------\n")
  cat("\n Simulation length:" , object@Info$Length,"\n")
  cat("\n Standard deviation of the error term:", object@Info[["ErrorTermSD"]],"\n")
  cat("\n Lags used:", object@Info$Lags, "\n")
  cat("\n Coefficients:\n")
  print(round(coefficients , 6) , digits = 4)  
  cat("\n")
  cat("\n Elapsed Time:" , round(as.double(object@Info[["ElapsedTime"]] , units = "secs" ) , digits = 3)  , "seconds\n")
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("plot" , signature(x= "HARModel", y = "missing"), 
          function(x, legend.loc = "topright",
                   col = 2:1, lwd= 2, lty = c(1,2),
                   main = NULL,
                   legend.names = c("Realized Measure" , "Fitted values"), yaxis.right = FALSE, ...){
  vY = x@Model$model$`mData[, 1]`
  vFitted.Val = x@Model$fitted.values
  vFitted.Val = xts(vFitted.Val, order.by = x@Info$dates)
  if(is.null(main)) main = paste("Observed vs. fitted based on model: ", x@Info$type)
  plot(cbind(vFitted.Val, vY), main = main, col = col, yaxis.right = yaxis.right, ...)
  addLegend(legend.loc = legend.loc, legend.names = legend.names,
                   col = col, lwd = lwd)
  
})

setMethod("plot" , signature(x = "HARForecast", y = "missing"), 
          function(x, legend.loc = "topright",
                   legend.names = c("Forecasted Values", "Realized Measure"), main =NULL,
                   col = 2:1, lwd= 2, lty = c(1,2), yaxis.right = FALSE,...){
  vForecastComp = x@Data$`ForecastComparison`
  vRollingForecastplot = xts(x@Forecast[1,], index(vForecastComp))
  if(is.null(main)) main = paste("Observed vs. forecasted based on model: ", x@Info$type)
  plot(cbind(vRollingForecastplot,vForecastComp), col = col, main = main, yaxis.right = yaxis.right, ...)
  addLegend(legend.loc = legend.loc, legend.names = legend.names, col = col, lwd = lwd)
  
})

setMethod("plot" , signature(x = "HARSim" , y = "missing"),
          function(x , length = "ALL" , ctrl = "start", main = "Simulated RV", ...){
  vY = xts(x@Simulation , order.by = as.Date(1:length(x@Simulation), origin = "1970/01/01"))
  if(length == "ALL"){
    plot(vY , main = main, ...)
    }
  else if (ctrl == "start" & is(length,"numeric")){
    plot(vY[1:length], main = main, ...)
  }else if(ctrl == "end" & is(length,"numeric")){
    plot(vY[(length(vY) - length):length(vY)], main = main, ... )
  }
  
})

setMethod("coef" , signature(object = "HARModel") , function(object){
  vCoef = object@Model$coefficients
  return(vCoef)
})

setMethod("coef" , signature(object = "HARForecast") , function(object){
  vCoef = object@Model@Model$coefficients
  return(vCoef)
})

setMethod("coef" , signature(object = "HARSim") , function(object){
  vCoef = object@Info$Coefficients
  return(vCoef)
})


setGeneric("uncmean", function(object)
standardGeneric("uncmean")
)

setMethod("uncmean" , signature(object = "HARModel") , function(object){
  if(!object@Info[["type"]]=="HAR"){
    print("Unconditional mean is only implemented for HAR type")
    return(NULL)
  }
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[-1]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setMethod("uncmean" , signature(object = "HARForecast") , function(object){
  if(!object@Info[["type"]]=="HAR"){
    print("Unconditional mean is only implemented for HAR type")
    return(NULL)
  }
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[2:length(vCoef)]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setMethod("uncmean" , signature(object = "HARSim") , function(object){
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[-1]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setGeneric("SandwichNeweyWest", function(object , lags)
  standardGeneric("SandwichNeweyWest")
)

setMethod("SandwichNeweyWest" , signature(object = "HARModel") , function(object , lags = 5){
  if(!object@Info$type=="HAR"){
    print("SandwichNeweyWest is only implemented for HAR type,  you can use sandwich::NeweyWest() on the lm submodel")
    return(NULL)
  }
  if(missingArg(lags)){
    lags = 5
  }
  cat("\n-------------------Newey-West Standard errors----------------\n")
  mVarCovar = sandwich::NeweyWest(object@Model , lags)
  coefficients = object@Model$coefficients
  cat(paste("\n Newey-West Standard errors using a lag order of ", lags , ":\n", sep=""))
  mPrint = rbind(coefficients , "Standard errors" =  diag(mVarCovar), 
                 "T-Statistics" = coefficients/diag(mVarCovar), "P-values" = dt(coefficients/diag(mVarCovar) , 
                                                                                df=object@Model$df.residual))
  print(round(mPrint,6 ) ,5 )
  cat("\n-------------------------------------------------------------\n")
  return("HACmatrix" = mVarCovar)
})

setGeneric("forc", function(object, WhichStep = 1)
standardGeneric("forc")  
)

setMethod("forc", signature(object = "HARForecast"), function(object, WhichStep = 1){
  vDates = object@Data$ForecastDates + WhichStep-1
  vForc = xts(object@Forecast[WhichStep,], vDates)
  return(vForc)
})

setGeneric("forecastres", function(object)
standardGeneric("forecastres")
)


setMethod("forecastres", signature(object = "HARForecast"), function(object){
  vRes = object@Forecast[1,] - object@Data$ForecastComparison
  return(vRes)
})


setGeneric("qlike", function(object)
  standardGeneric("qlike")
)

setMethod("qlike", signature(object = "HARModel"), function(object){
  RM = object@Model$model[,1] #extract RM from model
  FV  = object@Model$fitted.values
  qLike = RM/FV - log(RM / FV) - 1
  return(qLike) 
}
)

setMethod("qlike", signature(object = "HARForecast"), function(object){
  RM = object@Data$ForecastComparison  #extract observed RM
  FV  = object@Forecast[1,]
  qLike = RM/FV - log(RM / FV) - 1
  return(qLike) 
}
)

###Methods that work with "lm" objects that I thought may be useful.
###Wrappers from the HARModel to "lm"
setMethod("logLik" , signature(object = "HARModel") , function(object, ...){
  out = logLik(object@Model, ...)
  return(out)
})

setMethod("confint" , signature(object = "HARModel") , function(object, parm, level = 0.95, ...){
  out = confint(object@Model, parm = parm, level = level, ...)
  return(out)
})

setMethod("residuals" , signature(object = "HARModel") , function(object, ...){
  out = residuals(object@Model, ...)
  return(out)
})

setMethod("summary" , signature(object = "HARModel") , function(object, ...){
  out = summary(object@Model,...)
  out$call =  as.name("lm(y ~ x)")
  return(out)
})

