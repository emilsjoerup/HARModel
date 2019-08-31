setClass("HARModel", slots = c("model", "info"))
setClass("HARForecast", slots = c("model", "forecast", "info", "data" ))
setClass("HARSim", slots = c("simulation", "info"))

setMethod("show", signature(object = "HARModel") , function(object) {
  coefficients = coef(object)
 
  cat("\n----------------------------------HARModel-------------------------------\n")
  cat("\n Observations used:", length(object@model$res),"\n", "Maximum lags",max(object@info$periods , object@info$periodsJ), "\n")
  if(object@info$type == "HARJ"){
    cat("\n Jump observations used", max(object@info$periodsJ), "\n")
  }
  cat(paste("\n Specification:"))
  cat("\n Type:" , object@info$type)
  cat("\n Lags:", object@info$periods)
  if(object@info$type == "HARJ"){
    cat("\n Lags used for the jump component:", object@info$periodsJ)  
  }
  if(object@info$type == "HARQ"){
    cat("\n Realized Quarticity lags:", object@info$periodsRQ)  
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
  show(object@model)
  cat("\n-------------------------------HARForecast-------------------------------\n")
  cat("\n Forecast specification:" , object@info[["windowType"]], "\n")
  cat("\n Length of rolls:" , dim(object@forecast)[1],"\n")
  cat("\n Rolls performed:" , dim(object@forecast)[2],"\n")
  cat("\n Elapsed Time:" , round(as.double(object@info[["elapsedTime"]] , units = "secs" ) , digits = 3)  , "seconds\n")
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("show" , signature(object = "HARSim") , function(object){
  coefficients = coef(object)
  cat("\n----------------------------------HARSim---------------------------------\n")
  cat("\n Simulation length:" , object@info$Length,"\n")
  cat("\n Standard deviation of the error term:", object@info[["errorTermSD"]],"\n")
  cat("\n Lags used:", object@info$periods, "\n")
  cat("\n Coefficients:\n")
  print(round(coefficients , 6) , digits = 4)  
  cat("\n")
  cat("\n Elapsed Time:" , round(as.double(object@info[["elapsedTime"]] , units = "secs" ) , digits = 3)  , "seconds\n")
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("plot" , signature(x= "HARModel", y = "missing"), 
          function(x, legend.loc = "topright",
                   col = 1:2, lwd= 2, lty = c(1,2),
                   main = NULL,
                   legend.names = c("Realized Measure" , "Fitted values"), yaxis.right = FALSE, ...){
  
  vY = xts(x@model$model$`mData[, 1]`, order.by = x@info$dates)
  vFitted.Val = xts(as.numeric(x@model$fitted.values), order.by = index(vY))
  if(is.null(main)) {
    main = paste("Observed vs. fitted based on model: ", x@info$type)
    h = x@info$h
    if(h!=1){
      main = paste(main, "\n h=", h)
    }
  }
  plot(vY, main = main, col = col[1], yaxis.right = yaxis.right, ...)
  lines(vFitted.Val, col = col[2])
  addLegend(legend.loc = legend.loc, legend.names = legend.names,
                   col = col, lwd = lwd)
  
})

setMethod("plot" , signature(x = "HARForecast", y = "missing"), 
          function(x, legend.loc = "topright",
                   legend.names = c("Realized Measure","Forecasted Values"), main =NULL,
                   col = 1:2, lwd= 2, lty = c(1,2), yaxis.right = FALSE,...){
  vForecastComp = x@data$`forecastComparison`
  vRollingForecastplot = getForc(x)
  if(is.null(main)) main = paste("Observed vs. forecasted based on model: ", x@info$type)
  plot(vForecastComp, col = col[1], main = main, yaxis.right = yaxis.right, ...)
  lines(vRollingForecastplot, col = col[2])
  addLegend(legend.loc = legend.loc, legend.names = legend.names, col = col, lwd = lwd)
  
})

setMethod("plot" , signature(x = "HARSim" , y = "missing"),
          function(x , length = "ALL" , ctrl = "start", main = "Simulated RV", ...){
  vY = xts(x@simulation , order.by = as.Date(1:length(x@simulation), origin = "1970/01/01"))
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
  vCoef = object@model$coefficients
  return(vCoef)
})

setMethod("coef" , signature(object = "HARForecast") , function(object){
  vCoef = object@model@model$coefficients
  return(vCoef)
})

setMethod("coef" , signature(object = "HARSim") , function(object){
  vCoef = object@info$Coefficients
  return(vCoef)
})


setGeneric("uncmean", function(object)
  standardGeneric("uncmean")
)

setMethod("uncmean" , signature(object = "HARModel") , function(object){
  if(object@info[["type"]]!="HAR"){
    print("Unconditional mean is only implemented for HAR type")
    return(NULL)
  }
  vCoef = coef(object)
  uncmean = vCoef[1]/(1-sum(vCoef[-1]))
  names(uncmean) = ""
  return(c("Unconditional Mean" = uncmean))
})

setMethod("uncmean" , signature(object = "HARForecast") , function(object){
  if(!object@info[["type"]]=="HAR"){
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

setGeneric("sandwichNeweyWest", function(object , lags)
  standardGeneric("sandwichNeweyWest")
)

setMethod("sandwichNeweyWest" , signature(object = "HARModel") , function(object , lags = 5){
  if(missingArg(lags)){
    lags = 5
  }
  cat("\n-------------------Newey-West Standard errors----------------\n")
  mVarCovar = sandwich::NeweyWest(object@model , lags)
  coefficients = coef(object)
  cat(paste("\n Newey-West Standard errors using a lag order of ", lags , ":\n", sep=""))
  mPrint = rbind(coefficients , "Standard errors" =  diag(mVarCovar), 
                 "T-Statistics" = coefficients/diag(mVarCovar), "P-values" = dt(coefficients/diag(mVarCovar) , 
                                                                                df=object@model$df.residual))
  print(round(mPrint,6 ) ,5 )
  cat("\n-------------------------------------------------------------\n")
  return("HACmatrix" = mVarCovar)
})

setGeneric("getForc", function(object, whichStep = 1)
  standardGeneric("getForc")  
)

setMethod("getForc", signature(object = "HARForecast"), function(object, whichStep = 1){
  vDates = object@data$forecastDates + whichStep-1
  vForc = xts(object@forecast[whichStep,], vDates)
  return(vForc)
})

setGeneric("forecastRes", function(object)
  standardGeneric("forecastRes")
)


setMethod("forecastRes", signature(object = "HARForecast"), function(object){
  vRes = getForc(object) - object@data$forecastComparison
  return(vRes)
})


setGeneric("qlike", function(object)
  standardGeneric("qlike")
)

setMethod("qlike", signature(object = "HARModel"), function(object){
  RM = object@model$model[,1] #extract RM from model
  FV  = object@model$fitted.values
  qLike = RM/FV - log(RM / FV) - 1
  return(qLike) 
}
)

setMethod("qlike", signature(object = "HARForecast"), function(object){
  RM = object@data$forecastComparison  #extract observed RM
  FV  = getForc(object)
  qLike = RM/FV - log(RM / FV) - 1
  return(qLike) 
}
)

###Methods that work with "lm" objects that I thought may be useful.
###Wrappers from the HARModel to "lm"
setMethod("logLik" , signature(object = "HARModel") , function(object, ...){
  out = logLik(object@model, ...)
  return(out)
})

setMethod("confint" , signature(object = "HARModel") , function(object, parm, level = 0.95, ...){
  out = confint(object@model, parm = parm, level = level, ...)
  return(out)
})

setMethod("residuals" , signature(object = "HARModel") , function(object, ...){
  out = residuals(object@model, ...)
  return(out)
})

setMethod("summary" , signature(object = "HARModel") , function(object, ...){
  out = summary(object@model,...)
  out$call =  as.name("lm(y ~ x)")
  return(out)
})

setMethod("fitted.values", signature(object = "HARModel"), function(object, ...){
  out = fitted.values(object@model, ...)
  return(out)
})
