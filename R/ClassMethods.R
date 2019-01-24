setClass("HARModel", representation(Model = "lm", Info = "list", Data = "list"))
setClass("HARForecast" , representation(Model = "HARModel" , Forecast = "matrix", Info = "list" , Data = "list"))
setClass("HARSim" , representation(Simulation = "numeric" , Info = "list"))

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
  print(round(coefficients , 6) , digits = 5)  
  cat("\n")
  cat("\n Elapsed Time:" , round(as.double(object@Info[["ElapsedTime"]] , units = "secs" ) , digits = 3)  , "seconds\n")
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("show" , signature(object = "HARForecast") , function(object){
  
  cat("\n First model estimated:")
  show(object@Model)
  
  cat("\n-------------------------------HARForecast-------------------------------\n")
  cat("\n Forecast specification:\n")
  cat("\n Length of rolls:" , dim(object@Forecast)[2],"\n")
  cat("\n Rolls performed:" , dim(object@Forecast)[1],"\n")
  
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
  print(round(coefficients , 6) , digits = 5)  
  cat("\n")
  cat("\n Elapsed Time:" , round(as.double(object@Info[["ElapsedTime"]] , units = "secs" ) , digits = 3)  , "seconds\n")
  cat("\n-------------------------------------------------------------------------\n")
})

setMethod("plot" , signature(x= "HARModel", y = "missing"), 
          function(x, legend.loc = "topright",
                   col = 2:1, lwd= 2, 
                   main = "Realized measure vs. Fitted values",
                   legend.names = c("Realized Measure" , "Fitted values"), ...){
  
  vY = x@Data$`RealizedMeasure`
  vFitted.Val = xts(x@Model$fitted.values , order.by = index(vY))
  p1 = plot(cbind(vFitted.Val, vY), main = main, col = col, ...)
  p1 = addLegend(legend.loc = legend.loc, legend.names = legend.names,
                   col = col, lwd = lwd, ...)
  p1
})

setMethod("plot" , signature(x = "HARForecast", y = "missing"), 
          function(x, legend.loc = "topright",
                   main = "Observed vs. forecasted",
                   legend.names = c("Realized Measure" , "Forecasted Values"),
                   col = 2:1, lwd= 2, ...){
  
  vForecastComp = x@Data$`ForecastComparison`
  vRollingForecastplot = xts(x@Forecast[1,], index(vForecastComp))
  p1 = plot(cbind(vRollingForecastplot,vForecastComp), col = col, main = main, ...)
  p1 = addLegend(legend.loc = legend.loc, legend.names = legend.names, col = col, lwd = lwd, ...)
  p1
})

setMethod("plot" , signature(x = "HARSim" , y = "missing"),
          function(x , length = "ALL" , ctrl = "start", main = "Simulated RV", ...){
  vY = xts(x@Simulation , order.by = as.Date(1:length(x@Simulation), origin = "1970/01/01"))
  if(length == "ALL"){
    p1 = plot(vY , main = main, ...)
    }
  else if (ctrl == "start" & is(length,"numeric")){
    p1 = plot(vY[1:length], main = main, ...)
  }else if(ctrl == "end" & is(length,"numeric")){
    p1 = plot(vY[(length(vY) - length):length(vY)], main = main, ... )
  }
  p1
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
  mPrint = rbind(coefficients , "Standard errors" =  diag(mVarCovar) , 
                 "T-Statistics" = coefficients/diag(mVarCovar), "P-values" = dt(coefficients/diag(mVarCovar) , 
                                                                                df=object@Model$df.residual))
  print(round(mPrint,6 ) ,5 )
  cat("\n-------------------------------------------------------------\n")
  return("HACmatrix" = mVarCovar)
})