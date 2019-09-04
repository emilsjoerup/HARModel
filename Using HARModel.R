#Using HAR model Package
rm(list = ls())

library(HARModel)
data("SP500RM")

SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV

#######################Estimation#######################
HARFit = HAREstimate(SP500rv, periods = c(1,5,22))

CHARFit = HAREstimate(SP500rv, BPV = SP500bpv, periods = c(1,5,22), type = "CHAR")

CHARQFit = HAREstimate(SP500rv, BPV = SP500bpv, RQ = SP500rq, 
                       periods = c(1,5,22), periodsRQ = c(1), type = "CHARQ")

HARJFit = HAREstimate(RM = SP500rv, BPV = SP500bpv, 
                      periods = c(1,5,22), periodsJ = c(1),  type = "HARJ" )

HARQFit = HAREstimate(RM = SP500rv, RQ = SP500rq, 
                      periods = c(1,5,22 ), periodsRQ = c(1),  type = "HARQ")

FullHARQFit = HAREstimate(RM = SP500rv, RQ = SP500rq,
                          periods = c(1,5,22 ), periodsRQ = c(1,5,22), type = "HARQ")

HARQJFit = HAREstimate(RM = SP500rv, BPV = SP500bpv, RQ = SP500rq, 
                       periods = c(1,5,22), periodsJ = c(1,5,22), 
                       periodsRQ = c(1), type = "HARQ-J")

FullHARQJFit = HAREstimate(RM = SP500rv, BPV = SP500bpv, RQ = SP500rq, 
                           periods = c(1,5,22), periodsJ = c(1,5,22), 
                           periodsRQ = c(1,5,22), type = "HARQ-J")

plot(HARQFit)

summary(HARQFit) #Wrapper for the lm submodel

logLik(HARQFit) #Wrapper for the lm submodel

confint(HARQFit) #Wrapper for the lm submodel
#Q-like:
mean(qlike(HARFit))
mean(qlike(HARQFit))



#######################Forecasting#######################

HARForc = HARForecast(SP500rv, periods = c(1,5,22),
                      nRoll = 3096, nAhead = 50, type = "HAR" )

HARIForc = HARForecast(SP500rv, periods = c(1,5,22) ,
                       nRoll = 3096, nAhead = 1, type = "HAR", windowType = "increasing" )

HARJForc = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                       periodsJ = c(1) ,nRoll = 3096, 
                       nAhead = 50, type = "HARJ" )

HARJIForc = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                        periodsJ = c(1) ,nRoll = 3096,nAhead = 1,
                        type = "HARJ", windowType = "increasing" )

FULLHARJForc = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                           periodsJ = c(1,5,22) ,nRoll = 3096,nAhead = 1, type = "HARJ" )

HARQForc = HARForecast(SP500rv, RQ = SP500rq, periods = c(1,5,22),
                       periodsRQ = c(1), nRoll = 3096, nAhead = 1,type = "HARQ" )

FULLHARQForc = HARForecast(SP500rv, RQ= SP500rq, periods = c(1,5,22),
                           periodsRQ = c(1,5,22),
                           nRoll = 3096, nAhead = 1,type = "HARQ")

HARQJForc = HARForecast(SP500rv, RQ = SP500rq, BPV = SP500bpv,
                        periods = c(1,5,22), periodsJ = c(1), periodsRQ = c(1),
                        nRoll = 3096,nAhead = 1, type = "HARQ-J")

CHARForc = HARForecast(SP500rv, RQ = SP500rq, BPV = SP500bpv,
                       periods = c(1,5,22), periodsJ = c(1), periodsRQ = c(1),
                       nRoll = 3096,nAhead = 1, type = "CHAR")

CHARQForc = HARForecast(SP500rv, RQ = SP500rq, BPV = SP500bpv,
                        periods = c(1,5,22), periodsJ = c(1), periodsRQ = c(1),
                        nRoll = 3096,nAhead = 1, type = "CHARQ")

#Plot:
plot(HARForc)

#To get the forecasted series of a HARForecast object:
forecast = getForc(HARForc)

#To get the forecast residualse:
forcRes = forecastRes(HARForc)

##Select loss functions from updated table 4 of BPQ:

mean(forecastRes(HARJForc)^2)/mean(forecastRes(HARForc)^2)
mean(forecastRes(HARQForc)^2)/mean(forecastRes(HARForc)^2)
mean(forecastRes(CHARForc)^2)/mean(forecastRes(HARForc)^2)


mean(qlike(HARJForc))/mean(qlike(HARForc))
mean(qlike(HARQForc))/mean(qlike(HARForc))#Slightly different from BPQ due to how I do the averaging vs how they do it (Both schemes make sense!)
mean(qlike(CHARForc))/mean(qlike(HARForc))


#######################Simulation#######################
set.seed(123)
lSim = HARSimulate(len = 4.5e4,  periods = c(1,5,22),
                   coef = c(0.1, 0.38, 0.28, 0.28), errorTermSD = 0.01)

plot(lSim)

FitSim = HAREstimate(lSim@simulation, periods = c(1,5,22))

plot(FitSim)
