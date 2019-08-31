#Using HAR model Package
rm(list = ls())

library(HARModel)
data("SP500RM")

SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV

#######################Estimation#######################
FitHAR = HAREstimate(as.vector(SP500rv), periods = c(1,5,22))

FitCHAR = HAREstimate(SP500rv, BPV = SP500bpv, periods = c(1,5,22), type = "CHAR")

FitCHARQ = HAREstimate(SP500rv, BPV = SP500bpv, RQ = SP500rq, periods = c(1,5,22), periodsRQ = c(1), type = "CHARQ")

FitHARJ = HAREstimate(RM = SP500rv, BPV = SP500bpv,
                       periods = c(1,5,22), periodsJ = c(1),  type = "HARJ" )


FitHARQ = HAREstimate(RM = SP500rv, RQ = SP500rq, periods = c(1,5,22 ), 
                       periodsRQ = c(1),  type = "HARQ")

FitFullHARQ = HAREstimate(RM = SP500rv, RQ = SP500rq, periods = c(1,5,22 ),
                       periodsRQ = c(1,5,22), type = "HARQ")

FitHARQJ = HAREstimate(RM = SP500rv, BPV = SP500bpv, 
                        RQ = SP500rq, periods = c(1,5,22),
                        periodsJ = c(1,5,22), periodsRQ = c(1), type = "HARQ-J")

FitFullHARQJ = HAREstimate(RM = SP500rv, BPV = SP500bpv, 
                            RQ = SP500rq, periods = c(1,5,22),
                            periodsJ = c(1,5,22), periodsRQ = c(1,5,22), type = "HARQ-J")

plot(FitCHARQ)

summary(FitHAR)

logLik(FitHAR)

confint(FitHAR)
#Q-like:
mean(qlike(FitHAR))
mean(qlike(FitHARQ))



#######################Forecasting#######################

ForcHAR = HARForecast(SP500rv, periods = c(1,5,22) ,nRoll = 3096, nAhead = 50, type = "HAR" )
ForcHARI = HARForecast(SP500rv, periods = c(1,5,22) ,nRoll = 3096, 
                       nAhead = 1, type = "HAR", windowType = "increasing" )


ForcHARJ = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1) ,nRoll = 3096,
                            nAhead = 50, type = "HARJ" )
ForcHARJI = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1) ,nRoll = 3096,
                            nAhead = 1, type = "HARJ", windowType = "increasing" )

ForcFULLHARJ = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1,5,22) ,nRoll = 3096,
                            nAhead = 1, type = "HARJ" )
ForcHARQ = HARForecast(as.vector(SP500rv), RQ = SP500rq, periods = c(1,5,22), 
                            periodsRQ = c(1), nRoll = 3096, nAhead = 1,
                            type = "HARQ" )

ForcFULLHARQ = HARForecast(SP500rv, RQ= SP500rq, periods = c(1,5,22), 
                            periodsRQ = c(1,5,22), nRoll = 3096, nAhead = 1,
                            type = "HARQ" )

ForcHARQJ = HARForecast(as.numeric(SP500rv), RQ = SP500rq, BPV = SP500bpv,
                             periods = c(1,5,22), periodsJ = c(1), 
                             periodsRQ = c(1), nRoll = 3096,
                             nAhead = 1, type = "HARQ-J" )

ForcCHAR = HARForecast(as.numeric(SP500rv), RQ = SP500rq, BPV = SP500bpv,
                             periods = c(1,5,22), periodsJ = c(1), 
                             periodsRQ = c(1), nRoll = 3096,
                             nAhead = 1, type = "CHAR" )

ForcCHARQ = HARForecast(as.numeric(SP500rv), RQ = SP500rq, BPV = SP500bpv,
                            periods = c(1,5,22), periodsJ = c(1), 
                            periodsRQ = c(1), nRoll = 3096,
                            nAhead = 1, type = "CHARQ" )

##Select loss functions from updated table 4 of BPQ:

mean(forecastRes(ForcHAR)^2)/mean(forecastRes(ForcHAR)^2)
mean(forecastRes(ForcHARJ)^2)/mean(forecastRes(ForcHAR)^2)
mean(forecastRes(ForcCHAR)^2)/mean(forecastRes(ForcHAR)^2)


mean(qlike(ForcHAR))/mean(qlike(ForcHAR))
mean(qlike(ForcHARJ))/mean(qlike(ForcHAR))
mean(qlike(ForcCHAR))/mean(qlike(ForcHAR))


#######################Simulation#######################

set.seed(123)
lSim = HARSimulate(len = 4.5e4,  periods = c(1,5,22),
                   coef = c(0.1, 0.38, 0.28, 0.28), errorTermSD = 0.01)

plot(lSim)

FitSim = HAREstimate(lSim@Simulation, periods = c(1,5,22))

plot(FitSim)
