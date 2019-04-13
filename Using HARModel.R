#Using HAR model Package
rm(list = ls())

library(HARModel)
data("SP500RM")

SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV

#######################Estimation#######################
FitHAR = HARestimate(as.vector(SP500rv), periods = c(1,5,22))

FitCHAR = HARestimate(SP500rv, BPV = SP500bpv, periods = c(1,5,22), type = "CHAR")

FitCHARQ = HARestimate(SP500rv, BPV = SP500bpv, RQ = SP500rq, periods = c(1,5,22), periodsRQ = c(1), type = "CHARQ")

FitHARJ = HARestimate(RM = SP500rv, BPV = SP500bpv,
                       periods = c(1,5,22), periodsJ = c(1),  type = "HARJ" )


FitHARQ = HARestimate(RM = SP500rv, RQ = SP500rq, periods = c(1,5,22 ), 
                       periodsRQ = c(1),  type = "HARQ")

FitFullHARQ = HARestimate(RM = SP500rv, RQ = SP500rq, periods = c(1,5,22 ),
                       periodsRQ = c(1,5,22), type = "HARQ")

FitHARQJ = HARestimate(RM = SP500rv, BPV = SP500bpv, 
                        RQ = SP500rq, periods = c(1,5,22),
                        periodsJ = c(1,5,22), periodsRQ = c(1), type = "HARQ-J")

FitFullHARQJ = HARestimate(RM = SP500rv, BPV = SP500bpv, 
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

ForcHAR = HARforecast(SP500rv, periods = c(1,5,22) ,iNRoll = 3096, iNAhead = 50, type = "HAR" )
ForcHARI = HARforecast(SP500rv, periods = c(1,5,22) ,iNRoll = 3096, 
                       iNAhead = 1, type = "HAR", windowType = "increasing" )


ForcHARJ = HARforecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1) ,iNRoll = 3096,
                            iNAhead = 50, type = "HARJ" )
ForcHARJI = HARforecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1) ,iNRoll = 3096,
                            iNAhead = 1, type = "HARJ", windowType = "increasing" )

ForcFULLHARJ = HARforecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1,5,22) ,iNRoll = 3096,
                            iNAhead = 1, type = "HARJ" )
ForcHARQ = HARforecast(as.vector(SP500rv), RQ = SP500rq, periods = c(1,5,22), 
                            periodsRQ = c(1), iNRoll = 3096, iNAhead = 1,
                            type = "HARQ" )

ForcFULLHARQ = HARforecast(SP500rv, RQ= SP500rq, periods = c(1,5,22), 
                            periodsRQ = c(1,5,22), iNRoll = 3096, iNAhead = 1,
                            type = "HARQ" )

ForcHARQJ = HARforecast(as.numeric(SP500rv), RQ = SP500rq, BPV = SP500bpv,
                             periods = c(1,5,22), periodsJ = c(1), 
                             periodsRQ = c(1), iNRoll = 3096,
                             iNAhead = 1, type = "HARQ-J" )

ForcCHAR = HARforecast(as.numeric(SP500rv), RQ = SP500rq, BPV = SP500bpv,
                             periods = c(1,5,22), periodsJ = c(1), 
                             periodsRQ = c(1), iNRoll = 3096,
                             iNAhead = 1, type = "CHAR" )

ForcCHARQ = HARforecast(as.numeric(SP500rv), RQ = SP500rq, BPV = SP500bpv,
                            periods = c(1,5,22), periodsJ = c(1), 
                            periodsRQ = c(1), iNRoll = 3096,
                            iNAhead = 1, type = "CHARQ" )

##Select loss functions from updated table 4 of BPQ:

mean(forecastres(ForcHAR)^2)/mean(forecastres(ForcHAR)^2)
mean(forecastres(ForcHARJ)^2)/mean(forecastres(ForcHAR)^2)
mean(forecastres(ForcCHAR)^2)/mean(forecastres(ForcHAR)^2)


mean(qlike(ForcHAR))/mean(qlike(ForcHAR))
mean(qlike(ForcHARJ))/mean(qlike(ForcHAR))
mean(qlike(ForcCHAR))/mean(qlike(ForcHAR))


#######################Simulation#######################

set.seed(123)
lSim = HARsimulate(iLength = 4.5e4,  periods = c(1,5,22),
                   coef = c(0.1, 0.38, 0.28, 0.28), Sigma = 0.01)

plot(lSim)

FitSim = HARestimate(lSim@Simulation, periods = c(1,5,22))

plot(FitSim)
