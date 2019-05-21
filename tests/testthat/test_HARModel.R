# estimation tests
library(testthat)
library(HARModel)
context("HARModel tests")

test_that("estimation tests", {
  
  data("SP500RM", package = "HARModel")
  SP500rv = SP500RM$RV
  SP500rq = SP500RM$RQ
  SP500bpv = SP500RM$BPV
  SP500tpq = SP500RM$TPQ
  FitHAR = HARestimate(SP500rv, periods = c(1,5,22))
  
  FitTVHAR = HARestimate(SP500rv, periods = c(1,5,22), type = "TV-HAR")
  
  FitCHAR = HARestimate(SP500rv, BPV = SP500bpv, periods = c(1,5,22), type = "CHAR")
  
  FitCHARQ = HARestimate(SP500rv, BPV = SP500bpv, RQ = SP500tpq, periods = c(1,5,22), periodsRQ = c(1), type = "CHARQ")
  
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
  
  
  expect_equal(round(coef(FitHAR), digits = 7), 
               c("beta0" = 0.1123142, "beta1" = 0.2273436,"beta5" = 0.4903494, "beta22" = 0.1863766))
  expect_equal(coef(FitCHARQ), c("beta0" = -0.006432962, "beta1" = 0.583432070, "beta5" = 0.418886906, 
                                 "beta22" = 0.113096324, "beta_q1" = -0.541018770))
  expect_equal(coef(FitFullHARQJ), c("beta0" =-0.01978921, "beta1" = 0.56971635,  "beta5"= 0.46751653, "beta22" = 0.02578378,
                                     "beta_q1" = -0.32592648, "beta_q5" = 0.04371038, "beta_q22" = -0.11871874, "beta_j1" = -0.12456935, 
                                     "beta_j5" = -1.58884086, "beta_j22" = 1.67095035))
  expect_equal(mean(qlike(FitTVHAR)), 0.133913)
  expect_equal(mean(qlike(FitHARQ)), 0.13577384)
  expect_equal(mean(qlike(FitFullHARQ)), 0.13801212)
  expect_equal(mean(residuals(FitHARQ)^2), 2.35696246)
  expect_equal(mean(residuals(FitFullHARQ)^2), 2.35455696)
  expect_equal(mean(residuals(FitCHARQ)^2)/mean(residuals(FitHAR)^2), 0.9368255)
  expect_equal(mean(residuals(FitHARJ)^2)/mean(residuals(FitHAR)^2), 0.96833979)
  
}
)


test_that("forecasting tests", {
  
  data("SP500RM", package = "HARModel")
  SP500rv = SP500RM$RV
  SP500rq = SP500RM$RQ
  SP500bpv = SP500RM$BPV
  SP500tpq = SP500RM$TPQ
  iNAhead = 1
  iNRoll  = 3096
  
  ###Rolling window forecasting:
  
  ForcHAR = HARforecast(SP500rv, periods = c(1,5,22) ,iNRoll = iNRoll, iNAhead = iNAhead, type = "HAR" )
  
  
  ForcHARJ = HARforecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                         periodsJ = c(1) ,iNRoll = iNRoll,
                         iNAhead = iNAhead, type = "HARJ" )
  
  
  ForcCHAR = HARforecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22), 
                         iNRoll = iNRoll, iNAhead = iNAhead,
                         type = "CHAR")
  
  
  ForcFULLHARQ = HARforecast(SP500rv, RQ= SP500rq, periods = c(1,5,22), 
                             periodsRQ = c(1,5,22), iNRoll = iNRoll, iNAhead = iNAhead,
                             type = "HARQ")
  
  ###Increasing window forecasting:
  
  
  ForcHARI = HARforecast(SP500rv, periods = c(1,5,22) ,iNRoll = iNRoll, 
                         iNAhead = 1, type = "HAR", windowType = "increasing" )
  
  ForcHARQI = HARforecast(SP500rv, RQ = SP500rq, periods = c(1,5,22), 
                          periodsRQ = c(1), iNRoll = iNRoll, iNAhead = iNAhead,
                          type = "HARQ", windowType = "increasing")
  
  ForcCHARQI = HARforecast(SP500rv, BPV = SP500bpv, RQ = SP500tpq, periods = c(1,5,22), 
                          periodsRQ = c(1), iNRoll = iNRoll, iNAhead = iNAhead,
                          type = "CHARQ", windowType = "increasing")
  
  
  
  
  
  ###BPQ table for out of sample forecasts 
  #Rolling window:
  
  expect_equal(mean(forecastres(ForcHARJ)^2)/mean(forecastres(ForcHAR)^2), 0.91755259)
  expect_equal(mean(forecastres(ForcFULLHARQ)^2)/mean(forecastres(ForcHAR)^2), 0.79287259)
  expect_equal(mean(forecastres(ForcCHAR)^2) / mean(forecastres(ForcHAR)^2), 0.95830776)
  expect_equal(mean(qlike(ForcHARJ))/mean(qlike(ForcHAR)), 1.01164543)
  expect_equal(mean(qlike(ForcCHAR))/mean(qlike(ForcHAR)), 1.02005358)
  
  #Increasing window:

  expect_equal(mean(forecastres(ForcHARQI)^2) / mean(forecastres(ForcHARI)^2), 0.89439887)
  expect_equal(mean(qlike(ForcHARQI))/mean(qlike(ForcHARI)), 0.88091626)
  
})
