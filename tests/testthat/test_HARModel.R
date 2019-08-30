# estimation tests
library(testthat)
suppressPackageStartupMessages(library(HARModel))
context("HARModel tests")

test_that("estimation tests", {
  data("SP500RM")
  SP500rv = SP500RM$RV
  SP500rq = SP500RM$RQ
  SP500bpv = SP500RM$BPV
  SP500tpq = SP500RM$TPQ
  FitHAR = HAREstimate(SP500rv, periods = c(1,5,22))
  
  FitTVHAR = HAREstimate(SP500rv, periods = c(1,5,22), type = "TV-HAR")
  
  FitCHAR = HAREstimate(SP500rv, BPV = SP500bpv, periods = c(1,5,22), type = "CHAR")
  
  FitCHARQ = HAREstimate(SP500rv, BPV = SP500bpv, RQ = SP500tpq, periods = c(1,5,22), periodsRQ = c(1), type = "CHARQ")
  
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
  
  
  expect_equal(round(coef(FitHAR), digits = 7), 
               c("beta0" = 0.1123142, "beta1" = 0.2273436,"beta5" = 0.4903494, "beta22" = 0.1863766))
  expect_equal(coef(FitCHARQ), c("beta0" = -0.006432962, "beta1" = 0.583432070, "beta5" = 0.418886906, 
                                 "beta22" = 0.113096324, "beta_q1" = -0.541018770))
  expect_equal(coef(FitFullHARQJ), c("beta0" =-0.01978921, "beta1" = 0.56971635,  "beta5"= 0.46751653, "beta22" = 0.02578378,
                                     "beta_q1" = -0.32592648, "beta_q5" = 0.04371038, "beta_q22" = -0.11871874, "beta_j1" = -0.12456935, 
                                     "beta_j5" = -1.58884086, "beta_j22" = 1.67095035))
  expect_equal(mean(qlike(FitTVHAR)), 0.133913)
  expect_equal(mean(qlike(FitHARQ)), 0.13577366)
  expect_equal(mean(qlike(FitFullHARQ)), 0.13801195)
  expect_equal(mean(residuals(FitHARQ)^2), 2.35696246)
  expect_equal(mean(residuals(FitFullHARQ)^2), 2.35455696)
  expect_equal(mean(residuals(FitCHARQ)^2)/mean(residuals(FitHAR)^2), 0.9368255)
  expect_equal(mean(residuals(FitHARJ)^2)/mean(residuals(FitHAR)^2), 0.96833979)
  
}
)


test_that("forecasting tests", {
  data("SP500RM")
  SP500rv = SP500RM$RV
  SP500rq = SP500RM$RQ
  SP500bpv = SP500RM$BPV
  SP500tpq = SP500RM$TPQ
  nAhead = 1
  nRoll  = 3096
  
  ###Rolling window forecasting:
  
  ForcHAR = HARForecast(SP500rv, periods = c(1,5,22) ,nRoll = nRoll, nAhead = nAhead, type = "HAR" )
  
  
  ForcHARJ = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                         periodsJ = c(1) ,nRoll = nRoll,
                         nAhead = nAhead, type = "HARJ" )
  
  
  ForcCHAR = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22), 
                         nRoll = nRoll, nAhead = nAhead,
                         type = "CHAR")
  
  
  ForcFULLHARQ = HARForecast(SP500rv, RQ= SP500rq, periods = c(1,5,22), 
                             periodsRQ = c(1,5,22), nRoll = nRoll, nAhead = nAhead,
                             type = "HARQ")
  
  ###Increasing window forecasting:
  
  
  ForcHARI = HARForecast(SP500rv, periods = c(1,5,22) ,nRoll = nRoll, 
                         nAhead = 1, type = "HAR", windowType = "increasing" )
  
  ForcHARQI = HARForecast(SP500rv, RQ = SP500rq, periods = c(1,5,22), 
                          periodsRQ = c(1), nRoll = nRoll, nAhead = nAhead,
                          type = "HARQ", windowType = "increasing")
  
  ForcCHARQI = HARForecast(SP500rv, BPV = SP500bpv, RQ = SP500tpq, periods = c(1,5,22), 
                          periodsRQ = c(1), nRoll = nRoll, nAhead = nAhead,
                          type = "CHARQ", windowType = "increasing")
  
  
  
  
  
  ###BPQ table for out of sample forecasts 
  #Rolling window:
  
  expect_equal(mean(forecastRes(ForcHARJ)^2)/mean(forecastRes(ForcHAR)^2), 0.91751679)
  expect_equal(mean(forecastRes(ForcFULLHARQ)^2)/mean(forecastRes(ForcHAR)^2), 0.792906786)
  expect_equal(mean(forecastRes(ForcCHAR)^2) / mean(forecastRes(ForcHAR)^2), 0.95830776)
  expect_equal(mean(qlike(ForcHARJ))/mean(qlike(ForcHAR)), 1.01151469)
  expect_equal(mean(qlike(ForcCHAR))/mean(qlike(ForcHAR)), 1.02005358)
  
  #Increasing window:

  expect_equal(mean(forecastRes(ForcHARQI)^2) / mean(forecastRes(ForcHARI)^2), 0.89439887)
  expect_equal(mean(qlike(ForcHARQI))/mean(qlike(ForcHARI)), 0.88091626)
  
})
