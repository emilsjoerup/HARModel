\name{HARForecast}
\alias{HARForecast}
\title{
HAR forecasting
}
\description{
Rolling out of sample forecasting of a HAR model.
}
\usage{
HARForecast(RM, BPV= NULL, RQ = NULL , periods = c(1,5,22), 
            periodsJ = NULL, periodsRQ = NULL, nRoll=10 , nAhead=1 , type = "HAR",
            windowType = "rolling", insanityFilter = TRUE, h = 1)
}
\arguments{
  \item{RM}{
An \code{xts} object containing a realized measure of the integrated volatility.
}
  \item{BPV}{
A \code{numeric} containing the jump proportion of the realized measure used for \code{HARJ} and \code{HARQ-J} types.
  }
  \item{RQ}{
A \code{numeric} containing  the realized quarticity used for \code{HARQ} and \code{HARQ-J} types.
  }
\item{periods}{
A vector denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
\item{periodsJ}{
A \code{numeric} denoting which lags should be used in Jump estimation, if applicable.
  }
  \item{periodsRQ}{
A \code{numeric} denoting which lags should be used in Realized Quarticity estimation, if applicable.
  }
\item{nRoll}{
How many rolling forecasts should be performed.
}
\item{nAhead}{
The length of each rolling forecast.
}
\item{type}{
A \code{character} denoting which type of HAR model to estimate.
  }
\item{windowType}{
A \code{character} denoting which kind of window to use, either "rolling"/"fixed" or "increasing"/"expanding". 2-letter abbreviations can be used.
}
\item{insanityFilter}{
A \code{logical} denoting whether the insanity filter should be used for the forecasted values see Bollerslev, Patton & Quaedvlieg(2016) footnote 17. 
  }
  \item{h}{
A \code{integer} denoting the whether and how much to aggregate the realized variance estimator, if h = 5 the model is forecasting the weekly volatility and if h = 22, the model is forecasting the monthly volatility, the default of 1 designates no aggregation..
  }
}
\value{
A \code{\linkS4class{HARForecast}} object
}
\details{
Not all models in this package are 'complete', which means some models use AR(1) processes to forecast e.g. realized quarticity in order to construct more than one step ahead forecasts.

The maximumm lag of the continuous or quarticity data must be lower than the maximum of the realized measure lag vector, the other cases are not implemented.

The estimates for the \code{HARQ} and \code{HARQ-J} models differ slightly from the results of BPQ (2016). This is due to a small difference in the demeaning approach for the realized quarticity. Here, the demeaning is done with mean(RQ) over all periods.

If \code{h} is greater than 1, then \code{nAhead} must be one, as multi-period ahead forecasts have not been implemented.


}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196.\cr
Bollerslev, T., Patton, A., Quaedvlieg, R. 2016,  Exploiting the errors: A simple approach for improved volatility forecasting, \emph{Journal of Econometrics }, vol.192 , issue 1, 1-18. \cr
}
\author{
Emil Sjoerup
}
\seealso{
See Also \code{\link{HAREstimate}}
}
\examples{
#HAR of Corsi(2009)
#load data:
data("SP500RM")
SP500rv = SP500RM$RV

ForecastHAR = HARForecast(SP500rv, periods = c(1,5,22), nRoll =50,
                          nAhead = 50, type = "HAR")

#plot the forecasted series along with the actual realizations:
plot(ForecastHAR)

#Calculate the MSE:
mean(forecastRes(ForecastHAR)^2)

#Calculate the Q-like loss function:
mean(qlike(ForecastHAR))




#HARJ
#load data:
data("SP500RM")
SP500rv = SP500RM$RV
SP500bpv = SP500RM$BPV


ForecastHARJ = HARForecast(SP500rv, BPV = SP500bpv, periods = c(1,5,22),
                            periodsJ = c(1,5,22) ,nRoll = 50,
                            nAhead = 50, type = "HARJ")

#Show the model:
show(ForecastHARJ)

#Extract the forecasted series:
forc = getForc(ForecastHARJ)




#HARQ BPQ(2016)
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ

ForecastHARQ = HARForecast(SP500rv, RQ= SP500rq, periods = c(1,5,22), 
                            periodsRQ = c(1,5,22), nRoll = 50, nAhead = 50,
                            type = "HARQ")





#HARQ-J BPQ(2016) with weekly aggregation.
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV

ForecastHARQJ = HARForecast(SP500rv, RQ = SP500rq, BPV = SP500bpv,
                             periods = c(1,5,22), periodsJ = c(1,5,22), 
                             periodsRQ = c(1,5,22), nRoll = 50,
                             nAhead = 1, type = "HARQ-J", h = 5)

}
\keyword{ HAR }
\keyword{ forecast }
