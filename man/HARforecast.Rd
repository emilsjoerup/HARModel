\name{HARforecast}
\alias{HARforecast}
\title{
HAR forecasting
}
\description{
Does rolling out of sample forecasting of a HAR model.
}
\usage{
HARforecast(vRealizedMeasure, vAuxData = NULL, vLags = c(1,5,22),
vJumpLags = NULL, iNRoll=10, iNAhead=10, type = "HAR", HARQargs = list(demean = T))
}
\arguments{
  \item{vRealizedMeasure}{
An \code{xts} object containing a realized measure of the integrated volatility.
}
\item{vAuxData}{
A \code{numeric} or an \code{xts} object containing the jump proportion of the realized measure, or the realized quarticity used for \code{HARQ} and \code{Full-HARQ} types.
  }
\item{vLags}{
A vector denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
\item{vJumpLags}{
A \code{numeric} denoting which lags should be used in Jump estimation, if applicable.
  }
\item{iNRoll}{
How many rolling forecasts should be performed.
}
\item{iNAhead}{
The length of each rolling forecast.
}
\item{type}{
A \code{character} denoting which type of HAR model to estimate.
  }
\item{HARQargs}{
A \code{list} denoting the extra arguments of a \code{type} \code{HARQ} or \code{Full-HARQ} model. So far only the boolean \code{demean} is implemented, which determines whether the data should be demeaned as is the case in Bollerslev, Patton & Quaedvlieg(2016)
  }
}
\value{
Returns a \code{HARforecast} object containing
\item{1}{A \code{HARmodel} object estimated from vRealizedMeasure less the rolling forecasts}
\item{2}{\code{Forecast}, a data frame containing the forecasts. The first row is one-step ahead, the second row is two-step ahead etc. The first column contains the first row, the second column the second roll etc.}
\item{3}{\code{Info}, a list containing information about the forecasting, namely the elapsed time, the rolls done and the horizon of the rolls}
\item{4}{\code{Data}, a list containing \code{ForecastDates}, the dates over the rolling forecasts, \code{Observations} which is the in-sample observations, and \code{ForecastComparison}, containing the observations that are out of sample (for the first roll)}
}
\details{
Only the \code{HAR} type model is complete, which means to model more than one period ahead the realized quarticity and the jump component must be modeled separately. This is done by AR(1) processes in both cases.
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
See Also \code{\link{HARestimate}}
}
\examples{
set.seed(123)
#Simulate a HAR process:
HARSim = HARsimulate(iLength = 1500,
vLags = c(1,5,22), vCoef = c(0.5, 0.36 , 0.28 , 0.28 ))
#Do 200 rolling forecasts each of length 200 on the simulated process:
HARForecast = HARforecast(vRealizedMeasure = HARSim@Simulation , 
vLags = c(1, 5, 22), iNRoll = 200 , iNAhead = 200)
}
\keyword{ HAR }
\keyword{ forecast }
