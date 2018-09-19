\name{HARforecast}
\alias{HARforecast}
\title{
HAR forecasting
}
\description{
Does rolling out of sample forecasting of a HAR model.
}
\usage{
HARforecast(vRealizedMeasure, vLags, iNRoll = 10, iNAhead = 10)
}
\arguments{
  \item{vRealizedMeasure}{
An \code{xts} object containing a realized measure of the integrated volatility.
}
  \item{vLags}{
A vector denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
  \item{iNRoll}{
How many rolling forecasts should be performed.
}
  \item{iNAhead}{
The length of each rolling forecast.
}
}
\value{
Returns a \code{HARforecast}, which is a subclass of \code{HARmodel}, object containing
\item{1}{A \code{HARmodel} object estimated from vRealizedMeasure less the rolling forecasts}
\item{2}{Forecast data frame containing the forecasts. The first row is one-step ahead, the second row is two-step ahead etc. The first column contains the first row, the second column the second roll etc.}
\item{3}{Info, a list containing information about the forecasting, namely the elapsed time, the rolls done and the horizon of the rolls}
\item{4}{Data, a list containing \code{ForecastDates}, the dates over the rolling forecasts, \code{Observations} which is the in-sample observations, and \code{ForecastComparison}, containing the observations that are out of sample (for the first roll)}
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196 .\cr
}
\author{
Emil Sjoerup
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See Also \code{\link{HARestimate}}
}
\examples{
\dontrun{
set.seed(123)
#Simulate a HAR process:
HARSim = HARsimulate(iLength = 10000 , iBurnin = 1000,
vLags = c(1,5,22), vCoef = c(0.5, 0.36 , 0.28 , 0.28 ))
#Do 200 rolling forecasts each of length 200 on the simulated process:
HARForecast = HARforecast(vRealizedMeasure = HARSim@Simulation , 
vLags = c(1, 5, 22), iNRoll = 200 , iNAhead = 200)
#Plot the forecast:
plot(HARForecast)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ HAR }% use one of  RShowDoc("KEYWORDS")
\keyword{ forecasting }% __ONLY ONE__ keyword per line
