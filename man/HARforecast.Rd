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
\item{2}{ForecastMatrix , matrix containing the forecasts. The first row is one-step ahead, the second row is two-step ahead etc. The first column contains the first row, the second column the second roll etc.}
\item{3}{vObservations,  the data that is used in the HARmodel estimation of the in-sample estimation}
\item{4}{vForecastComp, an \code{xts} object containing the out of sample observations, which is useful for forecast accuracy assessment}
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
Empty until I find data to include in the package, sorry.
The usage of the function is rather simple, 
if you have problems, send me an e-mail. 
Will be updated.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
