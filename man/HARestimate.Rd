\name{HARestimate}
\alias{HARestimate}
\title{
HAR estimation
}
\description{
This function does HAR estimation in line with Corsi(2009). It is important to note that the maximum of the lag vector (22 in the standard case) observations are omitted because these are used to generate the data.
}
\usage{
HARestimate(vRealizedMeasure, vLags = c(1,5,22), 
iLagSE = 5 , show = TRUE)
}
\arguments{
  \item{vRealizedMeasure}{
An \code{xts} object containing a realized measure of the integrated volatility.
}
  \item{vLags}{
A vector denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
\item{iLagSE}{
An integer used to calculate the Newey-West standard errors.
}
\item{show}{
Logical to determine whether the output should be shown when done. Setting to FALSE makes Monte Carlo simulation and forecasting much faster.
}
}
\details{
Only \code{xts} objects are supported.
}
\value{
Returns a \code{HARModel} object which is an S4 class object containing:
\item{1}{\code{Model}, an \code{lm} object containing the estimation details. The \code{lm} object also contains the Newey-West variance covariance matrix, which is called \code{mVarCovar}. }
\item{2}{a \code{list} containing information used to estimate, namely the Newey-West lag order, the vector of lags , \code{vLags}, and \code{Dates}, which is a vector of the dates of the data (if the provided realized measure is an "xts" object).}
\item{3}{The List \code{Data} contains the realized measure provided.}

}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196 .\cr
}
\author{
Emil Sjoerup
}

\seealso{
See Also \code{\link{HARDataCreationC}}
}
\examples{
\dontrun{
set.seed(123)
#Simulate a HAR process:
HARSim = HARsimulate(iLength = 10000 , iBurnin = 1000, 
vLags = c(1,5,22), vCoef = c(0.5, 0.36 , 0.28 , 0.28 ))
#Estimate the HAR process:
HARModel = HARestimate(vRealizedMeasure = HARSim@Simulation,
vLags = c(1,5,22) , iLagSE = 5)
}
}
\note{
}
