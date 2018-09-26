\name{SandwichNeweyWest}
\alias{SandwichNeweyWest}

\title{
HAC Variance Covariance Matrix For \code{HARModel} Objects
}
\description{
A method for extracting the Newey-West variance-covariance matrix for a \code{HARModel} object
}
\usage{
SandwichNeweyWest(object, lags)
}

\arguments{
  \item{object}{ A \code{HARModel} object
}
  \item{lags}{An integer to denote the allowed autocorrelation.
}
}
\details{
Usual values for \code{lags} are 5 and \eqn{T^0.25} rounding to the nearest integer.
This function is merely a method to implement code from the \code{sandwich} package.
}
\value{
\item{HACMatrix}{The Variance-Covariance matrix}
}
\references{
\code{sandwich} package. https://CRAN.R-project.org/package=sandwich
}
\examples{
set.seed(123)
#Simulate a HAR process:
HARSim = HARsimulate(iLength = 10000, 
vLags = c(1,5,22), vCoef = c(0., 0.36 , 0.28 , 0.28 ) 
, dSigma = 0.001)
#Estimate the HAR process:
HARModel = HARestimate(vRealizedMeasure = HARSim@Simulation,
vLags = c(1,5,22))

SandwichNeweyWest(HARModel , lags = 5)
}

