\name{uncondmean}
\alias{uncondmean}

\title{
Unconditional mean
}
\description{
Extracts the unconditional mean of a HAR models and other objects
}
\usage{
uncondmean(object)
}

\arguments{
  \item{object}{
an object from the \code{HARModel} package
}
}
\value{\item{Unconditional Mean}{ a double containing the unconditional mean of the model or simulation}}
\author{Emil Sjoerup}
\note{
This is a method
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
# uncondmean can be used on both Simulation and
# estimated models
uncondmean(HARSim)
uncondmean(HARModel)
}
\keyword{HAR}
\keyword{Unconditional}
