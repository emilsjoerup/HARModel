\name{HARSimulate}
\alias{HARSimulate}
\title{
HAR simulation
}
\description{
Simulates a HAR model. From using the AR representation of the HAR model.
}
\usage{
HARSimulate(len=1500, periods = c(1, 5, 22), 
            coef = c(0.01, 0.36 ,0.28 , 0.28), errorTermSD = 0.001)
}
\arguments{
  \item{len}{
An \code{integer} determining the length of the simulated process.
}
  \item{periods}{
A \code{numeric} of lags for constructing the model, standard is c(1,5,22).
}
  \item{coef}{
A \code{numeric} of coefficients which will be used to simulate the process. 
}
  \item{errorTermSD}{
A \code{numeric} determining the standard deviation of the error term.
}

}
\value{
A \code{\linkS4class{HARSim}} object
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196.\cr
}
\author{
Emil Sjoerup
}
\examples{
set.seed(123)
#Simulate the process of size 10000
HARSim = HARSimulate(len = 10000, periods = c(1, 5, 22),
                     coef = c(0.01, 0.36, 0.28, 0.28), errorTermSD = 0.001)
HARFit = HAREstimate(HARSim@simulation, periods = c(1, 5, 22))

}
\keyword{ HAR }
\keyword{ simulation }
