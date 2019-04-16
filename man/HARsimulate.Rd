\name{HARsimulate}
\alias{HARsimulate}
\title{
HAR simulation
}
\description{
Simulates a HAR model. From using the AR representation of the HAR model.
}
\usage{
HARsimulate(iLength=1500, periods = c(1, 5, 22), 
coef = c(0.01, 0.36 ,0.28 , 0.28), Sigma = 0.001)
}
\arguments{
  \item{iLength}{
Integer length of the simulated process.
}
  \item{periods}{
Vector of lags for constructing the model, standard is c(1,5,22) which is in line with Corsi(2009)
}
  \item{coef}{
Coefficient vector which will be used to simulate the process. 
}
  \item{Sigma}{
Standard deviation of the error term.
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
HARSim = HARsimulate(iLength=10000, periods = c(1, 5, 22),
coef = c(0.01, 0.36, 0.28, 0.28), Sigma = 0.001)
HARFit = HARestimate(HARSim@Simulation, periods = c(1, 5, 22))

}
\keyword{ HAR }
\keyword{ simulation }
