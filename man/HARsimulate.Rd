\name{HARsimulate}
\alias{HARsimulate}
\title{
HAR simulation
}
\description{
Simulates a HAR model. From using the AR representation of the HAR model.
}
\usage{
HARsimulate(iLength=1500, vLags = c(1, 5, 22), 
vCoef = c(0.01, 0.36 ,0.28 , 0.28), dSigma = 0.001, show = TRUE)
}
\arguments{
  \item{iLength}{
Integer length of the simulated process.
}
  \item{vLags}{
Vector of lags for constructing the model, standard is c(1,5,22) which is in line with Corsi(2009)
}
  \item{vCoef}{
Coefficient vector which will be used to simulate the process. 
}
  \item{dSigma}{
Standard deviation of the error term.
}

\item{show}{
Logical to determine whether the output should be shown when done.
}
}
\value{
Returns an S4 object of type \code{HARSim} which contains:
\item{Simulation }{The simulated process}
\item{Info }{Information about the simulation (The input parameters) and the elapsed time.}

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
HARSim = HARsimulate(iLength=10000, vLags = c(1, 5, 22),
vCoef = c(0.01, 0.36, 0.28, 0.28), dSigma = 0.001)

HARFit = HARestimate(HARSim@Simulation, vLags = c(1, 5, 22))

}
\keyword{ HAR }
\keyword{ simulation }
