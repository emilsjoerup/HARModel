\name{HARMonteCarlo}
\alias{HARMonteCarlo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Monte Carlo simulation of HAR model
}
\description{
Does simulation and fast estimation of a HAR process. Note that only coefficient estimates are returned.
}
\usage{
HARMonteCarlo(iLength = 1000, vLags = c(1, 5, 22),
vCoef = c(1, 0.36, 0.28, 0.28), iBurnin = 100, dSigma = 1,
iLagsPlusOne = length(vLags) + 1)
}

\arguments{
  \item{iLength}{
Length of the simulated process.
}
  \item{vLags}{
Vector of lags for constructing the model, standard is c(1,5,22) which is in line with Corsi(2009)
}
  \item{vCoef}{
Coefficient vector which will be used to simulate the process. 
}
  \item{iBurnin}{
Integer to determine the length of burnin in the simulation. Increasing this value beyond 100 does not change much, but lowering it CAN be bad.
}
  \item{dSigma}{
Standard deviation of the error term.
}
  \item{iLagsPlusOne}{
Integer to denote the amount of lags + 1. This Parameter is used in a non-user-callable function. Setting the function to 4 (in the case of 3 lags) significantly speeds up MonteCarlo simulations.
}
}
\details{
Does fast estimation but only returns coefficient estimates. This function is created with the intention of making a function that is easy to use with the \code{parallel} package.
}
\value{
Returns a vector containing the parameter estimates.
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196 .\cr
}
\author{
Emil Sjoerup
}
\note{
Using this function once makes little sense but an example is given none-theless.
This function was made with parallel computing in mind and thus an example utilizing the \code{parallel} package is provided.
}


\seealso{
See Also \code{\link{HARestimate}} & \code{\link{HARsimulate}}
}
\examples{
\dontrun{
#A Single repetition:
HARMonteCarlo = HARMonteCarlo(iLength = 1000, vLags = c(1, 5, 22),
vCoef = c(1, 0.36, 0.28, 0.28), iBurnin = 100, dSigma = 1,
iLagsPlusOne = length(vLags) + 1)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MonteCarlo }% use one of  RShowDoc("KEYWORDS")
\keyword{ HAR }% __ONLY ONE__ keyword per line
