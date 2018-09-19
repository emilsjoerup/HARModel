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
Only does 1 repetition.
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
#A Single repetition:
HARMonteCarlo = HARMonteCarlo(iLength = 1000, vLags = c(1, 5, 22),
vCoef = c(1, 0.36, 0.28, 0.28), iBurnin = 100, dSigma = 1,
iLagsPlusOne = length(vLags) + 1)

#An example using parallel computing:
library(parallel)
iLength = 1000
iBurnin = 50
vLags = c(1,5,22)
vCoef = c(0.5, 0.36 , 0.28 , 0.28 )
iB = 1000
mAPar = matrix(numeric(iB*length(vCoef)) , ncol =length(vCoef) )
set.seed(1)
iJ = detectCores()
Cluster = makeCluster(iJ)
vSeeds = sample(1:44684, iJ)
for (j in 1:iJ) {
  foo = vSeeds[j]
  clusterExport(Cluster[j] , c("foo",
  "HARDataCreationC" , "HARsimulate" , "HARMonteCarlo"))
  clusterEvalQ(Cluster[j], set.seed(foo))
}


mAPar = do.call(rbind, parLapply(cl = Cluster , 1:iB, 
function(b, iLength = iLength , iBurnin = iBurnin,
vLags = vLags, vCoef = vCoef){
HARMonteCarlo(iLength = iLength , iBurnin = iBurnin, vLags = vLags, 
vCoef = vCoef , iLagsPlusOne = 4)
} , iLength = iLength , iBurnin = iBurnin,
vLags = vLags, vCoef = vCoef))


stopCluster(Cluster)

#Plotting the results:

par(mfrow = c(2,2), ps = 20 , cex = 0.6 , cex.main=1)
hist(mAPar[,1] , main="constant" , xlab="",
breaks = 50 , freq=T , col = "light green")
abline( v= vCoef[1] , col="red")


hist(mAPar[,2] , main="beta1" , xlab="" ,
breaks = 50 , freq=T,col = "light green")
abline(v= vCoef[2] , col="red")


hist(mAPar[,3] , main="beta2" , xlab="" , 
breaks = 50 , freq=T,col = "light green")
abline(v= vCoef[3] , col="red")


hist(mAPar[,4] , main="beta3" , xlab="" ,
breaks = 50 , freq = T,col = "light green")
abline(v= vCoef[4] , col="red")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
