\name{HARDataCreationC}
\alias{HARDataCreationC}
\title{
Data creation for HAR estimation and forecasting.
}
\description{
Constructs a matrix containing the realized measure and the lagged moving averages according to the supplied lag-vector.
}
\usage{
HARDataCreationC(vRealizedmeasure, vLags)
}

\arguments{
  \item{vRealizedmeasure}{
A vector containing a realized measure of the integrated volatility.
}
  \item{vLags}{
A vector denoting which lags should be used in the creation \code{c(1,5,22)} is in line with Corsi(2009).
}
}
\details{
The first \code{max(vLags)} rows contain NA's. Written in C++ for speed.
}
\value{
A matrix containing the data used for HAR estimation or forecasting.
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196 .\cr
}
\author{
Emil Sjoerup
}
\seealso{
See Also \code{\link{HARestimate}}
}
\examples{
\dontrun{
## Construct a non-negative vector of size 1000:
vData = rnorm(1000)^2
## Construct a lag-vector
vLags = c(1,5,22)
mExample = HARDataCreationC(vData , vLags)
## The dimensions of mExample is (978 , 4)
}
}

