\name{HARestimate}
\alias{HARestimate}
\title{
HAR estimation
}
\description{
This function does HAR estimation in line with Corsi(2009). It is important to note that the maximum of the lag vector (22 in the standard case) observations are omitted because these are used to generate the data.
}
\usage{
HARestimate(vRealizedMeasure, vLags, bStandardErrors = T, iLagSE = 5)
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
}
\details{
Only \code{xts} objects are supported.
}
\value{
Returns a \code{HARmodel}, which is a subclass of \code{lm}, object containing
\item{1}{An \code{lm} object with a rewritten terms section for nice formatting}
\item{2}{mVarCovar, a the variance-covariance matrix with the Newey-West standard errors}
\item{3}{dates, a \code{Date} object containing the dates of the data used for estimation}

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
Empty until I find data to include in the package, sorry.
The usage of the function is rather simple, 
if you have problems, send me an e-mail. 
Will be updated.
}
}
\note{
If \code{bStandardErrors} is \code{FALSE} , the class of the returned object is \code{lm}
}
