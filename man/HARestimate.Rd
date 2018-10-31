\name{HARestimate}
\alias{HARestimate}
\title{
HAR estimation
}
\description{
This function does HAR estimation in line with Corsi(2009). It is important to note that the maximum of the lag vector (22 in the standard case) observations are omitted because these are used to generate the data.
}
\usage{
HARestimate(vRealizedMeasure, vAuxData = NULL ,  vLags = c(1,5,22), 
vJumpLags = NULL  , type = "HAR" , HARQargs = list(demean = T) ,show=TRUE )
}
\arguments{
  \item{vRealizedMeasure}{
A \code{numeric} or an \code{xts} object containing a realized measure of the integrated volatility.
}
  \item{vAuxData}{
A \code{numeric} or an \code{xts} object containing the jump proportion of the realized measure, or the realized quarticity used for \code{HARQ} and \code{Full-HARQ} types.
  }
  \item{vLags}{
A \code{numeric} denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
  \item{vJumpLags}{
A \code{numeric} denoting which lags should be used in Jump estimation, if applicable.
  }
  \item{type}{
A \code{character} denoting which type of HAR model to estimate.
  }
  \item{HARQargs}{
A \code{list} denoting the extra arguments of a \code{type} \code{HARQ} or \code{Full-HARQ} model. So far only the boolean \code{demean} is implemented, which determines whether the data should be demeaned as is the case in Bollerslev, Patton & Quaedvlieg(2016)
  }
\item{show}{
Logical to determine whether the output should be shown when done.
}
}
\value{
Returns a \code{HARModel} object which is an S4 class object containing:
\item{1}{\code{Model}, an \code{lm} object containing the estimation details.}
\item{2}{\code{Info}, a \code{list} containing information used, the vector of lags , \code{vLags}, \code{Dates}, which is a vector of the dates of the data (if the provided realized measure is an "xts" object), and \code{ElapsedTime} , the time elapsed in seconds. }
\item{3}{Data, a \code{list} containing: \code{Realized Measure} contains the realized measure provided.
                                         Depending on the type of model \code{Jump Component} or \code{Realized Quarticiy} will be returned too, if applicable. }
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196.\cr
Bollerslev, T., Patton, A., Quaedvlieg, R. 2016,  Exploiting the errors: A simple approach for improved volatility forecasting, \emph{Journal of Econometrics }, vol.192 , issue 1, 1-18. \cr
}
\author{
Emil Sjoerup
}

%\seealso{
%}
\examples{
set.seed(123)
#Simulate a HAR process:
HARSim = HARsimulate(iLength = 1500, 
vLags = c(1,5,22), vCoef = c(0.5, 0.36 , 0.28 , 0.28 ))
#Estimate the HAR process:
HARModel = HARestimate(vRealizedMeasure = HARSim@Simulation,
vLags = c(1,5,22))
}