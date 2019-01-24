\name{HARforecast}
\alias{HARforecast}
\title{
HAR forecasting
}
\description{
Does rolling out of sample forecasting of a HAR model.
}
\usage{
HARforecast(vRealizedMeasure, vJumpComponent= NULL, vAuxData = NULL , vLags = c(1,5,22), 
            vJumpLags = NULL, vAuxLags = NULL, iNRoll=10 , iNAhead=10 , type = "HAR",
            InsanityFilter = TRUE ,  HARQargs = list(demean = T))
}
\arguments{
  \item{vRealizedMeasure}{
An \code{xts} object containing a realized measure of the integrated volatility.
}
  \item{vJumpComponent}{
A \code{numeric} containing the jump proportion of the realized measure used for \code{HARJ} and \code{HARQ-J} types.
  }
  \item{vAuxData}{
A \code{numeric} containing  the realized quarticity used for \code{HARQ} and \code{HARQ-J} types.
  }
\item{vLags}{
A vector denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
\item{vJumpLags}{
A \code{numeric} denoting which lags should be used in Jump estimation, if applicable.
  }
  \item{vAuxLags}{
A \code{numeric} denoting which lags should be used in Realized Quarticity estimation, if applicable.
  }
\item{iNRoll}{
How many rolling forecasts should be performed.
}
\item{iNAhead}{
The length of each rolling forecast.
}
\item{type}{
A \code{character} denoting which type of HAR model to estimate.
  }
\item{InsanityFilter}{
A \code{logical} denoting whether the insanity filter should be used for the forecasted values see Bollerslev, Patton & Quaedvlieg(2016) footnote 17. 
  }
\item{HARQargs}{
A \code{list} denoting the extra arguments of a \code{type} \code{HARQ} model. So far only the boolean \code{demean} is implemented, which determines whether the data should be demeaned as is the case in Bollerslev, Patton & Quaedvlieg(2016)
  }
}
\value{
A \code{\linkS4class{HARForecast}} object
}
\details{
Only the \code{HAR} type model is complete, which means to model more than one period ahead the realized quarticity and the jump component must be modeled separately. This is done by AR(1) processes in both cases.

The maximumm lag of the jump or auxiliary data must be lower than the maximum of the realized measure lag vector, the other cases are not implemented yet.
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196.\cr
<<<<<<< HEAD
=======
Bollerslev, T., Patton, A., Quaedvlieg, R. 2016,  Exploiting the errors: A simple approach for improved volatility forecasting, \emph{Journal of Econometrics }, vol.192 , issue 1, 1-18. \cr
>>>>>>> Development
}
\author{
Emil Sjoerup
}
\seealso{
See Also \code{\link{HARestimate}}
}
\examples{
#HAR of Corsi(2009)
#load data:
data("SP500RM")
SP500rv = SP500RM$RV

ForecastHAR = HARforecast(SP500rv, vLags = c(1,5,22), iNRoll =50,
                          iNAhead = 50, type = "HAR")


#HARJ
#load data:
data("SP500RM")
SP500rv = SP500RM$RV
SP500bpv = SP500RM$BPV
vJumpComponent = SP500rv - SP500bpv
vJumpComponent = ifelse(vJumpComponent>=0, vJumpComponent, 0)


ForecastHARJ = HARforecast(SP500rv, vJumpComponent = vJumpComponent, vLags = c(1,5,22),
                            vJumpLags = c(1,5,22) ,iNRoll = 50,
                            iNAhead = 50, type = "HARJ")



#HARQ BPQ(2016)
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ

ForecastHARQ = HARforecast(SP500rv, vAuxData= SP500rq, vLags = c(1,5,22), 
                            vAuxLags = c(1,5,22), iNRoll = 50, iNAhead = 50,
                            type = "HARQ")




#HARQ-J BPQ(2016)
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV
vJumpComponent = SP500rv - SP500bpv
vJumpComponent = ifelse(vJumpComponent>=0, vJumpComponent, 0)

ForecastHARQJ = HARforecast(SP500rv, vAuxData = SP500rq, vJumpComponent = vJumpComponent,
                             vLags = c(1,5,22), vJumpLags = c(1,5,22), 
                             vAuxLags = c(1,5,22), iNRoll = 50,
                             iNAhead = 50, type = "HARQ-J")

}
\keyword{ HAR }
\keyword{ forecast }
