\name{HARestimate}
\alias{HARestimate}
\title{
HAR estimation
}
\description{
This function does HAR estimation in line with Corsi(2009)..
}
\usage{
HARestimate(vRealizedMeasure, vJumpComponent = NULL, vAuxData = NULL, vLags = c(1,5,22),
            vJumpLags = NULL, vAuxLags = NULL, type = "HAR",
            InsanityFilter = TRUE, HARQargs = list(demean = TRUE ),show=TRUE )
}
\arguments{
  \item{vRealizedMeasure}{
A \code{numeric} containing a realized measure of the integrated volatility.
}
  \item{vJumpComponent}{
A \code{numeric} containing the jump proportion of the realized measure used for \code{HARJ} and \code{HARQ-J} types.
  }
  \item{vAuxData}{
A \code{numeric} containing  the realized quarticity used for \code{HARQ} and \code{HARQ-J} types.
  }
  \item{vLags}{
A \code{numeric} denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
  \item{vJumpLags}{
A \code{numeric} denoting which lags should be used in Jump estimation, if applicable.
  }
  \item{vAuxLags}{
A \code{numeric} denoting which lags should be used in Realized Quarticity estimation, if applicable.
  }
  \item{type}{
A \code{character} denoting which type of HAR model to estimate.
  }
  \item{InsanityFilter}{
A \code{logical} denoting whether the insanity filter should be used for the fitted values of the estimation see Bollerslev, Patton & Quaedvlieg(2016) footnote 17. 
  }
  \item{HARQargs}{
A \code{list} denoting the extra arguments of a \code{type} \code{HARQ} or \code{Full-HARQ} model. So far only the boolean \code{demean} is implemented, which determines whether the data should be demeaned as is the case in Bollerslev, Patton & Quaedvlieg(2016).
  }
\item{show}{
Logical to determine whether the output should be shown when done.
}
}
\value{
A \code{\linkS4class{HARModel}} object
}
\references{
Corsi, F. 2009, A Simple Approximate Long-Memory Model
of Realized Volatility, \emph{Journal of Financial Econometrics}, 174--196.\cr
Bollerslev, T., Patton, A., Quaedvlieg, R. 2016,  Exploiting the errors: A simple approach for improved volatility forecasting, \emph{Journal of Econometrics }, vol.192, issue 1, 1-18. \cr
}
\author{
Emil Sjoerup
}
\examples{
#Vanilla HAR from Corsi(2009)
#load data
data("SP500RM")
SP500rv = SP500RM$RV
#Estimate the HAR model:
FitHAR = HARestimate(vRealizedMeasure = SP500rv, vLags = c(1,5,22))

#HAR-J:
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500bpv = SP500RM$BPV

vJumpComponent = SP500rv - SP500bpv
vJumpComponent = ifelse(vJumpComponent>=0, vJumpComponent, 0)
#Estimate the HAR-J model:
FitHARJ = HARestimate(vRealizedMeasure = SP500rv, vJumpComponent = vJumpComponent, 
                       vLags = c(1,5,22), vJumpLags = c(1,5,22),  type = "HARJ" )

#HAR-Q of BPQ(2016)
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
#Estimate the HAR-Q model:
FitHARQ = HARestimate(vRealizedMeasure = SP500rv, vAuxData = SP500rq, vLags = c(1,5,22),
                       vAuxLags = c(1,5,22),  type = "HARQ", HARQargs = list(demean =TRUE))

#HARQ-J of BPQ(2016)
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV

vJumpComponent = SP500rv - SP500bpv
vJumpComponent = ifelse(vJumpComponent>=0, vJumpComponent, 0)
#Estimate the HARQ-J model:
FitHARQJ = HARestimate(vRealizedMeasure = SP500rv, vJumpComponent = vJumpComponent,
                        vAuxData = SP500rq, vLags = c(1,5,22),
                        vJumpLags = c(1,5,22), vAuxLags = c(1), 
                        type = "HARQ-J", HARQargs = list(demean = TRUE))


}