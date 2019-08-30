\name{HAREstimate}
\alias{HAREstimate}
\title{
HAR estimation
}
\description{
HAR estimation
}
\usage{
HAREstimate(RM, BPV = NULL, RQ = NULL, periods = c(1,5,22),
            periodsJ = NULL, periodsRQ = NULL, type = "HAR",
            insanityFilter = TRUE, h = 1)
}
\arguments{
  \item{RM}{
A \code{numeric} containing a realized measure of the integrated volatility.
}
  \item{BPV}{
A \code{numeric} containing the estimate of the continuous part of the integrated volatility used for \code{HARJ} and \code{HARQ-J} types.
  }
  \item{RQ}{
A \code{numeric} containing  the realized quarticity used for \code{HARQ} and \code{HARQ-J} types.
  }
  \item{periods}{
A \code{numeric} denoting which lags should be used in the estimation, standard of \code{c(1,5,22)} is in line with Corsi(2009).
}
  \item{periodsJ}{
A \code{numeric} denoting which lags should be used in Jump estimation, if applicable.
  }
  \item{periodsRQ}{
A \code{numeric} denoting which lags should be used in Realized Quarticity estimation, if applicable.
  }
  \item{type}{
A \code{character} denoting which type of HAR model to estimate.
  }
  \item{insanityFilter}{
A \code{logical} denoting whether the insanity filter should be used for the fitted values of the estimation see Bollerslev, Patton & Quaedvlieg(2016) footnote 17. 
  }
  \item{h}{
A \code{integer} denoting the whether and how much to aggregate the realized variance estimator, if h = 5 the model is for the weekly volatility and if h = 22, the model is for the monthly volatility, the default of 1 designates no aggregation.
  }
}
\value{
A \code{\linkS4class{HARModel}} object
}
\details{
The estimates for the HARQ and HARQ-J models differ slightly from the results of BPQ (2016). This is due to a small difference in the demeaning approach for the realized quarticity. Here, the demeaning is done with mean(RQ) over all periods.
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
FitHAR = HAREstimate(RM = SP500rv, periods = c(1,5,22))

#extract the estimated coefficients:
coef(FitHAR)
#plot the fitted values
plot(FitHAR)

#calculate the Q-like loss-function:
mean(qlike(FitHAR))




#HAR-J:
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500bpv = SP500RM$BPV

#Estimate the HAR-J model:
FitHARJ = HAREstimate(RM = SP500rv, BPV = SP500bpv, 
                       periods = c(1,5,22), periodsJ = c(1,5,22),  type = "HARJ" )

#Calculate the Q-like loss-function:
mean(qlike(FitHARJ))




#HAR-Q of BPQ(2016) with weekly aggregation
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
#Estimate the HAR-Q model:
FitHARQ = HAREstimate(RM = SP500rv, RQ = SP500rq, periods = c(1,5,22),
                       periodsRQ = c(1,5,22),  type = "HARQ", h = 5)
#Show the model:
show(FitHARQ)

#Extract the coefficients:
HARQcoef = coef(FitHARQ)




#HARQ-J of BPQ(2016) with monthly aggregation
#load data
data("SP500RM")
SP500rv = SP500RM$RV
SP500rq = SP500RM$RQ
SP500bpv = SP500RM$BPV

#Estimate the HARQ-J model:
FitHARQJ = HAREstimate(RM = SP500rv, BPV = SP500bpv,
                        RQ = SP500rq, periods = c(1,5,22),
                        periodsJ = c(1), periodsRQ = c(1), 
                        type = "HARQ-J", h = 22)
#show the model:
show(FitHARQJ)



}