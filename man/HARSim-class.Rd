\name{HARSim-class}
\alias{HARSim-class}
\docType{class}

\title{
Class for HAR simulation
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it}
\section{Slots}{
  \describe{
    \item{\code{Simulation}:}{Object of class \code{numeric} containing the simulated series}
    \item{\code{Info}:}{Object of class \code{list} cointaining:}
    \itemize{
    \item{\code{Length}:} Object of class \code{numeric}  containing the length of the simulated series
    \item{\code{Lags}:} Object of class \code{numeric} containing the lag-vector used for simulation
    \item{\code{Coefficients}:} Object of class \code{numeric} containing the coefficients used for simulation
    \item{\code{ErrorTermSD}:} Object of class \code{numeric} containing the standard error of the error term
    \item{\code{BurninLength}:} Object of class \code{numeric} containing the length of the burnin period in the simulation
    \item{\code{Elapsed Time}:} Object of class \code{difftime} containing the time elapsed in seconds}
    \item{\code{Data}:}{Object of class \code{list} containing:}
    }
}

\section{Methods}{
\itemize{
\item{\code{show}: Shows summary}
\item{\code{plot}: Plot the forecasted series and observed series as well as the residuals}
}}
\author{Emil Sjoerup}
\keyword{classes}