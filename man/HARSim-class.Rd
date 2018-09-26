\name{HARSim-class}
\alias{HARSim-class}
\alias{show,HARSim-method}
\alias{plot,HARSim,missing-method}
\alias{uncondmean,HARSim-method}
\alias{coef,HARSim-method}
\docType{class}
\title{
HARSim
}
\description{
Class for HARSim object
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
    \item{\code{Elapsed Time}:} Object of class \code{difftime} containing the time elapsed in seconds}
    }
}

\section{Methods}{
\itemize{
\item{\code{show}: Shows summary}
\item{\code{plot}: Plot the forecasted series and observed series as well as the residuals}
\item{\code{uncondmean}: Extracts the unconditional mean from the simulation}
\item{\code{coef}: Extracts the coefficients from the simulation}
}}
\author{Emil Sjoerup}
\keyword{classes}