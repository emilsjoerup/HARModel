\name{HARSim-class}
\alias{HARSim-class}
\alias{show,HARSim-method}
\alias{plot,HARSim,missing-method}
\alias{uncmean,HARSim-method}
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
    \item{\code{simulation}:}{Object of class \code{numeric} containing the simulated series}
    \item{\code{info}:}{Object of class \code{list} cointaining:}
    \itemize{
    \item{\code{len}:} Object of class \code{numeric}  containing the length of the simulated series
    \item{\code{periods}:} Object of class \code{numeric} containing the lag-vector used for simulation
    \item{\code{coefficients}:} Object of class \code{numeric} containing the coefficients used for simulation
    \item{\code{errorTermSD}:} Object of class \code{numeric} containing the standard error of the error term
    \item{\code{elapsedTime}:} Object of class \code{difftime} containing the time elapsed in seconds}
    }
}

\section{Methods}{
\describe{
\item{\code{show}:}{\code{signature(object = "HARSim")}: Shows summary}
\item{\code{plot}:}{\code{signature(x = "HARSim" , y = "missing")}: Plot the forecasted series and observed series as well as the residuals}
\item{\code{uncmean}:}{\code{signature(object = "HARSim")}: Extracts the unconditional mean from the simulation}
\item{\code{coef}:}{\code{signature(object = "HARSim")}: Extracts the coefficients from the simulation}
}}
\author{Emil Sjoerup}
\keyword{classes}