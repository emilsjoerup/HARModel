\name{HARForecast-class}
\alias{HARForecast-class}
\docType{class}

\title{
Class for HAR forecast
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it}
\section{Slots}{
  \describe{
    \item{\code{Model}:}{Object of class \code{HARModel}. see \link{HARModel} }
    \item{\code{Forecast}:}{Object of class \code{data.frame} containing the forecasted series}
    \item{\code{Info}:}{Object of class \code{list} cointaining:}
    \itemize{\item{\code{Elapsed Time}:} Object of class \code{difftime} containing the time elapsed in seconds
    \item{\code{Rolls}:} \code{double} Integer containing the amount of rolls done in the forecasting routine
    \item{\code{Horizon}:} \code{double} Integer containing the length of the horizon used for forecasting during each of the rolls}
    \item{\code{Data}:}{Object of class \code{list} containing:}
    \itemize{\item{\code{Dates}:}Object of type \code{Integer} or \code{Date} containing the indices of the forecasted series either in integer or date format
    \item{\code{Observations}:}Object of type \code{double} or \code{xts} containing the in-sample observations
    \item{\code{ForecastComparison}:} Object of type \code{double} or \code{xts} containing the observations kept out of sample
    }
    }
}

\section{Methods}{
\itemize{
\item{\code{show}: Shows summary}
\item{\code{plot}: Plot the forecasted series and observed series as well as the residuals}
}}
\author{Emil Sjoerup}
\keyword{classes}