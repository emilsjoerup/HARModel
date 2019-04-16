\name{HARForecast-class}
\alias{HARForecast-class}
\alias{show,HARForecast-method}
\alias{plot,HARForecast,missing-method}
\alias{uncmean,HARForecast-method}
\alias{coef,HARForecast-method}
\alias{qlike,HARForecast-method}
\alias{forecastres,HARForecast-method}
\alias{forecastres,ANY-method}
\alias{forecastres}
\docType{class}
\title{
HARForecast
}
\description{
Class for HARForecast object
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it}
\section{Slots}{
  \describe{
    \item{\code{Model}:}{Object of class \code{HARModel}. see \link{HARModel} }
    \item{\code{Forecast}:}{Object of class \code{matrix} containing the forecasted series}
    \item{\code{Info}:}{Object of class \code{list} cointaining:}
    \itemize{\item{\code{Elapsed Time}:} Object of class \code{difftime} containing the time elapsed in seconds
    \item{\code{Rolls}:} \code{Integer}  containing the amount of rolls done in the forecasting routine
    \item{\code{Horizon}:} \code{Integer} containing the length of the horizon used for forecasting during each of the rolls}
    \item{\code{Data}:}{Object of class \code{list} containing:}
    \itemize{\item{\code{Dates}:}Object of type \code{Integer} or \code{Date} containing the indices of the forecasted series either in integer or date format
    \item{\code{Observations}:}Object of type \code{numeric} or \code{xts} containing the in-sample observations
    \item{\code{ForecastComparison}:} Object of type \code{numeric} or \code{xts} containing the observations kept out of sample for the first roll
    }
    }
}

\section{Methods}{
\describe{
\item{\code{show}:}{\code{signature(object = "HARForecast")}: Shows summary}
\item{\code{plot}:}{\code{signature(x = "HARForecast", y = "missing")}: Plot the out of sample observed series with the forecasts overlayed}
\item{\code{uncmean}:}{\code{signature(object = "HARForecast")}: Extracts the unconditional mean from the Model}
\item{\code{coef}:}{\code{signature(object = "HARForecast")}: Extracts the coefficients from the first estimated Model}
\item{\code{qlike}:}{\code{signature(object = "HARForecast")}: Calculate the out of sample 'qlike' loss function for a HARForecast object}
\item{\code{forecastres}:}{\code{signature(object = "HARForecast")}: Retrieve the forecast residuals from HARForecast object}
}}
\author{Emil Sjoerup}
\keyword{classes}