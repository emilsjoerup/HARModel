\name{HARModel-class}
\alias{HARModel-class}
\alias{show,HARModel-method}
\alias{plot,HARModel,missing-method}
\alias{uncmean,HARModel-method}
\alias{uncmean,ANY-method}
\alias{uncmean}
\alias{coef,HARModel-method}
\alias{SandwichNeweyWest,HARModel-method}
\alias{SandwichNeweyWest,ANY-method}
\alias{SandwichNeweyWest}
\alias{qlike,HARModel-method}
\alias{qlike,ANY-method}
\alias{qlike}
\alias{logLik,HARModel-method}
\alias{confint,HARModel-method}
\alias{residuals,HARModel-method}
\alias{summary,HARModel-method}




\docType{class}
\title{
HARModel
}
\description{
Class for HARModel objects
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
  \describe{
    \item{\code{Model}:}{Object of class \code{lm}. Contains the linear model fitted.}
    \item{\code{Info}:}{Object of class \code{list} cointaining:}
    \itemize{\item{\code{Lags}:} \code{numeric} containing the lags used to create the model
    \item{\code{Dates}:} \code{Date} object containing the dates for which the estimation was done, only applicable if the Model was estimated using an "xts" object.}
    }
}
\section{Methods}{
\describe{
\item{\code{show}:}{\code{signature(object = "HARModel")} Shows summary}
\item{\code{plot}:}{\code{signature(x = "HARModel", y = "missing")}: Plots the observed values with fitted values overlayed}
\item{\code{uncmean}:}{\code{signature(object = "HARModel")}: Extracts the unconditional mean from the Model}
\item{\code{coef}:}{\code{signature(object = "HARModel")}: Extracts the coefficients from the Model}
\item{\code{SandwichNeweyWest}:}{\code{signature(object = "HARModel")}: Utilize the \code{sandwich} package to create newey west standard errors}
\item{\code{qlike}:}{\code{signature(object = "HARModel")}: Calculate the in sample 'qlike' loss function for a HARModel object}
\item{\code{logLik}:}{A wrapper for the "lm" subclass of the HARModel object}
\item{\code{confint}:}{A wrapper for the "lm" subclass of the HARModel object}
\item{\code{residuals}:}{A wrapper for the "lm" subclass of the HARModel object}
\item{\code{summary}:}{A wrapper for the "lm" subclass of the HARModel object}

}}
\author{Emil Sjoerup}
\keyword{classes}