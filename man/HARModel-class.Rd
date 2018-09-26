\name{HARModel-class}
\alias{HARModel-class}
\alias{show,HARModel-method}
\alias{plot,HARModel,missing-method}
\alias{uncondmean,HARModel-method}
\alias{coef,HARModel-method}
\alias{SandwichNeweyWest,HARModel-method}
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
    \item{\code{Model}:}{Object of class \code{lm}. Contains the linear model fitted to which \code{mVarCovar}, the Newey-West variance-covariance matrix, is appended }
    \item{\code{Info}:}{Object of class \code{list} cointaining:}
    \itemize{\item{\code{NWLagOrder}:} Integer denoting the used Lag order in the Newey-West Variance Covariance estimation.
    \item{\code{Lags}:} \code{numeric} containing the lags used to create the model
    \item{\code{Dates}:} \code{Date} object containing the dates for which the estimation was done, only applicable if the Model was estimated using an "xts" object.}
    \item{\code{Data}:}{\code{list} object containing the provided data.}
    }
}
\section{Methods}{
\itemize{
\item{\code{show}: Shows summary}
\item{\code{plot}: Plots the HARModel}
\item{\code{uncondmean}: Extracts the unconditional mean from the Model}
\item{\code{coef}: Extracts the coefficients from the Model}
\item{\code{SandwichNeweyWest}: Utilize the \code{sandwich} package to create newey west standard errors}
}}
\author{Emil Sjoerup}
\keyword{classes}