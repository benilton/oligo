\name{platformDesign-class}
\docType{class}
\alias{platformDesign-class}
\alias{mmindex,platformDesign-method}
\alias{names,platformDesign-method}
\alias{nProbes,platformDesign-method}
\alias{pmindex,platformDesign-method}
\alias{probeNames,platformDesign-method}
\alias{show,platformDesign-method}
\alias{$.platformDesign}

\title{Class "platformDesign"}
\description{A class to store array related information.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("platformDesign", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{featureInfo}:}{Object of class \code{"environment"}
      describing the array at the feature level.}
    \item{\code{featureTypeDescription}:}{Object of class \code{"list"}
      describing the feature types.}
    \item{\code{manufacturer}:}{Object of class \code{"character"}.}
    \item{\code{type}:}{Object of class \code{"character"} describing
      the array type.}
    \item{\code{nrow}:}{Object of class \code{"numeric"} describing the
      number of rows of the array.}
    \item{\code{ncol}:}{Object of class \code{"numeric"} describing the
      number of columns of the array.}
    \item{\code{nwells}:}{Object of class \code{"numeric"} describing
      the number of wells per arrays, not used, but allocated for future
    extensions.}
    \item{\code{lookup}:}{Object of class \code{"data.frame"} describing
    location on the array for multiple-wells chips. Not used for now.}
    \item{\code{indexes}:}{Object of class \code{"list"}}
    \item{\code{platforms}:}{Object of class \code{"character"}.}
  }
}
\section{Methods}{
  \describe{
    \item{mmindex}{\code{signature(object = "platformDesign")}: indexes
      for MMs.}
    \item{names}{\code{signature(x = "platformDesign")}: variable names
      available in the PDEnv.}
    \item{nProbes}{\code{signature(object = "platformDesign")}: number
      of features.}
    \item{pmindex}{\code{signature(object = "platformDesign")}: indexes
      of PMs.}
    \item{probeNames}{\code{signature(object = "platformDesign", subset
	= "characterOrNULL")}:
      feature names for PMs.}
    \item{show}{\code{signature(object = "platformDesign")}: basic
      description of the PDEnv.}
	 }
}
\keyword{classes}
