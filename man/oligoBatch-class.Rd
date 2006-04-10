\name{oligoBatch-class}
\docType{class}
\alias{oligoBatch-class}
\alias{allele,oligoBatch-method}
\alias{boxplot,oligoBatch-method}
\alias{featureIndex,oligoBatch-method}
\alias{geneNames,oligoBatch-method}
\alias{getPlatformDesign,oligoBatch-method}
\alias{hist,oligoBatch-method}
\alias{image,oligoBatch-method}
\alias{indexFeatureSetName,oligoBatch-method}
\alias{length,oligoBatch-method}
\alias{mmindex,oligoBatch-method}
\alias{mm<-,oligoBatch-method}
\alias{mm,oligoBatch-method}
\alias{ncol,oligoBatch-method}
\alias{npixels,oligoBatch-method}
\alias{nrow,oligoBatch-method}
\alias{platformDesignName,oligoBatch-method}
\alias{platform,oligoBatch-method}
\alias{pmindex,oligoBatch-method}
\alias{pm<-,oligoBatch-method}
\alias{pm,oligoBatch-method}
\alias{probeNames,oligoBatch-method}
\alias{sd,oligoBatch-method}
\alias{show,oligoBatch-method}
\alias{$.platformDesign, oligoBatch-method}

\title{Class "oligoBatch"}
\description{
  This is a class representation for oligonucleotide arrays at probe level. The main component are the intensities from multiple arrays of
  the same design file. It is an extension of the \code{\link{eSet-class}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("oligoBatch", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"}
      describing the manufacturer.}
    \item{\code{platform}:}{Object of class \code{"character"}
      describing the design file associated to the object.}
    \item{\code{assayData}:}{Object of class \code{"listOrEnv"} as
      described in \code{\link{eSet-class}}.}
    \item{\code{sampleNames}:}{Object of class \code{"character"} with
      the sample names.}
    \item{\code{reporterNames}:}{Object of class \code{"character"} with
    the reporter names.}
    \item{\code{description}:}{Object of class \code{"MIAME"}}
    \item{\code{notes}:}{Object of class \code{"character"} with notes.}
    \item{\code{annotation}:}{Object of class \code{"character"} with annotation.}
    \item{\code{history}:}{Object of class \code{"character"} with history.}
    \item{\code{reporterInfo}:}{Object of class \code{"data.frameOrNULL"}}
    \item{\code{phenoData}:}{Object of class \code{"phenoData"}}
  }
}
\section{Extends}{
Class \code{"eSet"}, directly.
Class \code{"annotatedDataset"}, by class \code{"eSet"}.
}
\section{Methods}{
  \describe{
    \item{allele}{\code{signature(object = "oligoBatch")}: extracts
      allele information for SNP arrays.}
    \item{boxplot}{\code{signature(x = "oligoBatch")}: produces boxplots.}
    \item{featureIndex}{\code{signature(object = "oligoBatch")}: Returns
    the indexes of features on the array.}
    \item{geneNames}{\code{signature(object = "oligoBatch")}: Returns
      the names of the (PM) features.}
    \item{getPlatformDesign}{\code{signature(object = "oligoBatch")}:
      Returns the Platform Design Environment associated to the object.}
    \item{hist}{\code{signature(x = "oligoBatch")}: Produces histograms.}
    \item{image}{\code{signature(x = "oligoBatch")}: Produces image plots.}
    \item{indexFeatureSetName}{\code{signature(object = "oligoBatch")}:
      Returns the indexes in the PDEnv for given features.}
    \item{length}{\code{signature(x = "oligoBatch")}: Returns the number
    of samples in the object.}
    \item{mmindex}{\code{signature(object = "oligoBatch")}: Returns the
      indexes for the MM features.}
    \item{mm<-}{\code{signature(object = "oligoBatch")}: Assigns a value
    to the MM matrix.}
    \item{mm}{\code{signature(object = "oligoBatch")}: Returns the MM matrix.}
    \item{ncol}{\code{signature(x = "oligoBatch")}: Returns the number
      of columns of the array.}
    \item{npixels}{\code{signature(object = "oligoBatch")}: Returns the
      matrix with the number of pixels used.}
    \item{nrow}{\code{signature(x = "oligoBatch")}: Returns the number
      of rows of the array.}
    \item{platformDesignName}{\code{signature(object = "oligoBatch")}:
      Returns the name of the PDEnv associated.}
    \item{platform}{\code{signature(object = "oligoBatch")}: Returns the
    type of the array.}
    \item{pmindex}{\code{signature(object = "oligoBatch")}: Returns the
      indexes for the PM features.}
    \item{pm<-}{\code{signature(object = "oligoBatch")}: Assigns a value
    to the PM matrix.}
    \item{pm}{\code{signature(object = "oligoBatch")}: Returns the PM matrix.}
    \item{probeNames}{\code{signature(object = "oligoBatch")}: Returns
      the names of the features.}
    \item{sd}{\code{signature(x = "oligoBatch")}: returns the matrix of
      standard deviations.}
    \item{show}{\code{signature(object = "oligoBatch")}: Display basic
      information about the array.}
	 }
}

\seealso{
  \code{\link{eSet-class}}
}

\keyword{classes}
