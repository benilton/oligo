\name{FeatureSet-class}
\docType{class}
\alias{FeatureSet-class}
\alias{allele,FeatureSet-method}
\alias{boxplot,FeatureSet-method}
\alias{featureIndex,FeatureSet-method}
\alias{geneNames,FeatureSet-method}
\alias{getPlatformDesign,FeatureSet-method}
\alias{hist,FeatureSet-method}
\alias{image,FeatureSet-method}
\alias{indexFeatureSetName,FeatureSet-method}
\alias{length,FeatureSet-method}
\alias{mmindex,FeatureSet-method}
\alias{mm<-,FeatureSet-method}
\alias{mm,FeatureSet-method}
\alias{ncol,FeatureSet-method}
\alias{npixels,FeatureSet-method}
\alias{nrow,FeatureSet-method}
\alias{platformDesignName,FeatureSet-method}
\alias{platform,FeatureSet-method}
\alias{pmindex,FeatureSet-method}
\alias{pm<-,FeatureSet-method}
\alias{pm,FeatureSet-method}
\alias{probeNames,FeatureSet-method}
\alias{sd,FeatureSet-method}
\alias{show,FeatureSet-method}
\alias{$.platformDesign, FeatureSet-method}
\alias{exprs, FeatureSet-method}
\alias{manufacturer, FeatureSet-method}
\alias{manufacturer<-, FeatureSet-method}
\alias{platform<-, FeatureSet-method}

\title{Class "FeatureSet"}
\description{
  This is a class representation for oligonucleotide arrays at probe level. The main component are the intensities from multiple arrays of
  the same design file. It is an extension of the \code{eSet-class}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FeatureSet", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"}
      describing the manufacturer.}
    \item{\code{platform}:}{Object of class \code{"character"}
      describing the design file associated to the object.}
    \item{\code{assayData}:}{Object of class \code{"listOrEnv"} as
      described in \code{eSet-class}.}
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
    \item{allele}{\code{signature(object = "FeatureSet")}: extracts
      allele information for SNP arrays.}
    \item{boxplot}{\code{signature(x = "FeatureSet")}: produces boxplots.}
    \item{featureIndex}{\code{signature(object = "FeatureSet")}: Returns
    the indexes of features on the array.}
    \item{geneNames}{\code{signature(object = "FeatureSet")}: Returns
      the names of the (PM) features.}
    \item{getPlatformDesign}{\code{signature(object = "FeatureSet")}:
      Returns the Platform Design Environment associated to the object.}
    \item{hist}{\code{signature(x = "FeatureSet")}: Produces histograms.}
    \item{image}{\code{signature(x = "FeatureSet")}: Produces image plots.}
    \item{indexFeatureSetName}{\code{signature(object = "FeatureSet")}:
      Returns the indexes in the PDEnv for given features.}
    \item{length}{\code{signature(x = "FeatureSet")}: Returns the number
    of samples in the object.}
    \item{mmindex}{\code{signature(object = "FeatureSet")}: Returns the
      indexes for the MM features.}
    \item{mm<-}{\code{signature(object = "FeatureSet")}: Assigns a value
    to the MM matrix.}
    \item{mm}{\code{signature(object = "FeatureSet")}: Returns the MM matrix.}
    \item{ncol}{\code{signature(x = "FeatureSet")}: Returns the number
      of columns of the array.}
    \item{npixels}{\code{signature(object = "FeatureSet")}: Returns the
      matrix with the number of pixels used.}
    \item{nrow}{\code{signature(x = "FeatureSet")}: Returns the number
      of rows of the array.}
    \item{platformDesignName}{\code{signature(object = "FeatureSet")}:
      Returns the name of the PDEnv associated.}
    \item{platform}{\code{signature(object = "FeatureSet")}: Returns the
      type of the array.}
    \item{platform<-}{\code{signature(object = "FeatureSet")}: Assigns a
      type for the array.}
    \item{pmindex}{\code{signature(object = "FeatureSet")}: Returns the
      indexes for the PM features.}
    \item{pm<-}{\code{signature(object = "FeatureSet")}: Assigns a value
    to the PM matrix.}
    \item{pm}{\code{signature(object = "FeatureSet")}: Returns the PM matrix.}
    \item{probeNames}{\code{signature(object = "FeatureSet")}: Returns
      the names of the features.}
    \item{sd}{\code{signature(x = "FeatureSet")}: returns the matrix of
      standard deviations.}
    \item{show}{\code{signature(object = "FeatureSet")}: Display basic
      information about the array.}
    \item{manufacturer}{\code{signature(object = "FeatureSet")}: returns
    the array manufacturer.}
    \item{manufacturer<-}{\code{signature(object = "FeatureSet")}: assigns
      the array manufacturer.}
    \item{exprs}{\code{signature(object = "FeatureSet")}: returns
      the matrix of intensities.}
%%    \item{initialize}{\code{signature(object = "FeatureSet")}:
%%      Initializes the object.}
  }
}

\seealso{
  \code{eSet-class}
}

\keyword{classes}
