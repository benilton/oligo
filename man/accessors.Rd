\name{Accessors}
\alias{npixels}
\alias{type}
\alias{platform}
\alias{platform<-}
\alias{manufacturer}
\alias{manufacturer<-}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for FeatureSet objects.}
\description{
  Set of accessors for \code{FeatureSet} objects.
}
\usage{
npixels(object)
type(object)
platform(object)
platform(object)<-value
manufacturer(object)
manufacturer(object)<-value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{FeatureSet} object.}
  \item{value}{string to be assigned to the platform slot.}
}
\value{
  \code{exprs} returns the matrix of intensities associated to the
  \code{FeatureSet} object.

  \code{npixels} returns the matrix of number of pixels associated to
  the \code{FeatureSet} object, if any.

  \code{sd} returns the matrix of standard-deviations associated to the
  \code{FeatureSet} object, if any.

  \code{platform} returns a string determing the platform of the array (SNP/expression).
}
\examples{
## Do not run
## FeatureSetObject <- read.celfiles(list.celfiles())
## intensities <- exprs(FeatureSetObject)
}
\keyword{methods}% at least one, from doc/KEYWORDS
