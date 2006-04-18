\name{probeMatch-methods}
\alias{pm}
\alias{mm}
\alias{pm<-}
\alias{mm<-}
\alias{pmindex}
\alias{mmindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for perfect match (pm) and mismatch (mm) probes.}
\description{
  Methods for perfect match and mismatch probes.
}
\usage{
pm(object, genenames=NULL)
pmindex(object)
pm(object)<-value
mm(object, genenames=NULL)
mmindex(object)
mm(object)<-value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{FeatureSet object}
  \item{genenames}{character vector with the featureset names of
    interest}
  \item{value}{new values to replace the PM/MM matrices.}
}
\value{
  \code{pm} and \code{mm} will return respectively the PM and MM
  matrices associated to the FeatureSet object.

  \code{pmindex} and \code{mmindex} will return vectors containing the
  indexes for the PM and MM probes.
}
\examples{
## Do not run
## FeatureSetObject <- read.celfiles(list.celfiles())
## pms <- pm(FeatureSetObject)
}
\keyword{methods}% at least one, from doc/KEYWORDS
