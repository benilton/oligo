\name{probeNames}
\alias{probeNames}
\alias{probesetNames}
\alias{probeNames,FeatureSet-method}
\alias{probesetNames,FeatureSet-method}
\alias{probeNames,ExonFeatureSet-method}
\alias{probeNames,GeneFeatureSet-method}
\alias{probeNames,DBPDInfo-method}
\title{Accessor to feature names}
\description{
  Accessors to featureset names.
}
\usage{
probeNames(object, subset = NULL)
probesetNames(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{FeatureSet} or \code{DBPDInfo}}
  \item{subset}{not implemented yet.}
}
\value{
  \code{probeNames} returns a string with the probeset names for *each probe*
  on the array. \code{probesetNames}, on the other hand, returns the
  *unique probeset names*.
}
\keyword{manip}
