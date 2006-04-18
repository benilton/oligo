\name{probeNames}
\alias{probeNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Feature Set Names}
\description{
  The function returns the feature set names for all PM probes.
}
\usage{
probeNames(object, subset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A \code{FeatureSet} object.}
  \item{subset}{a \code{characterOrNULL} object. The subset to be used (not implemented yet)}
}
\value{
  A character vector with the feature set names. Restricted to PM only.
}
\section{Warning}{This function returns the feature set names
  \bold{only} for PM probes.}
\keyword{manip}% at least one, from doc/KEYWORDS
