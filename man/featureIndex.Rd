\name{featureIndex}
\alias{featureIndex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Indexes of features on the array}
\description{
  Returns the indexes of features on the array.
}
\usage{
featureIndex(object, which=c("both","pm","mm"), genenames=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{FeatureSet} object.}
  \item{which}{A string defining what type of probes should have the
    indexes returned.}
  \item{genenames}{A character vector descbribing what feature sets the
    results should be restricted to.}
}
\seealso{\code{\link{pmindex}}, \code{\link{mmindex}}}
\keyword{manip}% at least one, from doc/KEYWORDS
