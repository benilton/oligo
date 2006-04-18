\name{bg.correct}
\alias{bg.correct.rma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Background Correction}
\description{
  Background corrects probe intensities in an object of class \code{FeatureSet}.
}
\usage{
bg.correct.rma(object, bgtype = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{FeatureSet}.}
  \item{bgtype}{For now set to 2.}
}
\details{
  the model based correction used by the RMA expression measure.
}
\value{
  An \code{FeatureSet} for which the intensities have been background
  adjusted. For some methods (RMA), only PMs are corrected and the MMs
  remain the same.
}
\keyword{manip}% at least one, from doc/KEYWORDS
