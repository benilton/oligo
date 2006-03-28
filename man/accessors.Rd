\name{Accessors}
\alias{exprs}
\alias{npixels}
\alias{sd}
\alias{platform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for oligoBatch objects.}
\description{
  Set of accessors for \code{oligoBatch} objects.
}
\usage{
exprs(object)
npixels(object)
sd(object)
platform(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{oligoBatch} object.}
}
\value{
  \code{exprs} returns the matrix of intensities associated to the
  \code{oligoBatch} object.

  \code{npixels} returns the matrix of number of pixels associated to
  the \code{oligoBatch} object, if any.

  \code{sd} returns the matrix of standard-deviations associated to the
  \code{oligoBatch} object, if any.

  \code{platform} returns a string determing the platform of the array (SNP/expression).
}
\examples{
## Do not run
## oBatch <- read.celfiles(list.celfiles())
## intensities <- exprs(oBatch)
}
\keyword{methods}% at least one, from doc/KEYWORDS
