\name{Array dimensions}
\alias{nrow}
\alias{ncol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for the dimensions of the array}
\description{
  Accessors for the array dimensions.
}
\usage{
nrow(x)
ncol(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{oligoBatch object.}
}
\value{
  \code{ncol} and \code{nrow} return the number of columns and number of
  rows of the array associated to the oligoBatch object.
}
\examples{
## DO NOT RUN
## oBatch <- read.celfiles(list.celfiles())
## ncol(oBatch)
}
\keyword{methods}% at least one, from doc/KEYWORDS
