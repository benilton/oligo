\name{Array dimensions}
\alias{nrow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for the dimensions of the array}
\description{
  Accessors for the array dimensions.
}
\usage{
nrow(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{oligoBatch object.}
}
\value{
  \code{nrow} returns the number of columns and number of
  rows of the array associated to the oligoBatch object.
}
\seealso{\code{\link{ncol}}}
\examples{
## DO NOT RUN
## oBatch <- read.celfiles(list.celfiles())
## ncol(oBatch)
}
\keyword{methods}% at least one, from doc/KEYWORDS
