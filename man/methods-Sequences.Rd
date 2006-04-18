\name{Methods for Sequences}
\alias{pmSequence}
\alias{mmSequence}
\alias{sequenceDesignMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Methods to deal with sequences.}
\description{
  Methods to access and create design matrices for sequences.
}
\usage{
pmSequence(x)
mmSequence(x)
sequenceDesignMatrix(seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{FeatureSet object.}
  \item{seqs}{character vector containing sequences of same length.}
}
\value{
  \code{pmSequence} and \code{mmSequence} return a character vector with
  the sequences for the PM and MM probes respectively.

  \code{sequenceDesignMatrix} returns a matrix of indicator variables
  with \code{length(seqs)} lines and \code{3*nchar(seqs[1])}
  columns. The name of the columns use the format X.Y standing for the
  indicator variable of base X at position Y.
}
\examples{
## Do not run
##
## FeatureSetObject <- read.celfiles(list.celfiles())
## pmSeqs <- pmSequence(FeatureSetObject)
## dMatrix <- sequenceDesignMatrix(pmSeqs)
}
\keyword{methods}% at least one, from doc/KEYWORDS
