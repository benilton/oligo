\name{pmFragmentLength}
\alias{pmFragmentLength}
\alias{pmFragmentLength,AffySNPPDInfo-method}
\alias{pmFragmentLength,SnpFeatureSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Access the fragment length for PM probes.}
\description{
  Accessor to the fragment length for PM probes.
}
\usage{
pmFragmentLength(object, enzyme, type=c('snp', 'cn'))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{PDInfo} or \code{SnpFeatureSet} object.}
  \item{enzyme}{Enzyme to be used for query. If missing, all enzymes are
    used.}
  \item{type}{Type of probes to be used: 'snp' for SNP probes; 'cn' for
    Copy Number probes.}
}
\note{
  There is not a 1:1 relationship between probes and expected fragment
  length. For one enzyme, a given probe may be associated to multiple
  fragment lengths. Therefore, the number of rows in the data.frame may
  not match the number of PM probes and the \code{row} column should be
  used to match the fragment length with the PM matrix.
}
\value{
  A list of length equal to the number of enzymes used for
  digestion. Each element of the list is a data.frame containing:
  \itemize{
    \item \code{row}: the row used to link to the PM matrix;
    \item \code{length}: expected fragment length.
  }
}

\keyword{manip}
