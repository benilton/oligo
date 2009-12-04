\name{chromosome}
\alias{chromosome}
\alias{chromosome<-}
\alias{chromosome<-,AnnotatedDataFrame,character-method}
\alias{pmChr}
\alias{pmChr,FeatureSet-method}
\title{Accessor for chromosome information}
\description{
  Returns chromosome information.
}
\usage{
%- chromosome(object)
pmChr(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{TilingFeatureSet} or \code{SnpCallSet} object}
}
\details{
  \code{chromosome()} returns the chromosomal information for all probes
  and \code{pmChr()} subsets the output to the PM probes only (if a
  TilingFeatureSet object).
}
\value{
  Vector with chromosome information.
}
\keyword{manip}
