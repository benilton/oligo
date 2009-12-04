\name{pmPosition}
\alias{pmPosition}
\alias{pmPosition,FeatureSet-method}
\alias{pmPosition,TilingFeatureSet-method}
\alias{pmPosition,ExpressionPDInfo-method}
\alias{pmPosition,TilingPDInfo-method}
\alias{pmOffset}
\alias{pmOffset,AffySNPPDInfo-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessor to position information}
\description{
  \code{pmPosition} will return the genomic position for the
  (PM) probes.
}
\usage{
%position(object)
pmPosition(object)
pmOffset(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{AffySNPPDInfo}, \code{TilingFeatureSet} or
    \code{SnpCallSet} object}
}
\details{
  \code{pmPosition} will return genomic position for PM probes on a
  tiling array.

  \code{pmOffset} will return the offset information for PM probes on
  SNP arrays.
}
\keyword{manip}
