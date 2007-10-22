\name{pmPosition}
\alias{pmPosition}
\alias{pmOffset}
\alias{pmOffset,AffySNPPDInfo-method}
\alias{pmPosition,TilingFeatureSet-method}
\alias{position}
\alias{position,TilingFeatureSet-method}
\alias{position,SnpCallSet-method}
\alias{pmPosition,ExpressionPDInfo-method}
\alias{pmPosition,TilingPDInfo-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessor to position information}
\description{
  \code{pmPosition} will return the genomic position for the
  (PM) probes.
}
\usage{
position(object)
pmPosition(object)
pmOffset(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{AffySNPPDInfo}, \code{TilingFeatureSet} or
    \code{SnpCallSet} object}
}
\details{
  \code{position} will return genomic position for all probes on a
  \code{FeatureSet} or \code{SnpCallSet} object.

  \code{pmPosition} will return genomic position for PM probes on a
  tiling array.

  \code{pmOffset} will return the offset information for PM probes on
  SNP arrays.
}
\value{
}
\references{}
\author{}
\note{}
\seealso{}
\examples{
}

\keyword{manip}
