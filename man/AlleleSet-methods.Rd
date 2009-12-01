\name{getA}
\alias{getA}
\alias{getM}
\alias{getA,AlleleSet-method}
\alias{getM,AlleleSet-method}
\alias{getM,TilingFeatureSet-method}
\alias{getA,TilingFeatureSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute average log-intensities / log-ratios}
\description{
  Methods to compute average log-intensities and log-ratios across
  alleles, within strand.
}
\usage{
getA(object)
getM(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{AlleleSet} or \code{TilingFeatureSet} object.}
}

\details{

  For SNP data, SNPRMA summarizes the SNP information into 4 quantities (log2-scale):
  \itemize{
    \item{antisenseThetaA}{antisense allele A}
    \item{antisenseThetaB}{antisense allele B}
    \item{senseThetaA}{sense allele A}
    \item{senseThataB}{sense allele B}
  }
  
  The average log-intensities are given by:
  \code{(antisenseThetaA+antisenseThetaB)/2} and
  \code{(senseThetaA+senseThetaB)/2}.

  The average log-ratios are given by:
  \code{antisenseThetaA-antisenseThetaB} and
  \code{senseThetaA-senseThetaB}.

  For Tiling data, \code{getM} and \code{getA} return the log-ratio and
  average log-intensities computed across channels:
  \code{M = log2(channel1)-log2(channel2)}
  \code{A = (log2(channel1)+log2(channel2))/2}
}
\value{
  A 3-dimensional array (SNP's x Samples x Strand) with the requested
  measure, when the input SNP data (50K, 250K).

  A 2-dimensional array (SNP's x Samples), when the input is from SNP
  5.0 and SNP 6.0 arrays.

  A 2-dimensional array if the input is from Tiling arrays.
}

\seealso{\code{\link[oligo]{snprma}}}
\keyword{manip}

