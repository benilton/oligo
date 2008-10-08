\name{getA}
\alias{getA}
\alias{getM}
\alias{getA,SnpQSet-method}
\alias{getM,SnpQSet-method}
\alias{getA,SnpCnvQSet-method}
\alias{getM,SnpCnvQSet-method}
\alias{getM,TilingQSet-method}
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
  \item{object}{\code{SnpQSet} object.}
}

\details{

  SNPRMA summarizes the SNP information into 4 quantities (log2-scale):
  \item{antisenseThetaA}{antisense allele A}
  \item{antisenseThetaB}{antisense allele B}
  \item{senseThetaA}{sense allele A}
  \item{senseThataB}{sense allele B}

  The average log-intensities are given by:
  \code{(antisenseThetaA+antisenseThetaB)/2} and
  \code{(senseThetaA+senseThetaB)/2}.

  The average log-ratios are given by:
  \code{antisenseThetaA-antisenseThetaB} and
  \code{senseThetaA-senseThetaB}.
  
}
\value{
  A 3-dimensional array (SNP's x Samples x Strand) with the requested
  measure.

  If a Genomewide array (like Affy SNP 6.0) is used, the strand
  dimension described above is dropped.
}

\references{}
\author{}
\note{}
\seealso{\code{\link{snprma}}}
\examples{
}

\keyword{manip}

