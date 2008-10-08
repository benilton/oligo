\name{antisenseThetaA}
\alias{antisenseThetaA}
\alias{antisenseThetaB}
\alias{senseThetaA}
\alias{senseThetaB}
\alias{antisenseThetaA,SnpQSet-method}
\alias{antisenseThetaB,SnpQSet-method}
\alias{senseThetaA,SnpQSet-method}
\alias{senseThetaB,SnpQSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarized Measures for SNP}
\description{
  SNPRMA summarizes the intensities for each SNP into 4
  quantities. These measures are at the log2 scale and are used to
  create log-ratios and log-averages, which are used by the genotype
  algorithm and copy number tool.
}

\usage{
antisenseThetaA(object)
antisenseThetaB(object)
senseThetaA(object)
senseThetaB(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpQSet} object}
}

\details{
  For a given SNP, an Affymetrix SNP array usually contains multiple
  features for each combination allele x strand, on which the summaries
  are computed.
}
\value{
  A matrix with R rows (SNP's) x C columns (Samples).
}

\references{}
\author{}
\note{
}

\seealso{\code{\link{snprma}}, \code{\link{crlmm}}}
\examples{
}

\keyword{manip}
