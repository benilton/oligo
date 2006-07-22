\name{snprma}
\alias{snprma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust multi-array average for SNP arrays}
\description{
  Performs robust multi-array average (RMA) method for SNP Arrays.
}
\usage{
snprma(oBatch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oBatch}{SnpFeatureSet object}
}
\details{
  This is an implementation of RMA for SNP arrays. It first preprocess
  the data, correcting for fragment length and average intensities and
  then summarizes at the SNP/Allele/Strand level.
}
\value{
  SnpQSet object.
}
\author{Benilton S Carvalho, Rafael A Irizarry}
\seealso{\code{\link{rma}}}
\examples{
## DO NOT RUN
## snpData <- read.celfiles(list.celfiles())
## snpQSet <- snprma(snpData)
}
\keyword{methods}% at least one, from doc/KEYWORDS

