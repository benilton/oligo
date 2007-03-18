\name{cnConfidence}
\alias{cnConfidence}
\alias{copyNumber}
\alias{cnConfidence<-}
\alias{copyNumber<-}
\alias{cnConfidence,SnpCopyNumberSet-method}
\alias{copyNumber,SnpCopyNumberSet-method}
\alias{cnConfidence<-,SnpCopyNumberSet,matrix-method}
\alias{copyNumber<-,SnpCopyNumberSet,matrix-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for SnpCopyNumberSet objects}
\description{
  Accessors and replacement methods for SnpCopyNumberSet objects.
}
\usage{
copyNumber(object)
cnConfidence(object)
copyNumber(object)<-value
cnConfidence(object)<-value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpCopyNumberSet} object}
  \item{value}{matrix}
}
\details{
}
\value{
  Matrix with R rows (SNP's) and C columns (Samples) containing copy
  number estimates and confidences.
}
\references{}
\author{}
\note{}
\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
