\name{calls}
\alias{calls}
\alias{calls,SnpCallSet-method}
\alias{calls<-}
\alias{calls<-,SnpCallSet,matrix-method}
\alias{callsConfidence}
\alias{callsConfidence,SnpCallSet-method}
\alias{callsConfidence<-}
\alias{callsConfidence<-,SnpCallSet,matrix-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for SnpCallSet}

\description{Accessors and replacement methods for genotype calls and
  associated confidences.}

\usage{
calls(object)
callsConfidence(object)
calls(object)<-value
callsConfidence(object)<-value
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpCallSet} object.}
  \item{value}{\code{matrix}}
}
\details{
}
\value{
  Matrix with R rows (SNP's) and C columns (Samples) with genotype calls
  (1 - AA; 2 - AB; 3 - BB) or confidences.
}

\references{}
\author{}
\note{}

\seealso{\code{\link{crlmm}} }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
