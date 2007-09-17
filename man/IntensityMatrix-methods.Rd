\name{mm}
\alias{mm}
\alias{pm}
\alias{mm,FeatureSet-method}
\alias{pm,FeatureSet-method}
\alias{pm,SnpCnvFeatureSet-method}
\alias{mm<-}
\alias{pm<-}
\alias{mm<-,FeatureSet,matrix-method}
\alias{pm<-,FeatureSet,matrix-method}
\alias{mm<-,FeatureSet,BufferedMatrix-method}
\alias{pm<-,FeatureSet,BufferedMatrix-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessor to the PM/MM matrices}
\description{
  Accessor to the PM/MM matrices.
}
\usage{
mm(object, subset = NULL)
pm(object, subset = NULL, ...)
mm(object)<-value
pm(object)<-value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{FeatureSet} object }
  \item{subset}{not implemented yet}
  \item{value}{\code{matrix} or \code{BufferedMatrix} object}
  \item{...}{extra arguments}
}
\details{
  The intensity matrix is stored as a \code{BufferedMatrix} object. This
  should leave more memory for computations.
}
\value{
  PM or MM matrix as requested as a \code{BufferedMatrix} object.
}

\references{}
\author{}
\note{}
\seealso{}
\examples{
}

\keyword{manip}
