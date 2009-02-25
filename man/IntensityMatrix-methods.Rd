\name{mm}
\alias{mm}
\alias{pm}
\alias{mm,FeatureSet-method}
\alias{pm,FeatureSet-method}
%\alias{pm,SnpCnvFeatureSet-method}
\alias{mm<-}
\alias{pm<-}
\alias{mm<-,FeatureSet,matrix-method}
\alias{pm<-,FeatureSet,matrix-method}
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
\examples{
if (require(maqcExpression4plex) & require(pd.hg18.60mer.expr)){
xysPath <- system.file("extdata", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
pm(ngsExpressionFeatureSet)[1:10,]
}
}

\keyword{manip}
