\name{mm}
\alias{mm}
\alias{pm}
\alias{bg}
\alias{intensity}
\alias{intensity,FeatureSet-method}
\alias{mm,FeatureSet-method}
\alias{pm,FeatureSet-method}
\alias{bg,FeatureSet-method}
\alias{mm,TilingFeatureSet-method}
\alias{pm,TilingFeatureSet-method}
\alias{bg,TilingFeatureSet-method}
\alias{intensity<-}
\alias{mm<-}
\alias{pm<-}
\alias{bg<-}
\alias{intensity<-,FeatureSet-method}
\alias{mm<-,FeatureSet,matrix-method}
\alias{pm<-,FeatureSet,matrix-method}
\alias{bg<-,FeatureSet,matrix-method}
\alias{mm<-,FeatureSet,ff_matrix-method}
\alias{pm<-,FeatureSet,ff_matrix-method}
\alias{bg<-,FeatureSet,ff_matrix-method}
\alias{mm<-,TilingFeatureSet,array-method}
\alias{pm<-,TilingFeatureSet,array-method}
\alias{bg<-,TilingFeatureSet,array-method}

%% OK!

\title{Accessors and replacement methods for the intensity/PM/MM/BG matrices.}
\description{
  Accessors and replacement methods for the PM/MM/BG matrices.
}
\usage{
intensity(object)
mm(object, subset = NULL)
pm(object, subset = NULL, ...)
bg(object, subset = NULL)
mm(object)<-value
pm(object)<-value
bg(object)<-value
}

\arguments{
  \item{object}{\code{FeatureSet} object.}
  \item{subset}{Not implemented yet.}
  \item{value}{\code{matrix} object.}
  \item{...}{Extra arguments for future implementation.}
}

\details{
  For all objects but \code{TilingFeatureSet}, these methods will
  return matrices. In case of \code{TilingFeatureSet} objects, the
  value is a 3-dimensional array (probes x samples x channels).

  \code{intensity} will return the whole intensity matrix associated to
  the object. \code{pm}, \code{mm}, \code{bg} will return the respective
  PM/MM/BG matrix.
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
