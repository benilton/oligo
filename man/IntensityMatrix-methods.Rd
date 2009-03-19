\name{mm}
\alias{mm}
\alias{pm}
\alias{bg}
\alias{mm,FeatureSet-method}
\alias{pm,FeatureSet-method}
\alias{bg,FeatureSet-method}
\alias{mm<-}
\alias{pm<-}
\alias{bg<-}
\alias{mm<-,FeatureSet,matrix-method}
\alias{pm<-,FeatureSet,matrix-method}
\alias{bg<-,FeatureSet,matrix-method}

%% OK!

\title{Accessors and replacement methods for the PM/MM/BG matrices.}
\description{
  Accessors and replacement methods for the PM/MM/BG matrices.
}
\usage{
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
\examples{
if (require(maqcExpression4plex) & require(pd.hg18.60mer.expr)){
xysPath <- system.file("extdata", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
pm(ngsExpressionFeatureSet)[1:10,]
}
}

\keyword{manip}
