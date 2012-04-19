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
\alias{mm<-,FeatureSet,ANY,ANY,matrix-method}
\alias{pm<-,FeatureSet,ANY,ANY,matrix-method}
\alias{bg<-,FeatureSet,matrix-method}
\alias{mm<-,FeatureSet,ANY,ANY,ff_matrix-method}
\alias{pm<-,FeatureSet,ANY,ANY,ff_matrix-method}
\alias{bg<-,FeatureSet,ff_matrix-method}
\alias{mm<-,TilingFeatureSet,ANY,ANY,array-method}
\alias{pm<-,TilingFeatureSet,ANY,ANY,array-method}
\alias{bg<-,TilingFeatureSet,array-method}

%% OK!

\title{Accessors and replacement methods for the intensity/PM/MM/BG matrices.}
\description{
  Accessors and replacement methods for the PM/MM/BG matrices.
}
\usage{
intensity(object)
mm(object, subset = NULL, target='core')
pm(object, subset = NULL, target='core')
bg(object, subset = NULL)
mm(object, subset = NULL, target='core')<-value
pm(object, subset = NULL, target='core')<-value
bg(object)<-value
}

\arguments{
  \item{object}{\code{FeatureSet} object.}
  \item{subset}{Not implemented yet.}
  \item{value}{\code{matrix} object.}
  \item{target}{One of 'probeset', 'core', 'full', 'extended'. This is
  ignored if the array design is something other than Gene ST or Exon ST.}
}

\details{
  For all objects but \code{TilingFeatureSet}, these methods will
  return matrices. In case of \code{TilingFeatureSet} objects, the
  value is a 3-dimensional array (probes x samples x channels).

  \code{intensity} will return the whole intensity matrix associated to
  the object. \code{pm}, \code{mm}, \code{bg} will return the respective
  PM/MM/BG matrix.

  When applied to \code{ExonFeatureSet} or \code{GeneFeatureSet}
  objects, \code{pm} will return the PM matrix at the transcript level
  ('core' probes) by default. The user should set the \code{target}
  argument accordingly if something else is desired. The valid values
  are: 'probeset' (Exon and Gene arrays), 'core' (Exon and Gene arrays),
  'full' (Exon arrays) and 'extended' (Exon arrays).

  The \code{target} argument has no effects when used on designs other
  than Gene and Exon ST.
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
