\name{rma}
\alias{rma}
\alias{rma,FeatureSet-method}
\title{RMA}
\description{
  Robust Multichip Average algorithm.
}
\usage{
rma(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{FeatureSet} object}
  \item{\dots}{the dots}
}
\details{
  *NOT* recommended for SNP arrays, see \code{\link{snprma}}.
}
\value{
  \code{ExpressionSet} object.
}
\seealso{\code{\link{snprma}}}
\examples{
if (require(maqcExpression4plex) & require(pd.hg18.60mer.expr)){
xysPath <- system.file("xysFiles", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
summarized <- rma(ngsExpressionFeatureSet)
show(summarized)
}
}
\keyword{manip}

