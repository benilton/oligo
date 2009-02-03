\name{read.xysfiles}
\alias{read.xysfiles}
\title{Parser to XYS files}
\description{
  NimbleGen provides XYS files which are read by this function.
}
\usage{
read.xysfiles(filenames, pkgname = NULL, phenoData = NULL, featureData = NULL, experimentData = NULL, notes = NULL, verbose = FALSE)
}
\arguments{
  \item{filenames}{\code{character} vector with filenames.}
  \item{pkgname}{\code{character} vector with alternative PD Info
    package name}
  \item{phenoData}{\code{phenoData}}
  \item{featureData}{\code{featureData}}
  \item{experimentData}{\code{experimentData}}
  \item{notes}{\code{notes}}
  \item{verbose}{\code{verbose}}
}
\details{
  The function will read the XYS files provided by NimbleGen Systems and
  return an object of class FeatureSet.
}
\value{
  \item{\code{ExpressionFeatureSet}}{if Expresssion arrays}
  \item{\code{TilingFeatureSet}}{if Tiling arrays}
}
\seealso{\code{\link{list.xysfiles}}, \code{\link{read.celfiles}}}
\examples{
if (require(maqcExpression4plex) & require(pd.hg18.60mer.expr)){
xysPath <- system.file("xysFiles", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
}
}
\keyword{IO}
