\name{read.xysfiles}
\alias{read.xysfiles}
\title{Parser to XYS files}
\description{
  NimbleGen provides XYS files which are read by this function.
}
\usage{
read.xysfiles(..., filenames, pkgname, phenoData, featureData,
experimentData, notes, verbose=TRUE, sampleNames, checkType=TRUE)
}
\arguments{
  \item{...}{file names}
  \item{filenames}{\code{character} vector with filenames.}
  \item{pkgname}{\code{character} vector with alternative PD Info
    package name}
  \item{phenoData}{\code{phenoData}}
  \item{featureData}{\code{featureData}}
  \item{experimentData}{\code{experimentData}}
  \item{notes}{\code{notes}}
  \item{verbose}{\code{verbose}}
  \item{sampleNames}{\code{character} vector with sample names (usually
    better descriptors than the filenames)}
  \item{checkType}{\code{logical}. Check type of each file? This can be
    time consuming.}
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
xysPath <- system.file("extdata", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
}
}
\keyword{IO}
