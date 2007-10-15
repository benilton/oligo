\name{read.xysfiles}
\alias{read.xysfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parser to XYS files}
\description{
  NimbleGen provides XYS files which are read by this function.
}
\usage{
read.xysfiles(filenames, pkgname = NULL, phenoData = NULL, featureData = NULL, experimentData = NULL, notes = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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

\references{}
\author{}
\note{}
\seealso{\code{\link{read.celfiles}} for Affymetrix arrays}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
