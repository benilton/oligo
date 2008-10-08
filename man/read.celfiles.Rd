\name{read.celfiles}
\alias{read.celfiles}
\alias{cleanPlatformName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parser to CEL files}
\description{
  Reads CEL files.
}
\usage{

read.celfiles(filenames, pkgname, phenoData, featureData, experimentData, notes = "", verbose = TRUE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{a \code{character} vector with the CEL filenames.}
  \item{pkgname}{alternative data package to be loaded.}
  \item{phenoData}{\code{phenoData}}
  \item{featureData}{\code{featureData}}
  \item{experimentData}{\code{experimentData}}
  \item{notes}{\code{notes}}
  \item{verbose}{\code{logical}}
}
\value{
  \item{\code{ExpressionFeatureSet}}{if Expresssion arrays}
  \item{\code{ExonFeatureSet}}{if Exon arrays}
  \item{\code{SnpFeatureSet}}{if SNP arrays}
  \item{\code{TilingFeatureSet}}{if Tiling arrays}
}

\references{}
\author{}
\note{}
\seealso{\code{\link{read.xysfiles}} for NimbleGen arrays}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}
