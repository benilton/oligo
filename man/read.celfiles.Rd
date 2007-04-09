\name{read.celfiles}
\alias{read.celfiles}
\alias{cleanPlatformName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parser to CEL files}
\description{
  Reads CEL files.
}
\usage{
read.celfiles(filenames, pkgname = NULL, phenoData = new("AnnotatedDataFrame"), featureData = NULL, description = NULL, notes = "", verbose = FALSE, rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, tmpdir = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{a \code{character} vector with the CEL filenames.}
  \item{pkgname}{alternative data package to be loaded.}
  \item{phenoData}{\code{phenoData}}
  \item{featureData}{\code{featureData}}
  \item{description}{\code{description}}
  \item{notes}{\code{notes}}
  \item{verbose}{\code{logical}}
  \item{rm.mask}{\code{logical}}
  \item{rm.outliers}{\code{logical}}
  \item{rm.extra}{\code{logical}}
  \item{tmpdir}{path to the directory that will store the temporary
    files.}
}
\details{
  \code{tmpdir} should be a local directory where temporary files
  (created by \code{BufferedMatrix} package) are going to be stored. If
  this points to a network drive, the process will be very slow. If it
  points to a directory where the user does not have 'write'
  permissions, the user will get errors.
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
