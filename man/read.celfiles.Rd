\name{read.celfiles}
\alias{read.celfiles}
\alias{cleanPlatformName}
\title{Parser to CEL files}
\description{
  Reads CEL files.
}
\usage{
read.celfiles(..., filenames, pkgname, phenoData, featureData,
experimentData, notes, verbose = TRUE, sampleNames, rm.mask = FALSE,
rm.outliers = FALSE, rm.extra = FALSE, sd = FALSE, checkType = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{names of files to be read.}
  \item{filenames}{a \code{character} vector with the CEL filenames.}
  \item{pkgname}{alternative data package to be loaded.}
  \item{phenoData}{\code{phenoData}}
  \item{featureData}{\code{featureData}}
  \item{experimentData}{\code{experimentData}}
  \item{notes}{\code{notes}}
  \item{verbose}{\code{logical}}
  \item{sampleNames}{\code{character} vector with sample names (usually
    better descriptors than the filenames)}
  \item{rm.mask}{\code{logical}. Read masked?}
  \item{rm.outliers}{\code{logical}. Remove outliers?}
  \item{rm.extra}{\code{logical}. Remove extra?}
  \item{sd}{\code{logical}. Read SD?}
  \item{checkType}{\code{logical}. Check type of each file? This can be
    time consuming.}
}
\value{
  \item{\code{ExpressionFeatureSet}}{if Expresssion arrays}
  \item{\code{ExonFeatureSet}}{if Exon arrays}
  \item{\code{SnpFeatureSet}}{if SNP arrays}
  \item{\code{TilingFeatureSet}}{if Tiling arrays}
}
\seealso{\code{\link{list.celfiles}}, \code{\link{read.xysfiles}}}
\examples{
if(require(pd.mapping50k.xba240) & require(hapmap100kxba)){
celPath <- system.file("celFiles", package="hapmap100kxba")
celFiles <- list.celfiles(celPath, full.name=TRUE)
affySnpFeatureSet <- read.celfiles(celFiles)
}
}
\keyword{IO}
