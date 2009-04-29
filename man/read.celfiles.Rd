\name{read.celfiles}
\alias{read.celfiles}
\alias{read.celfiles2}
\alias{cleanPlatformName}
\title{Parser to CEL files}
\description{
  Reads CEL files.
}
\usage{
read.celfiles(..., filenames, pkgname, phenoData, featureData,
experimentData, notes, verbose = TRUE, sampleNames, rm.mask = FALSE,
rm.outliers = FALSE, rm.extra = FALSE, sd = FALSE, checkType = TRUE,
useAffyio = TRUE)

read.celfiles2(channel1, channel2, pkgname, phenoData, featureData,
experimentData, notes, verbose = TRUE, sampleNames, rm.mask = FALSE,
rm.outliers = FALSE, rm.extra = FALSE, sd = FALSE, checkType = TRUE,
useAffyio = TRUE)
}

\arguments{
  \item{...}{names of files to be read.}
  \item{filenames}{a \code{character} vector with the CEL filenames.}
  \item{channel1}{a \code{character} vector with the CEL filenames for
    the first 'channel' on a Tiling application}
  \item{channel2}{a \code{character} vector with the CEL filenames for
    the second 'channel' on a Tiling application}
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
  \item{useAffyio}{\code{logical}. Use 'affyio' instead of 'affxparser'
    to read in CEL files.}
}
\value{
  \item{\code{ExpressionFeatureSet}}{if Expresssion arrays}
  \item{\code{ExonFeatureSet}}{if Exon arrays}
  \item{\code{SnpFeatureSet}}{if SNP arrays}
  \item{\code{TilingFeatureSet}}{if Tiling arrays}
}
\details{
  When using 'affyio' to read in CEL files, the user can read compressed
  CEL files (CEL.gz). Additionally, 'affyio' is much faster than
  'affxparser'.

    The function guesses which annotation package to use from the header
  of the CEL file. The user can also provide the name of the annotaion
  package to be used (via the \code{pkgname} argument). If the
  annotation package cannot be loaded, the function returns an
  error. If the annotation package is not available from BioConductor,
  one can use the \code{pdInfoBuilder} package to build one.
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
