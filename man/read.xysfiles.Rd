\name{read.xysfiles}
\alias{read.xysfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parser to XYS files}
\description{
  NimbleGen provides XYS files which are read by this function.
}
\usage{
read.xysfiles(filenames, phenoData = new("AnnotatedDataFrame"), featureData = NULL, description = NULL, notes = "", verbose = FALSE, tmpdir = tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{\code{character} vector with filenames.}
  \item{phenoData}{\code{phenoData}}
  \item{featureData}{\code{featureData}}
  \item{description}{\code{description}}
  \item{notes}{\code{notes}}
  \item{verbose}{\code{verbose}}
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
