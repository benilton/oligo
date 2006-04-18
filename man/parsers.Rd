\name{File parsers}
\alias{read.xysfiles}
\alias{readonexysfile}
\alias{readxysHeader}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{XYS files parser.}
\description{
  These functions provide parsing capabilities to the oligo package,
  making it possible to read XYS files.
}
\usage{
read.xysfiles(filenames, phenoData = new("AnnotatedDataFrame"), description =
NULL, notes = "", verbose = FALSE)
readonexysfile(filename)
readxysHeader(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{character vector containing the names of the files to
    be read.}
  \item{filename}{string with the name of the file}
  \item{phenoData}{\code{phenoData} object.}
  \item{description}{\code{MIAME} object.}
  \item{notes}{notes.}
  \item{verbose}{logical. Verbose flag.}
}
\details{
  The functions \code{readonexysfile} and \code{readxysHeader} are not
  expected to be used by the final user. They are mainly functions for
  internal use.
}
\value{
  The functions return an \code{FeatureSet} object.
}

\seealso{ \code{\link{list.celfiles}}, \code{\link{read.celfiles}}, \code{\link{list.xysfiles}} }
\examples{
## Do not run
## FeatureSetObject <- read.xysfiles(list.xysfiles())
}
\keyword{methods}% at least one, from doc/KEYWORDS
