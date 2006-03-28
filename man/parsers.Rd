\name{File parsers}
\alias{read.xysfiles}
\alias{read.celfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to parse CEL and XYS files.}
\description{
  These functions provide parsing capabilities to the oligo package,
  making it possible to read CEL and XYS files.
}
\usage{
read.xysfiles(filenames, phenoData = new("phenoData"), description = NULL, notes = "", verbose = FALSE)
read.celfiles(filenames, phenoData = new("phenoData"), description = NULL, notes = "", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{character vector containing the names of the files to
  be read.}
  \item{phenoData}{\code{phenoData} object.}
  \item{description}{\code{MIAME} object.}
  \item{notes}{notes.}
  \item{verbose}{logical. Verbose flag.}
}
\value{
  The functions return an \code{oligoBatch} object.
}
\examples{
## Do not run
## oBatch <- read.celfiles(list.celfiles())
}
\keyword{methods}% at least one, from doc/KEYWORDS
