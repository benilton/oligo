\name{File Parsers}
\alias{read.celfiles}
\alias{read.affybatch}
\alias{stuffForXYSandCELreaders}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Affymetrix CEL files parser}
\description{
  Read Affymetrix CEL files into an oligoBatch.
}
\usage{
read.celfiles(filenames,
              pdenv=TRUE,
              arrayType=NULL,
              pkgname=NULL,
              sd=FALSE,
              npixels=FALSE,
              phenoData=new("phenoData"),
              description=NULL,
              notes="",
              verbose=FALSE,
              compress=FALSE,
              rm.mask=FALSE,
              rm.outliers=FALSE,
              rm.extra=FALSE)
read.affybatch(..., filenames = character(0), phenoData = new("phenoData"),
          description = NULL, notes = "", compress=getOption("BioC")$affy$compress.cel,
          rm.mask = FALSE, rm.outliers = FALSE, rm.extra = FALSE, verbose =FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{a character vector containing the file names.}
  \item{pdenv}{should load the standard Platform Design Environment (pdenv)
    associated to the CEL files?}
  \item{arrayType}{type of array: "SNP" or "expression".}
  \item{pkgname}{string specifying the alternative PDEnv to be
    loaded. (see Details)}
  \item{sd}{should return standard deviations?}
  \item{npixels}{should return number of pixels?}
  \item{phenoData}{a \code{\link[Biobase]{phenoData}} object.}
  \item{description}{a \code{\link[Biobase]{MIAME}} object.}
  \item{notes}{notes.}
  \item{verbose}{verbosity flag}
  \item{compress}{are the CEL files compressed?}
  \item{rm.mask}{should the spots marked as 'MASKS' set to \code{NA}?}
  \item{rm.outliers}{should the spots marked as 'OUTLIERS' set to \code{NA}?}
  \item{rm.extra}{if \code{TRUE}, overrides what is in \code{rm.mask}
    and \code{rm.oultiers}}
  \item{\dots}{extra arguments to be passed to \code{\link{list}}}
}

\details{
  \code{read.celfiles} is an Affymetrix CEL file parser.

  If \code{pdenv} is set to \code{FALSE}, then \code{arrayType} is required and
  expected to assume either "SNP" or "expression". This will allow for
  the correct assignment of the class of the object
  (affysnpBatch/affyexprsBatch).

  If \code{pkgname} is not \code{NULL}, then \code{pdenv} is
  automatically set to \code{FALSE}.

  The arguments \code{pdenv}, \code{arrayType}, \code{pkgname} and
  \code{uniquecdf} are meant to be changed only by "power-users".
}

\value{
  Either an \code{affysnpBatch} or an \code{affyexprsBatch} object,
  depending on the type of the array (SNP or expression).
}

\author{Benilton Carvalho \email{bcarvalh@jhsph.edu}, Rafael A. Irizarry, Ben Bolstad.}

\seealso{ \code{\link{list.celfiles}}, \code{\link{read.xysfiles}}, \code{\link{list.xysfiles}} }
\examples{
## NOT RUN
##
## oligoBatch <- read.celfiles(list.celfiles())
}
\keyword{manip}% at least one, from doc/KEYWORDS
