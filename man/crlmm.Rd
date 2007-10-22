\name{crlmm}
\alias{crlmm}
\alias{justCRLMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genotype Calls}
\description{
  Performs genotype calls via CRLMM (Corrected Robust Linear Model with
  Maximum-likelihood based distances).
}
\usage{
crlmm(object, correction = NULL, recalibrate = TRUE, minLLRforCalls =
c(5, 1, 5), verbose = TRUE,  correctionFile = NULL, prefix =
"tmp.crlmm.", balance = 1.5)
justCRLMM(filenames, batch_size = 40000, minLLRforCalls = c(5, 1, 5),
recalibrate = TRUE, balance = 1.5, phenoData = NULL, verbose = TRUE,
pkgname = NULL, tmpdir=tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpQSet} object}
  \item{filenames}{character vector with the filenames.}
  \item{batch_size}{integer defining how many SNPs should be processed
    at a time.}
  \item{correction}{The output of the EM algorithm}
  \item{recalibrate}{Logical - should recalibration be performed?}
  \item{prefix}{String defining the prefix to be used with the temporary files.}
  \item{balance}{Control parameter to balance homozygotes and
    heterozygotes calls.}
  \item{minLLRforCalls}{Minimum thresholds for genotype calls.}
  \item{verbose}{Logical.}
  \item{correctionFile}{A filename.}
  \item{phenoData}{\code{phenoData} object or \code{NULL}}
  \item{pkgname}{alt. pdInfo package to be used}
  \item{tmpdir}{Directory where temporary files are going to be stored at.}
}

\details{The \code{correctionFile} is a string (eg, "outputEM.rda")
  pointing to a filename. If the file does not exist, \code{crlmm} will
  save a file with that name containing the results of the EM
  algorithm. If the file exists, its content is loaded and used by
  CRLMM. The \code{correctionFile} is meant to save time if \code{crlmm}
  is be run multiple times, as the EM algorithm does not need to be run
  everytime.

  The \code{justCRLMM} method is more efficient in terms of memory. It
  uses the CEL files directly, taking SNPs by batch.
}

\value{
  \code{SnpCallSetPlus} object.
}

\references{}
\author{}
\note{}

\seealso{}
\examples{
## crlmmResults <- justCRLMM(list.celfiles())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

