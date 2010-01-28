\name{summarize}
\alias{summarize}
\alias{summarize-methods}
\alias{summarize,matrix-method}
\alias{summarize,big.matrix-method}
\alias{backgroundCorrect}
\alias{backgroundCorrect-methods}
\alias{backgroundCorrect,matrix-method}
\alias{backgroundCorrect,big.matrix-method}
\alias{normalize}
\alias{normalize-methods}
\alias{normalize,matrix-method}
\alias{normalize,big.matrix-method}

\title{
  Tools for microarray preprocessing
}
\description{
  Preprocess microarray data. Includes background correction,
  normalization and summarization methods.
}
\usage{
backgroundCorrect(object, method="rma", copy=TRUE, verbose=TRUE)
summarize(object, probes=rownames(object), method="medianpolish", verbose=TRUE)
normalize(object, method="quantile", copy=TRUE, verbose=TRUE)
}
\arguments{
  \item{object}{Object containing probe intensities to be preprocessed.}
  \item{method}{String determining which method to use at that
    preprocessing step.}
  \item{probes}{Character vector that identifies the name of the probes represented
    by the rows of \code{object}.}
  \item{copy}{Logical flag determining if data must be copied before
    processing (TRUE), or if data can be overwritten (FALSE).}
  \item{verbose}{Logical flag for verbosity.}
}
\details{
  Number of rows of \code{object} must match the length of
  \code{probes}.
  Currently, only the following methods are implemented:
  - backgroundCorrectiong: rma.background
  - normalize: quantile
  - summarization: median-polish
}
\value{
  Expression matrix with \code{length(unique(probes))} rows and
  \code{ncol(object)} columns.
}

\examples{
ns <- 100
nps <- 1000
np <- 10
intensities <- matrix(rnorm(ns*nps*np, 8000, 400), nc=ns)
ids <- rep(as.character(1:nps), each=np)
bgCorrected <- backgroundCorrect(intensities)
normalized <- normalize(bgCorrected)
expression <- summarize(normalized, probes=ids)
intensities[1:20, 1:3]
expression[1:20, 1:3]
}
\keyword{manip}
