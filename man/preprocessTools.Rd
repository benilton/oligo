\name{summarize}
\alias{backgroundCorrectionMethods}
\alias{normalizationMethods}
\alias{summarizationMethods}
\alias{summarize}
\alias{summarize-methods}
\alias{summarize,matrix-method}
\alias{summarize,ff_matrix-method}
\alias{backgroundCorrect}
\alias{backgroundCorrect-methods}
\alias{backgroundCorrect,matrix-method}
\alias{backgroundCorrect,ff_matrix-method}
\alias{backgroundCorrect,FeatureSet-method}
\alias{normalize}
\alias{normalize-methods}
\alias{normalize,matrix-method}
\alias{normalize,ff_matrix-method}
\alias{normalize,FeatureSet-method}
\alias{normalizeToTarget}
\alias{normalizeToTarget-methods}
\alias{normalizeToTarget,matrix-method}
\alias{normalizeToTarget,ff_matrix-method}

\title{
  Tools for microarray preprocessing.
}

\description{
  These are tools to preprocess microarray data. They include background correction,
  normalization and summarization methods.
}

\usage{
backgroundCorrectionMethods()
normalizationMethods()
summarizationMethods()
backgroundCorrect(object, method=backgroundCorrectionMethods(), copy=TRUE, extra, subset=NULL, target='core', verbose=TRUE)
summarize(object, probes=rownames(object), method="medianpolish", verbose=TRUE, ...)
normalize(object, method=normalizationMethods(), copy=TRUE, subset=NULL,target='core', verbose=TRUE, ...)
normalizeToTarget(object, target, method="quantile", copy=TRUE, verbose=TRUE)
}

\arguments{
  \item{object}{Object containing probe intensities to be preprocessed.}
  \item{method}{String determining which method to use at that
    preprocessing step.}
  \item{target}{Vector with the target distribution}
  \item{probes}{Character vector that identifies the name of the probes represented
    by the rows of \code{object}.}
  \item{copy}{Logical flag determining if data must be copied before
    processing (TRUE), or if data can be overwritten (FALSE).}
  \item{subset}{Not yet implemented.}
  \item{target}{One of the following values: 'core', 'full', 'extended',
  'probeset'. Used only with Gene ST and Exon ST designs.}
  \item{extra}{Extra arguments to be passed to other methods.}
  \item{verbose}{Logical flag for verbosity.}
  \item{\dots}{Arguments to be passed to methods.}
}

\details{
  Number of rows of \code{object} must match the length of
  \code{probes}.
}

\value{
  \code{backgroundCorrectionMethods} and \code{normalizationMethods}
  will return a character vector with the methods implemented currently.

  \code{backgroundCorrect}, \code{normalize} and
  \code{normalizeToTarget} will return a matrix with same dimensions as
  the input matrix. If they are applied to a FeatureSet object, the PM
  matrix will be used as input.

  The \code{summarize} method will return a matrix with
  \code{length(unique(probes))} rows and \code{ncol(object)} columns.
}

\examples{
ns <- 100
nps <- 1000
np <- 10
intensities <- matrix(rnorm(ns*nps*np, 8000, 400), nc=ns)
ids <- rep(as.character(1:nps), each=np)
bgCorrected <- backgroundCorrect(intensities)
normalized <- normalize(bgCorrected)
summarizationMethods()
expression <- summarize(normalized, probes=ids)
intensities[1:20, 1:3]
expression[1:20, 1:3]
target <- rnorm(np*nps)
normalizedToTarget <- normalizeToTarget(intensities, target)

if (require(oligoData) & require(pd.hg18.60mer.expr)){
  ## Example of normalization with real data
  data(nimbleExpressionFS)
  boxplot(nimbleExpressionFS, main='Original')
  for (mtd in normalizationMethods()){
    message('Normalizing with ', mtd)
    res <- normalize(nimbleExpressionFS, method=mtd, verbose=FALSE)
    boxplot(res, main=mtd)
  }
}

}
\keyword{manip}
