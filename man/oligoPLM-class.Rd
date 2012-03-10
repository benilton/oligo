\name{oligoPLM-class}
\Rdversion{1.1}
\docType{class}
\alias{oligoPLM-class}
\alias{annotation,oligoPLM-method}
\alias{boxplot,oligoPLM-method}
\alias{coef,oligoPLM-method}
\alias{coefs.probe,oligoPLM-method}
\alias{geometry,oligoPLM-method}
\alias{image,oligoPLM-method}
\alias{manufacturer,oligoPLM-method}
\alias{method,oligoPLM-method}
\alias{ncol,oligoPLM-method}
\alias{nprobes,oligoPLM-method}
\alias{nprobesets,oligoPLM-method}
\alias{residuals,oligoPLM-method}
\alias{residualSE,oligoPLM-method}
\alias{se,oligoPLM-method}
\alias{se.probe,oligoPLM-method}
\alias{show,oligoPLM-method}
\alias{weights,oligoPLM-method}

\title{Class \code{"oligoPLM"}}
\description{
  A class to represent Probe Level Models.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{fitProbeLevelModel(FeatureSetObject)}, where
  \code{FeatureSetObject} is an object obtained through
  \code{read.celfiles} or \code{read.xysfiles}, representing intensities
  observed for different probes (which are grouped in probesets or
  meta-probesets) across distinct samples.
}
\section{Slots}{
  \describe{
    \item{\code{chip.coefs}:}{\code{"matrix"} with chip/sample effects -
    probeset-level}
    \item{\code{probe.coefs}:}{\code{"numeric"} vector with probe effects}
    \item{\code{weights}:}{\code{"matrix"} with weights - probe-level}
    \item{\code{residuals}:}{\code{"matrix"} with residuals - probe-level}
    \item{\code{se.chip.coefs}:}{\code{"matrix"} with standard errors
      for chip/sample coefficients}
    \item{\code{se.probe.coefs}:}{\code{"numeric"} vector with standard
      errors for probe effects}
    \item{\code{residualSE}:}{scale - residual standard error}
    \item{\code{geometry}:}{array geometry used for plots}
    \item{\code{method}:}{\code{"character"} string describing method
      used for PLM}
    \item{\code{manufacturer}:}{\code{"character"} string with manufacturer name}
    \item{\code{annotation}:}{\code{"character"} string with the name of
      the annotation package}
    \item{\code{narrays}:}{\code{"integer"} describing the number of arrays}
    \item{\code{nprobes}:}{\code{"integer"} describing the number of
      probes before summarization}
    \item{\code{nprobesets}:}{\code{"integer"} describing the number of
      probesets after summarization}
  }
}
\section{Methods}{
  \describe{
    \item{annotation}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to annotation slot}
    \item{boxplot}{\code{signature(x = "oligoPLM")}: boxplot method}
    \item{coef}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to coef slot}
    \item{coefs.probe}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to coefs.probe slot}
    \item{geometry}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to geometry slot}
    \item{image}{\code{signature(x = "oligoPLM")}: image method}
    \item{manufacturer}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to manufacturer slot}
    \item{method}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to method slot}
    \item{ncol}{\code{signature(x = "oligoPLM")}: accessor/replacement
      method to ncol slot}
    \item{nprobes}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to nprobes slot}
    \item{nprobesets}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to nprobesets slot}
    \item{residuals}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to residuals slot}
    \item{residualSE}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to residualSE slot}
    \item{se}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to se slot}
    \item{se.probe}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to se.probe slot}
    \item{show}{\code{signature(object = "oligoPLM")}: show method}
    \item{weights}{\code{signature(object = "oligoPLM")}:
      accessor/replacement method to weights slot}
  }
}


\author{This is a port from Ben Bolstad's work implemented in the
  affyPLM package. Problems with the implementation in oligo should be
  reported to the package's maintainer.}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  Berkeley.}

\seealso{\code{\link[oligo]{rma}}, \code{\link[oligo]{summarize}}}

\examples{
## TODO: review code and fix broken
\dontrun{
if (require(oligoData)){
  data(nimbleExpressionFS)
  fit <- fitProbeLevelModel(nimbleExpressionFS)
  image(fit)
  NUSE(fit)
  RLE(fit)
}
}
}
\keyword{classes}
