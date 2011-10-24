\name{PLMset-class}
\docType{class}
\alias{PLMset-class}
\alias{PLMset}
\alias{show,PLMset-method}
\alias{NUSE}
\alias{RLE}
\alias{NUSE,PLMset-method}
\alias{RLE,PLMset-method}
\alias{coefs}
\alias{se}
\alias{coefs,PLMset-method}
\alias{se,PLMset-method}

\title{Class PLMset}
\description{This is a class representation for Probe level Linear
  Models fitted to Affymetrix GeneChip probe level data.}
\section{Objects from the Class}{
  Objects can be created using the function \code{\link{fitPLM}}}
\section{Slots}{
  \describe{
    \item{\code{probe.coefs}:}{Object of class "matrix". Contains model
      coefficients related to probe effects.}
    \item{\code{se.probe.coefs}:}{Object of class "matrix". Contains
      standard error estimates for the probe coefficients.}  
    \item{\code{chip.coefs}:}{Object of class "matrix". Contains model
      coefficients related to chip (or chip level) effects for each fit.}
    \item{\code{se.chip.coefs}:}{Object of class "matrix". Contains
      standard error estimates for the chip coefficients.}
    \item{\code{const.coefs}:}{Object of class "matrix". Contains model
      coefficients related to intercept effects for each fit.}
    \item{\code{se.const.coefs}:}{Object of class "matrix". Contains
      standard error estimates for the intercept estimates}
    \item{\code{model.description}:}{Object of class "character". This
      string describes the probe level model fitted.}
    \item{\code{weights}:}{List of objects of class "matrix". Contains probe
      weights for each fit. The matrix has columns for chips and rows
      are probes.}
    \item{\code{phenoData}:}{Object of class "phenoData" This is an
      instance of class \code{phenoData} containing the patient
      (or case) level data. The columns of the pData slot of this
      entity represent variables and the rows represent patients or cases.}
    \item{\code{annotation}}{A character string identifying the
      annotation that may be used for the \code{ExpressionSet} instance.}
    \item{\code{experimentData}:}{Object of class "MIAME". For
      compatibility with previous version of this class description can
      also be a "character". The class \code{characterOrMIAME} has been
      defined just for this.}
    \item{\code{cdfName}:}{A character string giving the name of the
      cdfFile.}
    \item{\code{nrow}:}{Object of class "numeric". Number of rows in chip.}
    \item{\code{ncol}:}{Object of class "numeric". Number of cols in
chip.}
      \item{\code{narrays}:}{Object of class "numeric". Number of arrays
	used in model fit.}
      \item{\code{normVec}:}{Object of class "matrix". For storing
	normalization vector(s). Not currentl used}
     \item{\code{varcov}:}{Object of class "list". A list of
      variance/covariance matrices.}
    \item{\code{residualSE}:}{Object of class "matrix". Contains residual
      standard error and df.}
    \item{\code{residuals}:}{List of objects of class "matrix". Contains
      residuals from model fit (if stored).}
    \item{\code{model.call}:}{Object of class "call"}
    \item{\code{manufacturer}:}{Manufacturer string (affymetrix/nimblegen)}
  }
}
\section{Methods}{
  \describe{
    \item{boxplot}{\code{signature(x = "PLMset")}: Boxplot of Normalized
      Unscaled Standard Errors (NUSE).}
    \item{NUSE}{\code{signature(x = "PLMset")} : Boxplot of Normalized
      Unscaled Standard Errors (NUSE) or NUSE values.}
    \item{RLE}{\code{signature(x = "PLMset")} : Relative Log Expression
      boxplot or values.}
    \item{coefs}{\code{signature(object = "PLMset")}: Extractor for chip
      coefficients}
    \item{se}{\code{signature(object = "PLMset")}: Extractor for the
	standard errors for chip coefficients}
 }
}
\note{This class is better described in the vignette.}
\author{B. M. Bolstad \email{bmb@bmbolstad.com}}
\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  Berkeley.}
\keyword{classes}
