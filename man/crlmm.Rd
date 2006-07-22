\name{crlmm}
\alias{crlmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Corrected Robust Linear Model with Maximum-Likelihood-Based Distance Classifier}
\description{
  The procedure calls genotypes using CRLMM.
}
\usage{
crlmm(object, correction = NULL, recalibrate = TRUE, minLLRforCalls = c(50, 40, 50), returnCorrectedM = FALSE, returnParams = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpQSet} object;}
  \item{correction}{List with correction parameters. See Details.}
  \item{recalibrate}{Logical flag determining if recalibration should be
  performed}
  \item{minLLRforCalls}{vector, with 3 elements, specifying the min LLR
    to be used for calls}
  \item{returnCorrectedM}{Logical flag specifying if corrected
    log-ratios are to be returned.}
  \item{returnParams}{Logical flag specifying if parameters are to be returned.}
  \item{verbose}{Logical flag for verbosity.}
}
\details{
  The \code{correction} argument is a list with the following elements:
  \code{f0} (scalar), \code{fs} (numeric vector), \code{pis} (numeric
  matrix) and \code{snr}.
}
\value{
  Returns a list with two elements:
  \item{calls}{matrix with calls (1 - AA; 2 - AB; 3 - BB)}
  \item{llr}{matrix with log-likelihood ratios}
}
\references{Benilton Carvalho, Terence P. Speed, and Rafael A. Irizarry,
"EXPLORATION, NORMALIZATION, AND GENOTYPE CALLS OF HIGH DENSITY
OLIGONUCLEOTIDE SNP ARRAY DATA" (July 2006). Johns Hopkins University,
Dept. of Biostatistics Working Papers. Working Paper 111.}
\author{Benilton S. Carvalho, Rafael A. Irizarry, Shin Lin}
\seealso{\code{\link{snprma}}, \code{\link{justsnprma}}}
\examples{
## DO NOT RUN
##
## files <- list.celfiles()
## snpFSet <- read.celfiles(files)
## snpQSet <- snprma(snpFSet)
## theCalls <- crlmm(snpQSet)
}
\keyword{methods}% at least one, from doc/KEYWORDS
