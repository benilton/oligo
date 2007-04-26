\name{crlmm}
\alias{crlmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Genotype Calls}
\description{
  Performs genotype calls via CRLMM (Corrected Robust Linear Model with
  Maximum-likelihood based distances).
}
\usage{
crlmm(object, correction = NULL, recalibrate = TRUE, minPforCalls = c(.5, .1, .5), verbose = TRUE, correctionFile = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{\code{SnpQSet} object}
  \item{correction}{The output of the EM algorithm}
  \item{recalibrate}{Logical - should recalibration be performed?}
  \item{minPforCalls}{Minimum thresholds for genotype calls.}
  \item{verbose}{Logical.}
  \item{correctionFile}{A filename.}
}

\details{The \code{correctionFile} is a string (eg, "outputEM.rda")
  pointing to a filename. If the file does not exist, \code{crlmm} will
  save a file with that name containing the results of the EM
  algorithm. If the file exists, its content is loaded and used by
  CRLMM. The \code{correctionFile} is meant to save time if \code{crlmm}
  is be run multiple times, as the EM algorithm does not need to be run
  everytime.}

\value{
  \code{SnpCallSet} object.
}

\references{}
\author{}
\note{}

\seealso{}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classif}

