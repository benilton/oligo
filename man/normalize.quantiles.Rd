\name{Quantile Normalization}
\alias{normalize.quantiles}
\alias{normalize.oligoBatch.quantiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Quantile Normalization Functions}
\description{
  Functions to perform quantile normalization
}
\usage{
normalize.quantiles(x, copy = TRUE)
normalize.oligoBatch.quantiles(obatch, type=c("separate","pmonly","mmonly","together"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A matrix of intensities where each column corresponds to a
    chip and each row is a probe.}
  \item{copy}{Make a copy of matrix before normalizing. Usually safer to
    work with a copy.}
  \item{obatch}{An \code{oligoBatch} object.}
  \item{type}{A string specifying how the normalization should be applied.}
}
\details{
  This method is based upon the concept of a quantile-quantile plot
  extended to n dimensions. No special allowances are made for
  outliers. If you make use of quantile normalization either through
  'rma' or 'expresso' please cite Bolstad et al, Bioinformatics (2003).

  The type arguement should be one of "separate", "pmonly", "mmonly",
  "together" which indicates whether to normalize only one probe type
  (PM,MM) or both together or separately.
}
\value{
  A normalized \code{oligoBatch}.
}
\references{
  Bolstad, B (2001) Probe Level Quantile Normalization of High Density
  Oligonucleotide Array Data. Unpublished manuscript
  \url{http://bmbolstad.com/stuff/qnorm.pdf}

  Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003) A
  Comparison of Normalization Methods for High Density Oligonucleotide
  Array Data Based on Bias and Variance. Bioinformatics 19(2) ,pp
  185-193.
  \url{http://bmbolstad.com/misc/normalize/normalize.html}
}
\keyword{methods}% at least one, from doc/KEYWORDS
