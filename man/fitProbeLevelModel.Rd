\name{fitProbeLevelModel}
\alias{fitProbeLevelModel}

\title{
  Tool to fit Probe Level Models.
}
\description{
  Fits robust Probe Level linear Models to all the (meta)probesets
  in an \code{\link[oligoClasses:FeatureSet]{FeatureSet}}. This is carried out
  on a (meta)probeset by (meta)probeset basis.
}
\usage{
fitProbeLevelModel(object, target = "core", subset, method = "plm", verbose = TRUE, S4 = TRUE)
}

\arguments{
  \item{object}{\code{\link[oligoClasses:FeatureSet]{FeatureSet}} object.}
  \item{target}{character vector describing the summarization
    target. Valid values are: 'probeset', 'core' (Gene/Exon), 'full'
    (Exon), 'extended' (Exon).}
  \item{subset}{subset to be passed down to \code{\link{getProbeInfo}}}
  \item{method}{summarization method to be used.}
  \item{verbose}{verbosity flag.}
  \item{S4}{return final value as an S4 object (\code{oligoPLM}) if
  \code{TRUE}. If \code{FALSE}, final value is returned as a \code{list}.}
}

\value{\code{fitProbeLevelModel} returns an \code{\link{oligoPLM}}
  object, if \code{S4=TRUE}; otherwise, it will return a list.}

\note{
  This is the initial port of \code{fitPLM} to oligo. Some features
  found on the original work by Ben Bolstad (in the affyPLM package) may
  not be yet available. If you found one of this missing
  characteristics, please contact Benilton Carvalho.
}

\author{This is a simplified port from Ben Bolstad's work implemented in
  the affyPLM package. Problems with the implementation in oligo should
  be reported to Benilton Carvalho.}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  Berkeley.}

\seealso{\code{\link[oligo]{rma}}, \code{\link[oligo]{summarizationMethods}}}

\examples{
if (require(oligoData)){
  data(affyExpressionFS)
  fit <- fitProbeLevelModel(affyExpressionFS)
  image(fit)
  NUSE(fit)
  RLE(fit)
}
}

\keyword{manip}
