\name{fitPLM}
\alias{fitPLM}
\title{Fit a Probe Level Model to Affymetrix Genechip Data.}
\description{
  This function converts an
  \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}} into an 
  \code{\link{PLMset}} by fitting a specified robust linear model to the
  probe level data.
}
\usage{
fitPLM(object,model=PM ~ -1 + probes +samples,
       variable.type=c(default="factor"),
       constraint.type=c(default="contr.treatment"),
       subset=NULL,
       background=TRUE, normalize=TRUE, background.method="RMA.2",
       normalize.method="quantile", background.param=list(),
       normalize.param=list(), output.param=NULL,
       model.param=NULL,
       verbosity.level=0, target="core")
}

\arguments{
  \item{object}{an \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}} object}
  \item{model}{A formula describing the model to fit. This is slightly
    different from the standard method of specifying formulae in R. Read
    the description below}
  \item{variable.type}{a way to specify whether variables in the model
    are factors or standard variables}
  \item{constraint.type}{should factor variables sum to zero or have
    first variable set to zero (endpoint constraint)}
  \item{subset}{a vector with the names of probesets to be used. If NULL
  then all probesets are used.}
  \item{normalize}{logical value. If \code{TRUE} normalize data using
    quantile normalization}
  \item{background}{logical value. If \code{TRUE} background correct
    using RMA background correction}
  \item{background.method}{name of background method to use.}
  \item{normalize.method}{name of normalization method to use.}
  \item{background.param}{A list of parameters for background routines}
  \item{normalize.param}{A list of parameters for normalization
    routines}
  \item{output.param}{A list of parameters controlling optional output
    from the routine.}
  \item{model.param}{A list of parameters controlling model procedure}
  \item{verbosity.level}{An integer specifying how much to print
    out. Higher values indicate more verbose. A value of 0 will print
    nothing}
  \item{target}{Level of summarization to be used when handling ST array
  (Exon/Gene). Allowed values for such array types are: \code{core} (Exon/Gene),
  \code{full} (Exon), \code{extended} (Exon), \code{probeset} (Exon/Gene)}
}
\details{
  This function fits robust Probe Level linear Models to all the probesets in
  an \code{\link[oligoClasses:FeatureSet-class]{FeatureSet}}. This is carried out
  on a probeset by probeset basis. The user has quite a lot of control
  over which model is used and what outputs are stored. For more details
  please read the vignette and the original documentation in the affyPLM
  package.
}

\note{
  This is the initial port of \code{fitPLM} to oligo. Some features
  found on the original work by Ben Bolstad (in the affyPLM package) may
  not be yet available. If you found one of this missing
  characteristics, please contact Benilton Carvalho.
}

\value{
  An \code{\link{PLMset}}
}

\author{This is a port from Ben Bolstad's work implemented in the
  affyPLM package. Problems with the implementation in oligo should be
  reported to Benilton Carvalho <Benilton.Carvalho@cancer.org.uk>}

\references{Bolstad, BM (2004) \emph{Low Level Analysis of High-density
    Oligonucleotide Array Data: Background, Normalization and
    Summarization}. PhD Dissertation. University of California,
  Berkeley.}

\seealso{\code{\link[oligo]{rma}}, \code{\link[oligo]{summarize}}}

\examples{
if (require(oligoData)){
  data(nimbleExpressionFS)
  fit <- fitPLM(nimbleExpressionFS)
  image(fit)
  NUSE(fit)
  RLE(fit)
}
}

\keyword{manip}
