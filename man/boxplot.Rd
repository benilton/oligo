\name{boxplot-methods}
\docType{methods}
\alias{boxplot-methods}
\alias{boxplot,ExpressionSet-method}
\alias{boxplot,FeatureSet-method}
\title{Boxplot}
\description{
  Boxplot for observed (log-)intensities in a FeatureSet-like object
  (ExpressionFeatureSet, ExonFeatureSet, SnpFeatureSet,
  TilingFeatureSet) and ExpressionSet.
}

\usage{
\S4method{boxplot}{FeatureSet}(x, which=c("pm", "mm", "bg", "both", "all"), transfo=log2, nsample=10000, \dots)

\S4method{boxplot}{ExpressionSet}(x, which, transfo=identity, nsample=10000, \dots)
}

\arguments{
  \item{x}{a \code{FeatureSet}-like object or \code{ExpressionSet}
    object.}
  \item{which}{character defining what probe types are to be used in the
  plot.}
  \item{transfo}{a function to transform the data before plotting. See
    'Details'.}
  \item{nsample}{number of units to sample and build the plot.}
  \item{\dots}{arguments to be passed to the default boxplot method.}
}

\details{
  The 'transfo' argument will set the transformation to be used. For raw
  data, 'transfo=log2' is a common practice. For summarized data (which
  are often in log2-scale), no transformation is needed (therefore 'transfo=identity').
}

\note{
  The boxplot methods for \code{FeatureSet} and \code{Expression} use a
  sample (via \code{sample}) of the probes/probesets to produce the
  plot. Therefore, the user interested in reproducibility is advised to
  use \code{set.seed}.
  }

\seealso{\code{\link{hist}}, \code{\link{image}}, \code{\link{sample}}, \code{\link{set.seed}}}

\keyword{hplot}
\keyword{methods}
