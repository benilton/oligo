\name{hist-methods}
\docType{methods}
\alias{hist-methods}
\alias{hist,FeatureSet-method}
\alias{hist,ExpressionSet-method}
\title{Density estimate}
\description{
  Plot the density estimates for each sample
}
\usage{

\S4method{hist}{FeatureSet}(x, transfo=log2, which=c("pm", "mm", "bg", "both", "all"),
                   nsample=10000, ...)

\S4method{hist}{ExpressionSet}(x, transfo=identity, nsample=10000, ...)

}
\arguments{
  \item{x}{\code{FeatureSet} or \code{ExpressionSet} object}
  \item{transfo}{a function to transform the data before plotting. See
    'Details'.}
  \item{nsample}{number of units to sample and build the plot.}
  \item{which}{set of probes to be plotted ("pm", "mm", "bg", "both",
    "all").}
  \item{\dots}{arguments to be passed to \code{matplot}}
}

\details{
  The 'transfo' argument will set the transformation to be used. For raw
  data, 'transfo=log2' is a common practice. For summarized data (which
  are often in log2-scale), no transformation is needed (therefore 'transfo=identity').
}

\note{
  The hist methods for \code{FeatureSet} and \code{Expression} use a
  sample (via \code{sample}) of the probes/probesets to produce the
  plot (unless nsample > nrow(x)). Therefore, the user interested in reproducibility is advised to
  use \code{set.seed}.
  }

\keyword{hplot}
\keyword{methods}
