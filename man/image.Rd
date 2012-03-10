\name{image}
\alias{image}
\docType{methods}
\alias{image-methods}
\alias{image,FeatureSet-method}
\alias{image,PLMset-method}
\title{Display a pseudo-image of a microarray chip}
\description{
  Produces a pseudo-image (\code{graphics::image}) for each sample.
}
\usage{

\S4method{image}{FeatureSet}(x, which, transfo=log2, ...)

\S4method{image}{PLMset}(x, which=0,
                   type=c("weights","resids", "pos.resids","neg.resids","sign.resids"),
                   use.log=TRUE, add.legend=FALSE, standardize=FALSE,
                   col=NULL, main, ...)

}
\arguments{
  \item{x}{\code{FeatureSet} object}
  \item{which}{integer indices of samples to be plotted (optional).}
  \item{transfo}{function to be applied to the data prior to plotting.}
  \item{type}{Type of statistics to be used.}
  \item{use.log}{Use log.}
  \item{add.legend}{Add legend.}
  \item{standardize}{Standardize residuals.}
  \item{col}{Colors to be used.}
  \item{main}{Main title.}
  \item{\dots}{parameters to be passed to \code{image}}
}

\examples{
if(require(oligoData) & require(pd.hg18.60mer.expr)){
  data(nimbleExpressionFS)
  par(mfrow=c(1, 2))
  image(nimbleExpressionFS, which=4)
##  fit <- fitPLM(nimbleExpressionFS)
##  image(fit, which=4)
  plot(1) ## while fixing fitPLM TODO
}
}

\keyword{hplot}
