\name{image-methods}
\docType{methods}
\alias{image-methods}
\alias{image,FeatureSet-method}
\title{Display a pseudo-image of a microarray chip}
\description{
  Produces a pseudo-image (\code{graphics::image}) for each sample.
}
\usage{

\S4method{image}{FeatureSet}(x, which, transfo=log2, ...)

}
\arguments{
  \item{x}{\code{FeatureSet} object}
  \item{which}{integer indices of samples to be plotted (optional).}
  \item{transfo}{function to be applied to the data prior to plotting.}
  \item{\dots}{parameters to be passed to \code{image}}
}
% Add examples
\keyword{hplot}
