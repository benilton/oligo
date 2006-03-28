\name{image}
\alias{image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Array Image}
\description{
  Creates an image for each sample.
}
\usage{
image(x, transfo=log, col=gray(c(0:64)/64), xlab="", ylab="", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an \code{oligoBatch} object.}
  \item{transfo}{Transformation to be applied to the data.}
  \item{col}{Colors to be used.}
  \item{ylab}{Title for the y axis.}
  \item{xlab}{Title for the x axis.}
  \item{\dots}{further arguments to be passed to \code{image}}
}
\keyword{hplot}% at least one, from doc/KEYWORDS
