\name{boxplot}
\alias{boxplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot for \code{oligoBatch}}
\description{
  Creates boxplots for \code{oligoBatch} objects.
}
\usage{
boxplot(x, which=c("both", "pm", "mm"), range=0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{oligoBatch} object.}
  \item{which}{A string setting wich probes should be plotted.}
  \item{range}{this determines how far the plot whiskers extend out from
    the box.  If 'range' is positive, the whiskers extend to the most
    extreme data point which is no more than 'range' times the
    interquartile range from the box. A value of zero causes the
    whiskers to extend to the data extremes.}
  \item{\dots}{Other arguments to be passed to \code{boxplot}}
}
\seealso{\code{\link{boxplot}}}
\keyword{hplot}% at least one, from doc/KEYWORDS

