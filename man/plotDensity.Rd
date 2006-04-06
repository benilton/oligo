\name{Plot Densities}
\alias{plotDensity}
\alias{plotDensity.oligoBatch}
\alias{hist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot Densities}
\description{
  Plots the non-parametric density estimates using values contained in
  the columns of a matrix.
}
\usage{
plotDensity(mat, ylab = "density", xlab = "x", type = "l", col = 1:6,
...)
plotDensity.oligoBatch(x, col = 1:6, log=TRUE, which=c("both","pm","mm"),
ylab = "density", xlab = NULL, ...)
hist(x, which=c("both","pm","mm"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{A matrix containing the values to make densities in the
    columns.}
  \item{x}{An \code{oligoBatch} object}
  \item{log}{logical value. If \code{TRUE} the log of the intensities in the
    \code{oligoBatch} are plotted.}
  \item{ylab}{A title for the y axis.}
  \item{xlab}{A title for the x axis.}
  \item{type}{the type of the plot.}
  \item{col}{The colors to use fot the different arrays.}
  \item{which}{Feature type to be plotted.}
  \item{\dots}{graphical parameters can be given as arguments to \code{plot}}
}
\examples{
## Do not run
## oBatch <- read.celfiles(list.celfiles())
## plotDensity(exprs(oBatch))
}
\keyword{hplot}% at least one, from doc/KEYWORDS
