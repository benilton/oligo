\name{MAplot}
\alias{MAplot}
\alias{MAplot-methods}
\alias{MAplot,FeatureSet-method}
\alias{MAplot,TilingFeatureSet-method}
\alias{MAplot,PLMset-method}
\alias{MAplot,ExpressionSet-method}
\alias{MAplot,matrix-method}
\title{MA plots}
\description{
  Create MA plots using a reference array (if one channel) or using
  channel2 as reference (if two channel).
}

\usage{
MAplot(object, ...)

\S4method{MAplot}{FeatureSet}(object, what=pm, transfo=log2, groups,
       refSamples, which, pch=".", summaryFun=rowMedians,
       plotFun=smoothScatter, main="vs pseudo-median reference chip",
       pairs=FALSE, ...)

\S4method{MAplot}{TilingFeatureSet}(object, what=pm, transfo=log2, groups,
       refSamples, which, pch=".", summaryFun=rowMedians,
       plotFun=smoothScatter, main="vs pseudo-median reference chip",
       pairs=FALSE, ...)

\S4method{MAplot}{PLMset}(object, what=coefs, transfo=identity, groups,
       refSamples, which, pch=".", summaryFun=rowMedians,
       plotFun=smoothScatter, main="vs pseudo-median reference chip",
       pairs=FALSE, ...)

\S4method{MAplot}{matrix}(object, what=identity, transfo=identity,
       groups, refSamples, which, pch=".", summaryFun=rowMedians,
       plotFun=smoothScatter, main="vs pseudo-median reference chip",
       pairs=FALSE, ...)

\S4method{MAplot}{ExpressionSet}(object, what=exprs, transfo=identity,
       groups, refSamples, which, pch=".", summaryFun=rowMedians,
       plotFun=smoothScatter, main="vs pseudo-median reference chip",
       pairs=FALSE, ...)

}
\arguments{
  \item{object}{\code{FeatureSet}, \code{PLMset} or \code{ExpressionSet}
    object.}
  \item{what}{function to be applied on \code{object} that will extract
    the statistics of interest, from which log-ratios and average
    log-intensities will be computed.}
  \item{transfo}{function to transform the data prior to plotting.}
  \item{groups}{factor describing groups of samples that will be
    combined prior to plotting. If missing, MvA plots are done per sample.}
  \item{refSamples}{integers (indexing samples) to define which subjects
    will be used to compute the reference set. If missing, a
    pseudo-reference chip is estimated using \code{summaryFun}.}
  \item{which}{integer (indexing samples) describing which samples are
    to be plotted.}
  \item{pch}{same as \code{pch} in \code{plot}}
  \item{summaryFun}{function that operates on a matrix and returns a
    vector that will be used to summarize data belonging to the same
    group (or reference) on the computation of grouped-stats.}
  \item{plotFun}{function to be used for plotting. Usually
    \code{smoothScatter}, \code{plot} or \code{points}.}
  \item{main}{string to be used in title.}
  \item{pairs}{logical flag to determine if a matrix of MvA plots is to
    be generated}
  \item{...}{Other arguments to be passed downstream, like \code{plot} arguments.}
}
\details{
  MAplot will take the following extra arguments:
  \enumerate{
    \item \code{subset}: indices of elements to be plotted to reduce
    impact of plotting 100's thousands points (if pairs=FALSE only);
    \item \code{span}: see \code{\link{loess}};
    \item \code{family.loess}: see \code{\link{loess}};
    \item \code{addLoess}: logical flag (default TRUE) to add a loess
    estimate;
    \item \code{parParams}: list of params to be passed to par() (if
    pairs=TRUE only);
  }
}
\value{Plot}
\author{Benilton Carvalho - heavily based on Ben Bolstad's original
  MAplot function.}
\seealso{\code{\link{plot}}, \code{\link{smoothScatter}}}
\examples{
if(require(oligoData) & require(pd.hg18.60mer.expr)){
  data(nimbleExpressionFS)
  nimbleExpressionFS
  groups <- factor(rep(c('brain', 'UnivRef'), each=3))
  data.frame(sampleNames(nimbleExpressionFS), groups)
  MAplot(nimbleExpressionFS, pairs=TRUE, ylim=c(-.5, .5), groups=groups)
}

}
\keyword{methods}
\keyword{hplot}
\keyword{smooth}
\keyword{loess}



