\name{boxplot}
\alias{boxplot}
\alias{boxplot,FeatureSet-method}
\alias{boxplot,ExpressionSet-method}
\title{Boxplot}
\description{
  Boxplot for observed (log-)intensities in a FeatureSet-like object
  (ExpressionFeatureSet, ExonFeatureSet, SnpFeatureSet, TilingFeatureSet).
}
\usage{
boxplot(x, \dots)
\S4method{boxplot}{FeatureSet}(x, which=pmindex(x), transfo=log2,
range=0, \dots)
\S4method{boxplot}{ExpressionSet}(x, which=1:nrow(x), transfo=identity,
range=0, \dots)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{a \code{FeatureSet}-like object or \code{ExpressionSet}
    object.}
  \item{which}{an integer vector determining which rows of 'x' should be
    plotted. See 'Details'.}
  \item{transfo}{a function to transform the data before plotting. See
    'Details'.}
  \item{range}{this determines how far the plot whiskers extend out from
    the box.  If 'range' is positive, the whiskers extend to the most
    extreme data point which is no more than 'range' times the
    interquartile range from the box. A value of zero causes the
    whiskers to extend to the data extremes.}
  \item{\dots}{arguments to be passed to plot}
}

\details{
  The 'which' argument should be used to subset the object to be
  plotted. For example, if the user wants to plot the PM probes, he
  should use 'which=pmindex(x)', if MM probes are to be plotted
  'which=mmindex(x)'. If all probes are to be plotted
  'which=1:nrow(x)'. Note that pmindex/mmindex options will not work for
  summarized data.

  The 'transfo' argument will set the transformation to be used. For raw
  data, 'transfo=log2' is a common practice. For summarized data (which
  are often in log2-scale), no transformation is needed (therefore 'transfo=identity').
  }
\seealso{\code{\link{hist}}, \code{\link{image}}}
%% Add examples
\keyword{hplot}
