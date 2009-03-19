\name{rma}
\alias{rma}
\alias{rma,ExpressionFeatureSet-method}
\alias{rma,ExonFeatureSet-method}
\alias{rma,GeneFeatureSet-method}
\title{RMA - Robust Multichip Average algorithm}
\description{
  Robust Multichip Average methodology. This will convert an
  (Expression/Exon/Gene)FeatureSet object to an ExpressionSet object by
  using RMA strategy.
}
\usage{
rma(object, ...)
}

\arguments{
  \item{object}{\code{FeatureSet} object}
  \item{\dots}{Extra arguments.}
}
\value{
  \code{ExpressionSet} object.
}

\details{
  This function computes the RMA (Robust Multichip Average) expression
  measure described in Irizarry et al
  Biostatistics (2003).

  Note that this expression measure is given to you in log base 2
  scale. This differs from most of the other expression measure methods.

  For Exon ST and Gene ST arrays, the user should be aware that the
  summarization is performed to the *probeset* level. The ExpressionSet
  returned when either Exon/Gene-FeatureSet objects are passed contain
  extra annotation on the featureData slot that the user should take
  into account for exon/gene-level analyses.
}

\references{
  Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie
  M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summa
  ries of Affymetrix GeneChip probe level data Nucleic Acids Research
  31(4):e15

  Bolstad, B.M., Irizarry R. A., Astrand M., and Speed, T.P. (2003),
  A Comparison of Normalization Methods for High Density O
  ligonucleotide Array Data Based on Bias and
  Variance. Bioinformatics 19(2):185-193
  
  Irizarry, RA, Hobbs, B, Collin, F, Beazer-Barclay, YD,
  Antonellis, KJ, Scherf, U, Speed, TP (2003) Exploration,
  Normalizati
  on, and Summaries of High Density Oligonucleotide Array Probe
  Level Data. Biostatistics .Vol. 4, Number 2: 249-264
}
      
\seealso{\code{\link{snprma}}}
\examples{
if (require(maqcExpression4plex) & require(pd.hg18.60mer.expr)){
xysPath <- system.file("extdata", package="maqcExpression4plex")
xysFiles <- list.xysfiles(xysPath, full.name=TRUE)
ngsExpressionFeatureSet <- read.xysfiles(xysFiles)
summarized <- rma(ngsExpressionFeatureSet)
show(summarized)
}
}
\keyword{manip}

