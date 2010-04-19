\name{rma-methods}
\docType{methods}
\alias{rma}
\alias{rma-methods}
\alias{rma,ExonFeatureSet-method}
\alias{rma,ExpressionFeatureSet-method}
\alias{rma,GeneFeatureSet-method}
\alias{rma,SnpCnvFeatureSet-method}
\title{RMA - Robust Multichip Average algorithm}

\description{
  Robust Multichip Average preprocessing methodology. This strategy
  allows background subtraction, quantile normalization and
  summarization (via median-polish).
}


\section{Methods}{
\describe{

  \item{\code{signature(object = "ExonFeatureSet")}}{
    When applied to an \code{ExonFeatureSet} object, \code{rma} can
produce summaries at different levels: probeset (as defined in the PGF),
core genes (as defined in the core.mps file), full genes (as defined in
the full.mps file) or extended genes (as defined in the extended.mps
file). To determine the level for summarization, use the \code{target}
argument.
}

\item{\code{signature(object = "ExpressionFeatureSet")}}{
When used on an \code{ExpressionFeatureSet} object, \code{rma} produces
summaries at the probeset level (as defined in the CDF or NDF files,
depending on the manufacturer).
}

\item{\code{signature(object = "GeneFeatureSet")}}{
When applied to a \code{GeneFeatureSet} object, \code{rma} can produce
summaries at different levels: probeset (as defined in the PGF) and
'core genes' (as defined in the core.mps file). To determine the level for
summarization, use the \code{target} argument.
}

\item{\code{signature(object = "SnpCnvFeatureSet")}}{
If used on a \code{SnpCnvFeatureSet} object (ie., SNP 5.0 or SNP 6.0
arrays), \code{rma} will produce summaries for the CNV probes. Note that
this is an experimental feature for internal (and quick) assessment of
CNV probes. We recommend the use of the 'crlmm' package, which contains
a Copy Number tool specifically designed for these data.
}

}

}

\usage{
  \S4method{rma}{ExonFeatureSet}(object, background=TRUE, normalize=TRUE, subset=NULL, target="core")
  \S4method{rma}{ExpressionFeatureSet}(object, background=TRUE, normalize=TRUE, subset=NULL)
  \S4method{rma}{GeneFeatureSet}(object, background=TRUE, normalize=TRUE, subset=NULL, target="core")
  \S4method{rma}{SnpCnvFeatureSet}(object, background=TRUE, normalize=TRUE, subset=NULL)
}

\references{
  
  Rafael. A. Irizarry, Benjamin M. Bolstad, Francois Collin, Leslie
  M. Cope, Bridget Hobbs and Terence P. Speed (2003), Summaries of
  Affymetrix GeneChip probe level data Nucleic Acids Research 31(4):e15

  Bolstad, B.M., Irizarry R. A., Astrand M., and Speed, T.P. (2003), A
  Comparison of Normalization Methods for High Density O ligonucleotide
  Array Data Based on Bias and Variance. Bioinformatics 19(2):185-193
  
  Irizarry, RA, Hobbs, B, Collin, F, Beazer-Barclay, YD, Antonellis, KJ,
  Scherf, U, Speed, TP (2003) Exploration, Normalizati on, and Summaries
  of High Density Oligonucleotide Array Probe Level
  Data. Biostatistics. Vol. 4, Number 2: 249-264
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

\keyword{methods}
