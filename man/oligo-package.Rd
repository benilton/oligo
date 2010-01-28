\name{oligo-package}
\alias{oligo-package}
\docType{package}
\title{
The oligo package: a tool for low-level analysis of oligonucleotide arrays
}
\description{

  The \pkg{oligo} package provides tools to preprocess different
oligonucleotide arrays types: expression, tiling, SNP and exon
chips. The supported manufacturers are Affymetrix and NimbleGen.

It offers support to large datasets (when the \pkg{bigmemory} is loaded)
and can execute preprocessing tasks in parallel (if, in addition to
\pkg{bigmemory}, the \pkg{snow} package is also loaded).
}

\details{

  The package will read the raw intensity files (CEL for Affymetrix; XYS
  for NimbleGen) and allow the user to perform analyses starting at the
  feature-level.

Reading in the intensity files require the existence of data packages
that contain the chip specific information (X/Y coordinates; feature
types; sequence). These data packages packages are built using the
\pkg{pdInfoBuilder} package. 

For Affymetrix SNP arrays, users are asked to download the already built
annotation packages from BioConductor. This is because these packages
contain metadata that are not automatically created. The following
annotation packages are available:

50K Xba  - pd.mapping50kxba.240
50K Hind - pd.mapping50khind.240
250K Sty - pd.mapping250k.sty
250K Nsp - pd.mapping250k.nsp
GenomeWideSnp 5 (SNP 5.0) - pd.genomewidesnp.5
GenomeWideSnp 6 (SNP 6.0) - pd.genomewidesnp.6

For users interested in genotype calls for SNP 5.0 and 6.0 arrays, we
strongly recommend the use use the \pkg{crlmm} package, which implements
a more efficient version of CRLMM.
}
\author{
  Benilton Carvalho - \email{carvalho@bclab.org}
}

\references{
Carvalho, B.; Bengtsson, H.; Speed, T. P. & Irizarry, R. A. Exploration, Normalization, and Genotype Calls of High Density Oligonucleotide SNP Array Data. Biostatistics, 2006.
}

\keyword{ package }
