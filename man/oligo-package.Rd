\name{oligo-package}
\alias{oligo-package}
\docType{package}
\title{
The oligo package: a tool for analysis of oligonucleotide arrays
}
\description{
The \pkg{oligo} package handles oligonucleotide arrays: expression,
tiling, SNP and exon chips. The supported manufacturers are Affymetrix
and NimbleGen. The package provides tools for preprocessing.
}

\details{

  The package will read the raw intensity files (CEL for Affymetrix; XYS
  for NimbleGen) and allow the user to perform analyses starting at the
  feature-level.

  Reading in the intensity files require the existence of data packages
that contain the chip specific information (X/Y coordinates; feature
types; sequence). These data packages packages are build using the
\pkg{makePlatformDesign} package, except for the Affymetrix SNP arrays
that uses the data packages built via \pkg{pdInfoBuilder} package. Soon,
all the data packages will be built using the \pkg{pdInfoBuilder}
package.

}
\author{
  Benilton Carvalho - \email{bcarvalh@jhsph.edu}
}

\references{
Carvalho, B.; Bengtsson, H.; Speed, T. P. & Irizarry, R. A. Exploration, Normalization, and Genotype Calls of High Density Oligonucleotide SNP Array Data. Biostatistics, 2006.
}


\keyword{ package }
\examples{
## rawData <- read.celfiles(list.celfiles())
## preProc <- snprma(rawData)
## gtCalls <- crlmm(preProc, correctionFile="outEM.rda")
}
