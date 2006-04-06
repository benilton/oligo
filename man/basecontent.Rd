\name{Obtain the ATCG content of a gene}
\alias{basecontent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Boxplot}
\description{
  This function accepts a character vector representing the nucleotide
  sequences and computes the frequencies of each base (A, C, G, T).
}
\usage{
basecontent(seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{character vector.}
}
\details{
  The base frequencies are calculated separately for each element of
  'x'. The elements of 'x' can be in upper case, lower case or mixed.
}
\value{
  A matrix with 4 columns and 'length(x)' rows. The columns are names
  'A', 'C', 'T', 'G', and the values in each column are the counts of
  the corresponding bases in the elements of 'x'.
}
\note{Based on MPbasecontent in matchprobes}
\examples{
seqs <- c("AATCG", "gTCa", "aaaTTTT")
basecontent(seqs)
}
\keyword{methods}% at least one, from doc/KEYWORDS

