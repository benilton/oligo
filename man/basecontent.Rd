\name{basecontent}
\alias{basecontent}
\title{Sequence Base Contents}
\description{
  Function to compute the amounts of each nucleotide in a sequence.
}
\usage{
basecontent(seq)
}
\arguments{
  \item{seq}{character vector of length \code{n} containg a valid
    sequence (A/T/C/G)}
}
\value{
  \code{matrix} with \code{n} rows and 4 columns with the counts for
  each base.
}
\examples{
sequences <- c("ATATATCCCCG", "TTTCCGAGC")
basecontent(sequences)
}
\keyword{manip}
