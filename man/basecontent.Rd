\name{basecontent}
\alias{basecontent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sequence Base Contents}
\description{
  Function to compute the amounts of each nucleotide in a sequence.
}
\usage{
basecontent(seq)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seq}{character vector of length \code{n} containg a valid
    sequence (A/T/C/G)}
}
\details{
}
\value{
  \code{matrix} with \code{n} rows and 4 columns with the counts for
  each base.
}
\references{}
\author{}
\note{
}
\seealso{}
\examples{
sequences <- c("ATATATCCCCG", "TTTCCGAGC")
basecontent(sequences)
}

\keyword{manip}
