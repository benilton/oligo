\name{getCrlmmSummaries}
\alias{getCrlmmSummaries}
\title{Function to get CRLMM summaries saved to disk}
\description{
  This will read the summaries written to disk and return them to the
  user as a \code{SnpCallSetPlus} or \code{SnpCnvCallSetPlus} object.
}
\usage{
getCrlmmSummaries(tmpdir)
}
\arguments{
  \item{tmpdir}{directory where CRLMM saved the results to.}
}
\value{
  If the data were from SNP 5.0 or 6.0 arrays, the function will return
  a \code{SnpCnvCallSetPlus} object. It will return a \code{SnpCallSetPlus}
  object, otherwise.
}
%% Add examples
\keyword{manip}
