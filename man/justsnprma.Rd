\name{justsnprma}
\alias{justsnprma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Robust multi-array average for SNP arrays}
\description{
  Performs robust multi-array average (RMA) method for SNP Arrays. This
  implementation is recommended for a large number of CEL files and/or
  when memory is an issue.
}
\usage{
justsnprma(files)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{list of CEL files.}
}
\details{
  This implementation performs \code{snprma} on a 1-by-1 basis, trying
  to avoid problems with memory.
}
\value{
  A \code{SnpQSet} object
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{\code{\link{snprma}}}
\examples{
## Do NOT RUN!
## snpCELfiles <- list.celfiles()
## snpQSet <- justsnprma(snpCELfiles)

##---- Should be DIRECTLY executable !! ----
}
\keyword{methods}% at least one, from doc/KEYWORDS
