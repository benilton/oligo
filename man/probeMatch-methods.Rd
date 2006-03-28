\name{probeMatch-methods}
\alias{pm}
\alias{mm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Accessors for perfect match (pm) and mismatch (mm) probes.}
\description{
  Methods for perfect match and mismatch probes.
}
\usage{
pm(object)
pmindex(object)
mm(object)
mmindex(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{oligoBatch object}
}
\value{
  \code{pm} and \code{mm} will return respectively the PM and MM
  matrices associated to the oligoBatch object.

  \code{pmindex} and \code{mmindex} will return vectors containing the
  indexes for the PM and MM probes.
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## Do not run
## oBatch <- read.celfiles(list.celfiles())
## pms <- pm(oBatch)
}
\keyword{methods}% at least one, from doc/KEYWORDS
