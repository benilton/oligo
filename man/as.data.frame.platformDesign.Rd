\name{as.data.frame.platformDesign}
\alias{as.data.frame.platformDesign}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Coerce to a Data Frame}
\description{
  Coerces a Platform Design Environment (PDEnv) to a data frame.
}
\usage{
  as.data.frame.platformDesign(x, row.names = NULL, optional = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a Platform Design Environment}
  \item{row.names}{'NULL' or a character vector giving the row names for
    the data frame.  Missing values are not allowed.}
  \item{optional}{logical. If 'TRUE', setting row names and converting
    column names (to syntactic names) is optional.}
}

\details{
  Behaves like \code{as.data.frame}.
}

\value{
  A data frame containing the variables available in the PDEnv.
}

\author{Benilton Carvalho \email{bcarvalh@jhsph.edu}, Rafael A. Irizarry, Ben Bolstad.}

\seealso{ \code{\link{getPD}} }
\examples{
## NOT RUN
##
## FeatureSet <- read.celfiles(list.celfiles())
## theDataFrame <- as.data.frame(getPD(FeatureSet))
}
\keyword{manip}% at least one, from doc/KEYWORDS
