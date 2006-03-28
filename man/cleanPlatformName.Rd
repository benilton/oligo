\name{cleanPlatformName}
\alias{cleanPlatformName}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Naming PDEnvs}
\description{
  The function suggests a good names for a Platform Design Environment.
}
\usage{
cleanPlatformName(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
}
\details{
  The function removes underscore and dash symbols from the string and
  adds the string "pd" at the beginning.
}
\value{
  A string in lower case.
}
\examples{
cleanPlatformName("Mapping50K_Xba240.CDF")
}
\keyword{manip}% at least one, from doc/KEYWORDS
