\name{getContainer}
\Rdversion{1.1}
\alias{getContainer}
\alias{getContainer-methods}
\alias{getContainer,TilingFeatureSet-method}
\title{
  Get container information for NimbleGen Tiling Arrays.
}
\description{
  Get container information for NimbleGen Tiling Arrays. This is useful
  for better identification of control probes.
}
\usage{
getContainer(object, probeType)
}
\arguments{
  \item{object}{
    A \code{TilingFeatureSet} or \code{TilingFeatureSet2} object.
}
  \item{probeType}{
    String describing which probes to query ('pm', 'bg')
}
}
\value{
  'character' vector with container information.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
