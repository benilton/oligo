\name{getX}
\alias{getX}
\alias{getY}
\alias{getX-methods}
\alias{getX,FeatureSet-method}
\alias{getX,DBPDInfo-method}
\alias{getY-methods}
\alias{getY,FeatureSet-method}
\alias{getY,DBPDInfo-method}
\title{Accessors for physical array coordinates.}

\description{
  Accessors for physical array coordinates.
}

\usage{
getX(object, type)
getY(object, type)
}

\arguments{
  \item{object}{\code{FeatureSet} object}
  \item{type}{'character' defining the type of the probes to be
    queried. Valid options are 'pm', 'mm', 'bg'}
}
\value{
  A vector with the requested coordinates.
}
\examples{
\dontrun{
x <- read.celfiles(list.celfiles())
theXpm <- getX(x, "pm")
theYpm <- getY(x, "pm")
}
}
\keyword{manip}
