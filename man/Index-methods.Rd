\name{mmindex}
\alias{mmindex}
\alias{mmindex,FeatureSet-method}
\alias{mmindex,DBPDInfo-method}
\alias{pmindex}
\alias{pmindex,FeatureSet-method}
\alias{pmindex,DBPDInfo-method}
\alias{bgindex}
\alias{bgindex,FeatureSet-method}
\alias{bgindex,DBPDInfo-method}
%% TODO: add mmindex / bgindex?
\alias{pmindex,stArrayDBPDInfo-method}

%% This help page is ready.

\title{Accessors for PM, MM or background probes indices.}
\description{
  Extracts the indexes for PM, MM or background probes.
}
\usage{
mmindex(object, ...)
pmindex(object, ...)
bgindex(object, ...)
}

\arguments{
  \item{object}{\code{FeatureSet} or \code{DBPDInfo} object}
  \item{...}{Extra arguments, not yet implemented}
}

\details{
  The indices are ordered by 'fid', i.e. they follow the order that the
  probes appear in the CEL/XYS files.
}

\value{
  A vector of integers representing the rows of the intensity matrix
  that correspond to PM, MM or background probes.
}

\examples{
## How pm() works
\dontrun{
x <- read.celfiles(list.celfiles())
pms0 <- pm(x)
pmi <- pmindex(x)
pms1 <- exprs(x)[pmi,]
identical(pms0, pms1)
}
}
\keyword{manip}
