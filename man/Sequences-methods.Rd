\name{mmSequence}
\alias{mmSequence}
\alias{mmSequence,DBPDInfo-method}
\alias{mmSequence,FeatureSet-method}
\alias{mmSequence,AffySNPPDInfo-method}
\alias{pmSequence}
\alias{pmSequence,DBPDInfo-method}
\alias{pmSequence,FeatureSet-method}
\alias{pmSequence,AffySNPPDInfo-method}
\alias{bgSequence}
\alias{bgSequence,ExonFeatureSet-method}
\alias{bgSequence,FeatureSet-method}
\alias{bgSequence,GeneFeatureSet-method}
\alias{bgSequence,DBPDInfo-method}

\title{Probe Sequeces}
\description{
  Accessor to the (PM/MM/background) probe sequences.
}
\usage{
mmSequence(object)
pmSequence(object, ...)
bgSequence(object, ...)
}
\arguments{
  \item{object}{\code{FeatureSet}, \code{AffySNPPDInfo} or \code{DBPDInfo} object}
  \item{...}{additional arguments}
}
\value{
  A DNAStringSet containing the PM/MM/background probe sequence
  associated to the array.
}

\keyword{manip}
