\name{getBaseProfile}
\alias{getBaseProfile}

\title{
Compute and plot nucleotide profile.
}
\description{
Computes and, optionally, lots nucleotide profile, describing the sequence effect on intensities.
}
\usage{
getBaseProfile(coefs, probeLength = 25, plot = FALSE, ...)
}

\arguments{
  \item{coefs}{affinity spline coefficients.}
  \item{probeLength}{length of probes}
  \item{plot}{logical. Plots profile?}
  \item{\dots}{arguments to be passed to matplot.}
}
\value{
  Invisibly returns a matrix with estimated effects.
}

\keyword{manip}

