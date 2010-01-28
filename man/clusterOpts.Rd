\name{setCluster}
\alias{setCluster}
\alias{getCluster}
\alias{delCluster}
\alias{oligoSamples}
\alias{oligoProbesets}

\title{
  Cluster and large dataset management utilities.
}
\description{
  Tools to simplify management of clusters via 'snow' package and large
  dataset handling through the 'bigmemory' package.
}
\usage{
setCluster(...)
getCluster()
delCluster()
oligoSamples(n)
oligoProbesets(n)
}

\arguments{
  \item{\dots}{arguments to be passed to \code{makeCluster} in the
    'snow' package.}
  \item{n}{integer representing the maximum number of samples/probesets
    to be processed simultaneously on a compute node.}
}

\details{
  Some methods in the oligo package, like \code{backgroundCorrect},
  \code{normalize}, \code{summarize} and \code{rma} can use a cluster
  (set through 'snow' package). The use of cluster features is
  conditioned on the availability of the 'bigmemory' (used to
  provide shared objects across compute nodes) and 'snow' packages.

  To use a cluster, 'oligo' checks for three requirements: 1)
  'bigmemory' is loaded; 2) 'snow' is loaded; and 3) the 'cluster'
  option is set (e.g., via options(cluster=makeCluster(...)) or
  setCluster(...)).

  If only the 'bigmemory' package is available and loaded (in addition
  to 'oligo'), these methods will allow the user to analyze datasets
  that would not fit in RAM at the expense of performance.

  In the situations above (large datasets and cluster), 'oligo' uses the
  options \code{oligoSamples} and \code{oligoProbesets} to limit the
  amount of RAM used by the machine(s). For example, if oligoSamples is
  set to 100, steps like background correction and normalization process
  (in RAM) 100 samples simultaneously on each compute node. If
  oligoProbesets is set to 10K, then summarization processes 10K
  probesets at a time on each machine.
}

\section{Warning}{
  In both scenarios (large dataset and/or cluster use), there
  is a penalty in performance because data are written to disk (to
  either minimize memory footprint or share data across compute nodes).
}

\author{
  Benilton Carvalho <carvalho@bclab.org>
}
\keyword{manip}

