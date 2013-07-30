\name{getProbeInfo}
\alias{getProbeInfo}
\alias{availProbeInfo}

\title{
  Probe information selector.
}
\description{
  A tool to simplify the selection of probe information, so user does
  not need to use the SQL approaches.
}
\usage{
getProbeInfo(object, field, probeType = "pm", target = "core", sortBy = c("fid", "man_fsetid", "none"), ...)
}

\arguments{
  \item{object}{\code{FeatureSet} object.}
  \item{field}{character string with names of field(s) of interest to be obtained from database.}
  \item{probeType}{character string: 'pm' or 'mm'}
  \item{target}{Used only for Exon or Gene ST arrays: 'core', 'full', 'extended', 'probeset'.}
  \item{sortBy}{Field to be used for sorting.}
  \item{...}{Arguments to be passed to \code{\link{subset}}}
}

\value{
  A \code{data.frame} with the probe level information.
}

\author{Benilton Carvalho}

\note{
  The code allows for querying info on MM probes, however it has been
  used mostly on PM probes.
}

\examples{
if (require(oligoData)){
   data(affyGeneFS)
   availProbeInfo(affyGeneFS)
   probeInfo <- getProbeInfo(affyGeneFS, c('fid', 'x', 'y', 'chrom'))
   head(probeInfo)
   ## Selecting antigenomic background probes
   agenGene <- getProbeInfo(affyGeneFS, field=c('fid', 'fsetid', 'type'), target='probeset', subset= type == 'control->bgp->antigenomic')
   head(agenGene)
}
}

\keyword{manip}
