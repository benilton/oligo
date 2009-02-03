\name{justSNPRMA}
\alias{justSNPRMA}
\title{Summarization of SNP data}
\description{
  This function implements the SNPRMA method for summarization of SNP
  data. It works directly with the CEL files, saving memory.
}
\usage{
justSNPRMA(filenames, verbose = TRUE, phenoData = NULL, normalizeToHapmap = TRUE)
}
\arguments{
  \item{filenames}{character vector with the filenames.}
  \item{verbose}{logical flag for verbosity.}
  \item{phenoData}{a \code{phenoData} object or \code{NULL}}
  \item{normalizeToHapmap}{Normalize to Hapmap? Should always be TRUE,
    but it's kept here for future use.}
}
\value{
  \code{SnpQSet} or a \code{SnpCnvQSet}, depending on the array type.
}
%% Add Examples
\examples{
## snprmaResults <- justSNPRMA(list.celfiles())
}
\keyword{manip}

