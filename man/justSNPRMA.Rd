\name{justSNPRMA}
\alias{justSNPRMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarization of SNP data}
\description{
  This function implements the SNPRMA method for summarization of SNP
  data. It works directly with the CEL files, saving memory.
}
\usage{
justSNPRMA(filenames, tmpdir = getwd(), verbose = TRUE, phenoData = NULL, normalizeToHapmap = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filenames}{character vector with the filenames.}
  \item{tmpdir}{path to a directory where temporary files are going to
    be saved to. A local directory is preferred.}
  \item{verbose}{logical flag for verbosity.}
  \item{phenoData}{a \code{phenoData} object or \code{NULL}}
  \item{normalizeToHapmap}{Normalize to Hapmap? Should always be TRUE,
    but it's kept here for future use.}
}
\details{
  \code{tmpdir} is expected to be a local directory. A network directory
  will decrease the performance significantly.
}
\value{
  \code{SnpQSet} or a \code{SnpCnvQSet}, depending on the array type.
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
## snprmaResults <- justSNPRMA(list.celfiles())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

