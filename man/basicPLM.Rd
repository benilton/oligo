\name{basicPLM}
\alias{basicPLM}

\title{
  Simplified interface to PLM.
}
\description{
  Simplified interface to PLM.
}
\usage{
basicPLM(pmMat, pnVec, normalize = TRUE, background = TRUE, transfo =
  log2, method = c('plm', 'plmr', 'plmrr', 'plmrc'), verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pmMat}{Matrix of intensities to be processed.}
  \item{pnVec}{Probeset names}
  \item{normalize}{Logical flag: normalize?}
  \item{background}{Logical flag: background adjustment?}
  \item{transfo}{function: function to be used for data transformation
  prior to summarization.}
  \item{method}{Name of the method to be used for normalization. 'plm'
  is the usual PLM model; 'plmr' is the (row and column) robust version of PLM; 'plmrr'
  is the row-robust version of PLM; 'plmrc' is the column-robust version
  of PLM.}
  \item{verbose}{Logical flag: verbose.}
}
\value{
  A list with the following components:
  \item{Estimates}{A (length(pnVec) x ncol(pmMat)) matrix with probeset summaries.}
  \item{StdErrors}{A (length(pnVec) x ncol(pmMat)) matrix with standard errors of 'Estimates'.}
  \item{Residuals}{A (nrow(pmMat) x ncol(pmMat)) matrix of residuals.}
}
\author{
  Benilton Carvalho
}
\note{
  Currently, only RMA-bg-correction and quantile normalization are allowed.
}

\seealso{
  \code{\link[preprocessCore]{rcModelPLM}},
  \code{\link[preprocessCore]{rcModelPLMr}},
  \code{\link[preprocessCore]{rcModelPLMrr}},
  \code{\link[preprocessCore]{rcModelPLMrc}},
  \code{\link{basicRMA}}
}
\examples{
set.seed(1)
pms <- 2^matrix(rnorm(1000), nc=20)
colnames(pms) <- paste("sample", 1:20, sep="")
pns <- rep(letters[1:10], each=5)
res <- basicPLM(pms, pns, TRUE, TRUE)
res[['Estimates']][, 1:3]
res[['StdErrors']][, 1:3]
res[['Residuals']][, 1:15]
}
\keyword{manip}

