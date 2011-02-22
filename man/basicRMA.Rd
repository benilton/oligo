\name{basicRMA}
\Rdversion{1.1}
\alias{basicRMA}
\title{
  Simplified interface to RMA.
}
\description{
  Simple interface to RMA.
}
\usage{
basicRMA(pmMat, pnVec, normalize = TRUE, background = TRUE, bgversion = 2, destructive = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{pmMat}{
    Matrix of intensities to be processed.
}
  \item{pnVec}{
    Probeset names.
}
  \item{normalize}{
    Logical flag: normalize?
}
  \item{background}{
    Logical flag: background adjustment?
}
  \item{bgversion}{
    Version of background correction.
}
  \item{destructive}{
    Logical flag: use destructive methods?
}
  \item{verbose}{
    Logical flag: verbose.
}
  \item{\dots}{
    Not currently used.
}
}
\value{
  Matrix.
}
\examples{
set.seed(1)
pms <- 2^matrix(rnorm(1000), nc=20)
colnames(pms) <- paste("sample", 1:20, sep="")
pns <- rep(letters[1:10], each=5)
res <- basicRMA(pms, pns, TRUE, TRUE)
res[, 1:3]
}
\keyword{manip}
