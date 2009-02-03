\name{sequenceDesignMatrix}
\alias{sequenceDesignMatrix}
\title{Create design matrix for sequences}
\description{
  Creates design matrix for sequences.
}
\usage{
sequenceDesignMatrix(seqs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seqs}{\code{character} vector of 25-mers.}
}
\details{
  This assumes all sequences are 25bp long.

  The design matrix is often used when the objecive is to adjust
  intensities by sequence.
}
\value{
  Matrix with length(seqs) rows and 75 columns.
}
\examples{
genSequence <- function(x)
    paste(sample(c("A", "T", "C", "G"), 25, rep=TRUE), collapse="", sep="")
seqs <- sapply(1:10, genSequence)
X <- sequenceDesignMatrix(seqs)
Y <- rnorm(10, mean=12, sd=2)
Ydemean <- Y-mean(Y)
X[1:10, 1:3]
fit <- lm(Ydemean~X)
coef(fit)
}
\keyword{manip}
