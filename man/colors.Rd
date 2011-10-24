\name{darkColors}
\alias{darkColors}
\alias{seqColors}
\alias{seqColors2}
\alias{divColors}

\title{
Create set of colors, interpolating through a set of preferred colors.
}
\description{
Create set of colors, interpolating through a set of preferred colors.
}
\usage{
darkColors(n)
seqColors(n)
seqColors2(n)
divColors(n)
}

\arguments{
  \item{n}{integer determining number of colors to be generated}
}

\details{
\code{darkColors} is based on the Dark2 palette in RColorBrewer, therefore
useful to describe qualitative features of the data.

\code{seqColors} is based on Blues and generates a gradient of blues, therefore
useful to describe quantitative features of the data. \code{seqColors2}
behaves similarly, but it is based on OrRd (white-orange-red).

\code{divColors} is based on the RdBu pallete in RColorBrewer, therefore
useful to describe quantitative features ranging on two extremes.
}

\examples{
x <- 1:10
y <- 1:10
cols1 <- darkColors(10)
cols2 <- seqColors(10)
cols3 <- divColors(10)
cols4 <- seqColors2(10)
plot(x, y, col=cols1, xlim=c(1, 13), pch=19, cex=3)
points(x+1, y, col=cols2, pch=19, cex=3)
points(x+2, y, col=cols3, pch=19, cex=3)
points(x+3, y, col=cols4, pch=19, cex=3)
abline(0, 1, lty=2)
abline(-1, 1, lty=2)
abline(-2, 1, lty=2)
abline(-3, 1, lty=2)
}
\keyword{hplot}

