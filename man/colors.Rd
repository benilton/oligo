\name{darkColors}
\alias{darkColors}
\alias{seqColors}

\title{
Create set of colors, interpolating through a set of preferred colors.
}
\description{
Create set of colors, interpolating through a set of preferred colors.
}
\usage{
darkColors(n)
seqColors(n)
}

\arguments{
  \item{n}{integer determining number of colors to be generated}
}

\details{
darkColors is based on the Dark2 palette in RColorBrewer, therefore
useful to describe qualitative features of the data.

seqColors is based on Blues and generates a gradient of blues, therefore
useful to describe quantitative features of the data.
}

\examples{
x <- 1:10
y <- 1:10
cols1 <- darkColors(10)
cols2 <- seqColors(10)
plot(x, y, col=cols1, xlim=c(1, 11))
points(x+1, y, col=cols2)
}
\keyword{hplot}

