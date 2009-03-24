\name{getNgsColorsInfo}
\Rdversion{1.1}
\alias{getNgsColorsInfo}
\title{
  Helper function to extract color information for filenames on
  NimbleGen arrays.
}
\description{
  This function will (try to) extract the color information for
  NimbleGen arrays. This is useful when using \code{read.xysfiles2} to
  parse XYS files for Tiling applications.
}
\usage{
getNgsColorsInfo(path = ".", pattern1 = "_532", pattern2 = "_635", ...)
}

\arguments{
  \item{path}{path where to look for files}
  \item{pattern1}{pattern to match files supposed to go to the
    first channel}
  \item{pattern2}{pattern to match files supposed to go to the second channel}
  \item{\dots}{extra arguments for \code{list.xysfiles}}
}
\details{
  Many NimbleGen samples are identified following the pattern
  sampleID_532.XYS / sampleID_635.XYS.

  The function suggests sample names if all the filenames follow the
  standard above.
}
\value{
  A data.frame with, at least, two columns: 'channel1' and 'channel2'. A
  third column, 'sampleNames', is returned if the filenames follow the
  sampleID_532.XYS / sampleID_635.XYS standard.
}
\author{
  Benilton Carvalho <bcarvalh@jhsph.edu>
}
\keyword{manip}
