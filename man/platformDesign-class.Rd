\name{platformDesign-class}
\docType{class}
\alias{platformDesign-class}

\title{Class "platformDesign"}
\description{Container for chip information, this is environment-based
  and is being replaced by \linkS4class{DBPDInfo}.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("platformDesign", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{featureInfo}:}{Object of class \code{"environment"}
      containing vectors describing the array.}
    \item{\code{featureTypeDescription}:}{Object of class \code{"list"}}
    \item{\code{type}:}{Object of class \code{"character"} describing
      the type of the array ("expression", "exon", "SNP", "tiling"}
    \item{\code{nrow}:}{Object of class \code{"numeric"} providing the
      number of rows.}
    \item{\code{ncol}:}{Object of class \code{"numeric"} providing the
      number of columns}
    \item{\code{nwells}:}{Object of class \code{"numeric"} providing the
    number of wells - NimbleGen specific.}
    \item{\code{lookup}:}{Object of class \code{"data.frame"}: a lookup
      table for complex designs by NimbleGen.}
    \item{\code{indexes}:}{Object of class \code{"list"} ~~ }
    \item{\code{platforms}:}{Object of class \code{"character"} ~~ }
    \item{\code{manufacturer}:}{Object of class \code{"character"} ~~ }
    \item{\code{genomebuild}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{PDInfo}"}, directly.
}
\section{Methods}{
No methods defined with class "platformDesign" in the signature.
}
\references{}
\author{}
\note{}

\seealso{
}
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
