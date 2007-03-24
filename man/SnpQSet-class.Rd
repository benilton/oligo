\name{SnpQSet-class}
\docType{class}
\alias{SnpQSet-class}

\title{Class "SnpQSet"}
\description{Class to store Snp (Quantifications) Summaries}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("SnpQSet",
  assayData, senseThetaA, senseThetaB, antisenseThetaA, antisenseThetaB,
  phenoData, featureData, experimentData, annotation)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
No methods defined with class "SnpQSet" in the signature.
}
\references{}
\author{}
\note{}

\seealso{\linkS4class[Biobase]{eSet}, \linkS4class[Biobase]{VersionedBiobase}, \linkS4class[Biobase]{Versioned}}

\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
