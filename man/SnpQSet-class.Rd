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
Class \code{\link[Biobase:class.eSet]{eSet}}, directly.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "eSet", distance 2.
Class \code{link[Biobase:class.Versioned]{Versioned}}, by class "eSet", distance 3.
}
\section{Methods}{
No methods defined with class "SnpQSet" in the signature.
}
\references{}
\author{}
\note{}

\seealso{\code{\link[Biobase:class.eSet]{eSet}}, \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, \code{link[Biobase:class.Versioned]{Versioned}}}

\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
