\name{SnpCallSet-class}
\docType{class}
\alias{SnpCallSet-class}
\alias{SnpCopyNumberSet-class}
\alias{oligoSnpSet-class}

\title{"SnpSet" Extensions}

\description{Class to store genotype calls / copy number estimates (and
  confidences) obtained via SNP arrays.}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SnpCallSet", assayData, featureData, phenoData, experimentData, annotation, calls, callsConfidence, ...)}.
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
No methods defined with class "SnpCallSet" in the signature.
}
\references{}
\author{}
\note{}

\seealso{\link{Biobase}{eSet}, \link{Biobase}{VersionedBiobase}, \link{Biobase}{Versioned}
}
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
