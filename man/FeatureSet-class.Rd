\name{FeatureSet-class}
\docType{class}
\alias{FeatureSet-class}

\title{Classes for Feature level data}
\description{Virtual class to store feature level data.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("FeatureSet", assayData, manufacturer, platform, exprs, phenoData, featureData, experimentData, annotation, ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"}  }
    \item{\code{platform}:}{Object of class \code{"character"}  }
    \item{\code{assayData}:}{Object of class \code{"AssayData"}  }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}  }
    \item{\code{annotation}:}{Object of class \code{"character"}  }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}  }
  }
}
\section{Extends}{
Class \code{\linkS4class{eSet}}, directly.
Class \code{\linkS4class{VersionedBiobase}}, by class "eSet", distance 2.
Class \code{\linkS4class{Versioned}}, by class "eSet", distance 3.
}
\section{Methods}{
No methods defined with class "FeatureSet" in the signature.
}
\references{}
\author{}
\note{}

\seealso{\linkS4class{eSet}, \linkS4class{VersionedBiobase}, \linkS4class{Versioned}}

\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
