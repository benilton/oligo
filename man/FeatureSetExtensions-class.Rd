\name{ExonFeatureSet-class}
\docType{class}
\alias{ExpressionFeatureSet-class}
\alias{ExonFeatureSet-class}
\alias{SnpFeatureSet-class}
\alias{TilingFeatureSet-class}

\title{"FeatureSet" Extensions}
\description{Classes to store data from Expression/Exon/SNP/Tiling arrays at the feature level.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ExonFeatureSet", assayData, manufacturer, platform, exprs, phenoData, featureData, experimentData, annotation, ...)}.
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
Class \code{\linkS4class{FeatureSet}}, directly.
Class \code{\link[Biobase:class.eSet]{eSet}}, by class "FeatureSet", distance 2.
Class \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, by class "FeatureSet", distance 3.
Class \code{link[Biobase:class.Versioned]{Versioned}}, by class "FeatureSet", distance 4.
}
\section{Methods}{
No methods defined with class "ExonFeatureSet" in the signature.
}
\references{}
\author{}
\note{}


\seealso{\code{\link[Biobase:class.eSet]{eSet}}, \code{\link[Biobase:class.VersionedBiobase]{VersionedBiobase}}, \code{link[Biobase:class.Versioned]{Versioned}}}

\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
