\name{FeatureSet subclasses}
\docType{class}
\alias{ExpressionFeatureSet-class}
\alias{SnpFeatureSet-class}
\alias{TilingFeatureSet-class}

\title{Subclasses of FeatureSet.}
\description{Classes to store Feature Set information for specific types
of array.}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("ExpressionFeatureSet", manufacturer, platform, exprs,
  phenoData, experimentData, annotation)},
\code{new("SnpFeatureSet", manufacturer, platform, exprs,
  phenoData, experimentData, annotation)}, \code{new("TilingFeatureSet", manufacturer, platform, exprs, phenoData, experimentData, annotation)}.
Such objects are extensions of \code{FeatureSet} objects.
}
\section{Slots}{
	 \describe{
    \item{\code{manufacturer}:}{Object of class \code{"character"} ~~ }
    \item{\code{platform}:}{Object of class \code{"character"} ~~ }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
Class \code{"FeatureSet"}, directly.
Class \code{"eSet"}, by class \code{"FeatureSet"}.
}
\section{Methods}{
  Methods defined for \code{FeatureSet}.
}
\seealso{
  \code{\link{eSet-class}}, \code{\link{FeatureSet-class}}.
}
\examples{
##---- Should be DIRECTLY executable !! ----
}
\keyword{classes}
