\name{getNetAffx}
\alias{getNetAffx}
\alias{getNetAffx-methods}
\alias{getNetAffx,ExpressionSet-method}
\title{NetAffx Biological Annotations}
\description{
  Gets NetAffx Biological Annotations saved in the annotation package
  (Exon and Gene ST Affymetrix arrays).
}
\usage{
getNetAffx(object, type = "probeset")
}

\arguments{
  \item{object}{'ExpressionSet' object (eg., result of rma())}
  \item{type}{Either 'probeset' or 'transcript', depending on what type
  of summaries were obtained.}
}
\details{
  This retrieves NetAffx annotation saved in the (pd) annotation package
  - annotation(object). It is only available for Exon ST and Gene ST
  arrays.

  The 'type' argument should match the summarization target used to
  generate 'object'. The 'rma' method allows for two targets: 'probeset'
  (target='probeset') and 'transcript' (target='core', target='full', target='extended').
}
\value{
  'AnnotatedDataFrame' that can be used as featureData(object)
}
\author{
  Benilton Carvalho
}
\keyword{ classif }

