\name{oligo-defunct}
\alias{oligo-defunct}
\alias{fitPLM}
\alias{coefs}
\alias{resids}
\title{Defunct Functions in Package 'oligo'}
\description{The functions or variables listed here are no longer part of 'oligo'}
\usage{
fitPLM(...)
coefs(...)
resids(...)
}

\arguments{
	\item{...}{Arguments.}
}

\details{

  	\code{fitPLM} was replaced by \code{fitProbeLevelModel},
  	allowing faster execution and providing more specific
  	models. \code{fitPLM} was based in the code written by Ben
  	Bolstad in the \code{affyPLM} package. However, all the
  	model-fitting functions are now in the package
  	\code{preprocessCore}, on which \code{fitProbeLevelModel}
  	depends.
	
	\code{coefs} and \code{resids}, like \code{fitPLM}, were
	inherited from the \code{affyPLM} package. They were replaced
	respectively by \code{coef} and \code{residuals}, because this
	is how these statistics are called everywhere else in \code{R}.

}
\keyword{manip}
