\docType{package}
\name{faoutlier}
\alias{faoutlier}
\alias{faoutlier-package}
\title{Influential case detection methods for FA and SEM}
\description{
  Influential case detection methods for factor analysis
  and SEM
}
\details{
  Implements robust Mahalanobis methods, generalized Cook's
  distances, likelihood ratio tests, model implied
  residuals, and various graphical methods to help detect
  and summarize influential cases that can affect
  exploratory and confirmatory factor analyses. The package
  can also use \code{OpenMx} for evaluating confirmatory
  models, although this version is not yet available
  through CRAN directly. To obtain an \code{OpenMx}
  compatable version do the following:

  \describe{ \item{devtools}{If not installed, obtain the
  \code{devtools} package with
  \code{install.packages('devtools')}}
  \item{install}{Install the github maintained
  \code{faoutlier} package using
  \code{devools::install_github('faoutlier',
  username='philchalmers', branch='omxversion'}}
  \item{load}{Now load the new faoutlier normally with
  \code{library(faoutlier)}} }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{package}

