\name{robustMD}
\alias{print.robmah}
\alias{robustMD}
\alias{summary.robmah}
\title{Robust Mahalanobis}
\usage{
  robustMD(data, method = "mve", na.rm = TRUE, digits = 5)

  \method{print}{robmah} (x, ...)

  \method{summary}{robmah} (object, ...)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{method}{type of estimation for robust means and
  covariance (see \code{\link{cov.rob}}}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}

  \item{x}{an object of class \code{robmah}}

  \item{...}{additional parameters to be passed}

  \item{object}{an object of class \code{robmah}}
}
\description{
  Compute Mahalanobis distances using the robust computing
  methods found in the \code{MASS} package.
}
\examples{
\dontrun{
output <- robustMD(data)
output
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{gCD}}, \code{\link{obs.resid}},
  \code{\link{LD}}
}
\keyword{covariance}

