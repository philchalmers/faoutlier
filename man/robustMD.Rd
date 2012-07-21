\name{robustMD}
\alias{plot.robmah}
\alias{print.robmah}
\alias{robustMD}
\alias{summary.robmah}
\title{Robust Mahalanobis}
\usage{
  robustMD(data, method = "mve", na.rm = TRUE, digits = 5)

  \method{print}{robmah} (x, ...)

  \method{summary}{robmah} (object, gt = 0, ...)

  \method{plot}{robmah} (x, y = NULL, type = "xyplot",
    main, ...)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{method}{type of estimation for robust means and
  covariance (see \code{\link{cov.rob}})}

  \item{na.rm}{logical; remove rows with missing data?}

  \item{digits}{number of digits to round in the final
  result}

  \item{x}{an object of class \code{robmah}}

  \item{...}{additional parameters to be passed}

  \item{object}{an object of class \code{robmah}}

  \item{gt}{only print values with MD's greater than
  \code{gt}}

  \item{y}{empty parameter passed to \code{plot}}

  \item{type}{type of plot to display, can be either
  \code{'qqplot'} or \code{'xyplot'}}

  \item{main}{title for plot. If missing titles will be
  generated automatically}
}
\description{
  Obtain Mahalanobis distances using the robust computing
  methods found in the \code{MASS} package.
}
\examples{
\dontrun{
data(holzinger)
output <- robustMD(holzinger)
output
summary(output)
plot(output)
plot(output, type = 'qqplot')
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old
  and new ideas for data screening and assumption testing
  for exploratory and confirmatory factor analysis.
  \emph{Frontiers in Psychology, 3}, 1-21.
}
\seealso{
  \code{\link{gCD}}, \code{\link{obs.resid}},
  \code{\link{LD}}
}
\keyword{covariance}

