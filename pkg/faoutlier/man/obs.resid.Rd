\name{obs.resid}
\alias{obs.resid}
\title{Model predicted residual outliers}
\usage{
  obs.resid(data, nfact, na.rm = TRUE, digits = 5)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{nfact}{number of factors to extract}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}
}
\description{
  Compute model predicted residuals for each variable using
  regression estimated factor scores.
}
\examples{
\dontrun{
output <- obs.resid(data, 2)
output
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{gCD}}, \code{\link{LD}},
  \code{link{robustMD}}
}
\keyword{covariance}

