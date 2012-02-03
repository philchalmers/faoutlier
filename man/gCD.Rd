\name{gCD}
\alias{gCD}
\title{Generalized Cook's Distance}
\usage{
  gCD(data, model, na.rm = TRUE, digits = 5)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{model}{if a single numeric number declares number
  of factors to extract in exploratory factor ansysis. If
  \code{class(model)} is an OpenMx model then a
  confirmatory factor analysis is performed instead}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}
}
\description{
  Compute generalize Cook's distances (gCD's) for
  exploratory and confirmatory FA. Can return DFBETA matrix
  if requested.
}
\examples{
\dontrun{
output <- gCD(data, 2)
output
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{LD}}, \code{\link{obs.resid}},
  \code{\link{robustMD}}
}
\keyword{cooks}

