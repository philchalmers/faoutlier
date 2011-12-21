\name{gCD}
\alias{gCD}
\title{Generalized Cook's Distance}
\usage{
  gCD(data, nfact, na.rm = TRUE, digits = 5)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{nfact}{number of factors to extract}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}
}
\description{
  Compute generalize Cook's distances (gCD's) for
  exploratory FA and SEM. Also returns DFBETA matrix.
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
  \code{link{robustMD}}
}
\keyword{cooks}

