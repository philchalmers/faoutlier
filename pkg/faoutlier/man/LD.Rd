\name{LD}
\alias{LD}
\title{Likelihood Distance}
\usage{
  LD(data, nfact, na.rm = TRUE, digits = 5)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{nfact}{number of factors to extract}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}
}
\description{
  Compute likelihood distances between models when removing
  the \eqn{i_{th}} case.
}
\examples{
\dontrun{
output <- LD(data, 2)
output
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{gCD}}, \code{\link{obs.resid}},
  \code{link{robustMD}}
}
\keyword{cooks}

