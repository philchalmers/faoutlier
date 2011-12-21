\name{forward.search}
\alias{forward.search}
\title{Forward search algorithm for outlier detection}
\usage{
  forward.search(data, nfact, criteria = c("LD", "mah"),
    n.subsets = 1000, p.base = 0.4, na.rm = TRUE,
    digits = 5)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{nfact}{number of factors to extract}

  \item{criteria}{character strings indicating the forward
  search method. Can contain \code{'LD'} for log-liklihood
  distance, \code{'mah'} for Mahalanobis distance, or
  \code{'res'} for model implied residuals.}

  \item{n.subsets}{a scalar indicating how many samples to
  draw to find a homogenous starting base group}

  \item{p.base}{proportion of sample size to use as the
  base group}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}
}
\description{
  Description
}
\examples{
\dontrun{
output <- forward.search(data, 2)
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
\keyword{forward.search}

