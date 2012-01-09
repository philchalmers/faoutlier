\name{forward.search}
\alias{forward.search}
\title{Forward search algorithm for outlier detection}
\usage{
  forward.search(data, model, criteria = c("LD", "mah"),
    n.subsets = 1000, p.base = 0.4, na.rm = TRUE,
    digits = 5)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{model}{if a single numeric number declares number
  of factors to extract in exploratory factor ansysis. If
  \code{class(model)} is an OpenMx model then a
  confirmatory factor analysis is performed instead}

  \item{criteria}{character strings indicating the forward
  search method Can contain \code{'LD'} for log-liklihood
  distance, \code{'mah'} for Mahalanobis distance, or
  \code{'res'} for model implied residuals}

  \item{n.subsets}{a scalar indicating how many samples to
  draw to find a homogenous starting base group}

  \item{p.base}{proportion of sample size to use as the
  base group}

  \item{na.rm}{logical; remove cases with missing data?}

  \item{digits}{number of digits to round in the final
  result}
}
\description{
  The forward search algorithm begins by selecting a
  homogeneous subset of cases based on a maximum likelihood
  criteria and continues to add individual cases at each
  iteration given an acceptance criteria. By default the
  function add cases that contribute most to the likelihood
  function and that have the closest robust mahalanobis
  distance, however model implied residuals may be included
  as well.
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
  \code{\link{robustMD}}
}
\keyword{forward.search}

