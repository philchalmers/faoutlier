\name{gCD}
\alias{gCD}
\alias{plot.gCD}
\alias{print.gCD}
\title{Generalized Cook's Distance}
\usage{
  gCD(data, model, na.rm = TRUE, digits = 5)

  \method{print}{gCD} (x, head = 0.05, DFBETAS = FALSE,
    ...)

  \method{plot}{gCD} (x, y = NULL, main = "gCD plot",
    ylab = "Generalized Cook Distance", ...)
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

  \item{x}{an object of class \code{gCD}}

  \item{head}{a ratio of how many extreme gCD cases to
  display}

  \item{DFBETAS}{logical; attach DFBETA matrix attribute to
  returned result?}

  \item{...}{additional parameters to be passed}

  \item{y}{a \code{NULL} value ignored by the plotting
  function}

  \item{main}{the main title of the plot}

  \item{ylab}{the y label of the plot}
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

