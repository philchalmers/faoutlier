\name{LD}
\alias{LD}
\alias{plot.LD}
\alias{print.LD}
\title{Likelihood Distance}
\usage{
  LD(data, model, na.rm = TRUE, digits = 5)

  \method{print}{LD} (x, ncases = 10, ...)

  \method{plot}{LD} (x, y = NULL,
    main = "Likelihood Distance", type = c("p", "h"),
    ylab = "LD", ...)
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

  \item{x}{an object of class \code{LD}}

  \item{ncases}{number of extreme cases to display}

  \item{...}{additional parameters to be passed}

  \item{y}{a \code{NULL} value ignored by the plotting
  function}

  \item{type}{type of plot to use, default displayes points
  and lines}

  \item{main}{the main title of the plot}

  \item{ylab}{the y label of the plot}
}
\description{
  Compute likelihood distances between models when removing
  the \eqn{i_{th}} case.
}
\examples{
\dontrun{
data(holzinger)

###Exploratory
nfact <- 2
(LDresult <- LD(holzinger, nfact))
plot(LDresult)

###Confirmatory
manifests <- colnames(holzinger)
latents <- c("G")
model <- mxModel("One Factor",
     type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from=latents, to=manifests),
      mxPath(from=manifests, arrows=2),
     mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxData(cov(holzinger), type="cov", numObs=nrow(holzinger))
	  )

(LDresult2 <- LD(holzinger, model))
plot(LDresult2)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{gCD}}, \code{\link{obs.resid}},
  \code{\link{robustMD}}
}
\keyword{cooks}

