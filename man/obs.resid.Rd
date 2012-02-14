\name{obs.resid}
\alias{obs.resid}
\alias{plot.obs.resid}
\alias{print.obs.resid}
\title{Model predicted residual outliers}
\usage{
  obs.resid(data, model, na.rm = TRUE, digits = 5)

  \method{print}{obs.resid} (x, ...)

  \method{plot}{obs.resid} (x, y = NULL,
    main = "obs.resid plot", ylab = "Observed residuals",
    ...)
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

  \item{x}{an object of class \code{obs.resid}}

  \item{...}{additional parameters to be passed}

  \item{y}{a \code{NULL} value ignored by the plotting
  function}

  \item{main}{the main title of the plot}

  \item{ylab}{the y label of the plot}
}
\description{
  Compute model predicted residuals for each variable using
  regression estimated factor scores.
}
\examples{
\dontrun{
data(holzinger)

###Exploratory
nfact <- 2
(obs.resid.result <- obs.resid(holzinger, nfact))
plot(obs.resid.result)

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

(obs.resid.result2 <- obs.resid(holzinger, model))
plot(obs.resid.result2)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{gCD}}, \code{\link{LD}},
  \code{\link{robustMD}}
}
\keyword{covariance}

