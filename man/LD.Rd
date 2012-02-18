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
\details{
  Note that \code{LD} is not limited to confirmatory factor
  analysis using OpenMx, and can apply to nearly any
  modelbeing studied where detection of influential
  observations is important.
}
\examples{
\dontrun{
data(holzinger)
data(holzinger.outlier)

###Exploratory
nfact <- 3
(LDresult <- LD(holzinger, nfact))
(LDresult.outlier <- LD(holzinger.outlier, nfact))
plot(LDresult)
plot(LDresult.outlier)

###Confirmatory
manifests <- colnames(holzinger)
latents <- c("F1","F2","F3")
#specify model, mxData not necessary but useful to check if mxRun works
model <- mxModel("Three Factor",
      type="RAM",
      manifestVars = manifests,
      latentVars = latents,
      mxPath(from="F1", to=manifests[1:3]),
	     mxPath(from="F2", to=manifests[4:6]),
	     mxPath(from="F3", to=manifests[7:9]),
      mxPath(from=manifests, arrows=2),
      mxPath(from=latents, arrows=2,
            free=FALSE, values=1.0),
      mxData(cov(holzinger), type="cov", numObs=nrow(holzinger))
	  )

(LDresult <- LD(holzinger, model))
(LDresult.outlier <- LD(holzinger.outlier, model))
plot(LDresult)
plot(LDresult.outlier)
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

