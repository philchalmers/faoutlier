\name{forward.search}
\alias{forward.search}
\alias{plot.forward.search}
\alias{print.forward.search}
\title{Forward search algorithm for outlier detection}
\usage{
  forward.search(data, model, criteria = c("LD", "mah"),
    n.subsets = 1000, p.base = 0.4, na.rm = TRUE,
    digits = 5)

  \method{print}{forward.search} (x, stat = "LR", ...)

  \method{plot}{forward.search} (x, y = NULL, stat = "LR",
    main = "Forward Search", type = c("p", "h"),
    ylab = "obs.resid", ...)
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

  \item{x}{an object of class \code{forward.search}}

  \item{stat}{type of statistic to use. Could be 'LR',
  'RMR', or 'gCD' for the likelihood ratio, root mean
  square residual, or generalized Cook's disntance,
  respectively}

  \item{...}{additional parameters to be passed}

  \item{y}{a \code{null} value ignored by \code{plot}}

  \item{main}{the main title of the plot}

  \item{type}{type of plot to use, default displayes points
  and lines}

  \item{ylab}{the y label of the plot}
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
data(holzinger)

###Exploratory
nfact <- 2
(forward.search.result <- forward.search(holzinger, nfact))
plot(forward.search.result)

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

(forward.search.result2 <- forward.search(holzinger, model))
plot(forward.search.result2)
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

