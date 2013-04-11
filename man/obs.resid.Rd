\name{obs.resid}
\alias{obs.resid}
\alias{plot.obs.resid}
\alias{print.obs.resid}
\title{Model predicted residual outliers}
\usage{
  obs.resid(data, model, na.rm = TRUE, digits = 5, ...)

  \method{print}{obs.resid} (x, restype = "obs", ...)

  \method{plot}{obs.resid} (x, y = NULL,
    main = "Observed Residuals", type = c("p", "h"),
    restype = "obs", ...)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{model}{if a single numeric number declares number
  of factors to extract in exploratory factor analysis. If
  \code{class(model)} is a sem (semmod), or lavaan
  (character), then a confirmatory approach is performed
  instead}

  \item{na.rm}{logical; remove rows with missing data? Note
  that this is required for EFA analysis and \code{sem}
  fitted models}

  \item{digits}{number of digits to round in the final
  result}

  \item{x}{an object of class \code{obs.resid}}

  \item{restype}{type of residual used, either \code{'obs'}
  for observation value (inner product), \code{'res'} or
  \code{'std_res'} for unstandardized and standardized for
  each variable, respectively}

  \item{...}{additional parameters to be passed}

  \item{y}{a \code{NULL} value ignored by the plotting
  function}

  \item{main}{the main title of the plot}

  \item{type}{type of plot to use, default displays points
  and lines}
}
\description{
  Compute model predicted residuals for each variable using
  regression estimated factor scores.
}
\examples{
\dontrun{
data(holzinger)
data(holzinger.outlier)

#Exploratory
nfact <- 3
(ORresult <- obs.resid(holzinger, nfact))
(ORresult.outlier <- obs.resid(holzinger.outlier, nfact))
plot(ORresult)
plot(ORresult.outlier)

#-------------------------------------------------------------------
#Confirmatory with sem
model <- specifyModel()
   F1 -> Remndrs,    lam11
	  F1 -> SntComp,    lam21
	  F1 -> WrdMean,    lam31
	  F2 -> MissNum,    lam41
	  F2 -> MxdArit,    lam52
	  F2 -> OddWrds,    lam62
	  F3 -> Boots,      lam73
	  F3 -> Gloves,     lam83
	  F3 -> Hatchts,    lam93
	  F1 <-> F1,   NA,     1
	  F2 <-> F2,   NA,     1
	  F3 <-> F3,   NA,     1

(ORresult <- obs.resid(holzinger, model))
(ORresult.outlier <- obs.resid(holzinger.outlier, model))
plot(ORresult)
plot(ORresult.outlier)

#-------------------------------------------------------------------
#Confirmatory with lavaan
model <- 'F1 =~  Remndrs + SntComp + WrdMean
F2 =~ MissNum + MxdArit + OddWrds
F3 =~ Boots + Gloves + Hatchts'

(obs.resid2 <- obs.resid(holzinger, model, orthogonal=TRUE))
(obs.resid2.outlier <- obs.resid(holzinger.outlier, model, orthogonal=TRUE))
plot(obs.resid2)
plot(obs.resid2.outlier)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old
  and new ideas for data screening and assumption testing
  for exploratory and confirmatory factor analysis.
  \emph{Frontiers in Psychology, 3}, 1-21.
}
\seealso{
  \code{\link{gCD}}, \code{\link{LD}},
  \code{\link{robustMD}}
}
\keyword{covariance}

