\name{LD}
\alias{LD}
\alias{plot.LD}
\alias{print.LD}
\title{Likelihood Distance}
\usage{
  LD(data, model, na.rm = TRUE, digits = 5, ...)

  \method{print}{LD} (x, ncases = 10, ...)

  \method{plot}{LD} (x, y = NULL,
    main = "Likelihood Distance", type = c("p", "h"),
    ylab = "LD", absolute = FALSE, ...)
}
\arguments{
  \item{data}{matrix or data.frame}

  \item{model}{if a single numeric number declares number
  of factors to extract in exploratory factor analysis. If
  \code{class(model)} is a sem (semmod), lavaan
  (character), then a confirmatory approach is performed
  instead}

  \item{na.rm}{logical; remove rows with missing data? Note
  that this is required for EFA analysis}

  \item{digits}{number of digits to round in the final
  result}

  \item{x}{an object of class \code{LD}}

  \item{ncases}{number of extreme cases to display}

  \item{...}{additional parameters to be passed}

  \item{y}{a \code{NULL} value ignored by the plotting
  function}

  \item{type}{type of plot to use, default displays points
  and lines}

  \item{main}{the main title of the plot}

  \item{ylab}{the y label of the plot}

  \item{absolute}{logical; use absolute values instead of
  deviations?}
}
\description{
  Compute likelihood distances between models when removing
  the \eqn{i_{th}} case.
}
\details{
  Note that \code{LD} is not limited to confirmatory factor
  analysis and can apply to nearly any model being studied
  where detection of influential observations is important.
}
\examples{
\dontrun{
data(holzinger)
data(holzinger.outlier)

#Exploratory
nfact <- 3
(LDresult <- LD(holzinger, nfact))
(LDresult.outlier <- LD(holzinger.outlier, nfact))
plot(LDresult)
plot(LDresult.outlier)

#-------------------------------------------------------------------
#Confirmatory with sem
model <- specifyModel()
	  F1 -> Remndrs,    lam11
	  F1 -> SntComp,    lam21
	  F1 -> WrdMean,    lam31
	  F2 -> MissNum,    lam42
	  F2 -> MxdArit,    lam52
	  F2 -> OddWrds,    lam62
	  F3 -> Boots,      lam73
	  F3 -> Gloves,     lam83
	  F3 -> Hatchts,    lam93
	  F1 <-> F1,   NA,     1
	  F2 <-> F2,   NA,     1
	  F3 <-> F3,   NA,     1

(LDresult <- LD(holzinger, model))
(LDresult.outlier <- LD(holzinger.outlier, model))
plot(LDresult)
plot(LDresult.outlier)

#-------------------------------------------------------------------
#Confirmatory with lavaan
model <- 'F1 =~  Remndrs + SntComp + WrdMean
F2 =~ MissNum + MxdArit + OddWrds
F3 =~ Boots + Gloves + Hatchts'

(LDresult <- LD(holzinger, model, orthogonal=TRUE))
(LDresult.outlier <- LD(holzinger.outlier, model, orthogonal=TRUE))
plot(LDresult)
plot(LDresult.outlier)

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
  \code{\link{gCD}}, \code{\link{obs.resid}},
  \code{\link{robustMD}}
}
\keyword{cooks}

