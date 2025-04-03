#' Generalized Cook's Distance
#'
#' Compute generalize Cook's distances (gCD's) for exploratory
#' and confirmatory FA. Can return DFBETA matrix if requested.
#' If mirt is used, then the values will be associated with the unique response patterns instead.
#'
#' Note that \code{gCD} is not limited to confirmatory factor analysis and
#' can apply to nearly any model being studied
#' where detection of influential observations is important.
#'
#' @aliases gCD
#' @param data matrix or data.frame
#' @param model if a single numeric number declares number of factors to extract in
#'   exploratory factor analysis (requires complete dataset, i.e., no missing).
#'   If \code{class(model)} is a sem (semmod), or lavaan (character),
#'   then a confirmatory approach is performed instead
#' @param vcov_drop logical; should the variance-covariance matrix of the parameter
#'   estimates be based on the unique \code{data[-i, ]} models
#'   (Pek and MacCallum, 2011) or original \code{data}?
#' @param progress logical; display the progress of the computations in the console?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#'   \code{\link{LD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}, \code{\link{setCluster}}
#' @references
#'
#' Chalmers, R. P. & Flora, D. B. (2015). faoutlier: An R Package for Detecting
#'   Influential Cases in Exploratory and Confirmatory Factor Analysis.
#'   \emph{Applied Psychological Measurement, 39}, 573-574. \doi{10.1177/0146621615597894}
#'
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data
#' screening and assumption testing for exploratory and confirmatory factor analysis.
#'  \emph{Frontiers in Psychology, 3}, 1-21. \doi{10.3389/fpsyg.2012.00055}
#'
#' Pek, J. & MacCallum, R. C. (2011). Sensitivity Analysis in Structural
#'   Equation Models: Cases and Their Influence. Multivariate Behavioral Research,
#'   46(2), 202-228.
#'
#' @keywords cooks
#' @export gCD
#' @examples
#'
#' \dontrun{
#'
#' #run all gCD functions using multiple cores
#' setCluster()
#'
#' #Exploratory
#' nfact <- 3
#' (gCDresult <- gCD(holzinger, nfact))
#' (gCDresult.outlier <- gCD(holzinger.outlier, nfact))
#' plot(gCDresult)
#' plot(gCDresult.outlier)
#'
#' #-------------------------------------------------------------------
#' #Confirmatory with sem
#' model <- sem::specifyModel()
#'    F1 -> Remndrs,    lam11
#' 	  F1 -> SntComp,    lam21
#' 	  F1 -> WrdMean,    lam31
#' 	  F2 -> MissNum,    lam41
#' 	  F2 -> MxdArit,    lam52
#' 	  F2 -> OddWrds,    lam62
#' 	  F3 -> Boots,      lam73
#'	  F3 -> Gloves,     lam83
#' 	  F3 -> Hatchts,    lam93
#' 	  F1 <-> F1,   NA,     1
#' 	  F2 <-> F2,   NA,     1
#' 	  F3 <-> F3,   NA,     1
#'
#' (gCDresult2 <- gCD(holzinger, model))
#' (gCDresult2.outlier <- gCD(holzinger.outlier, model))
#' plot(gCDresult2)
#' plot(gCDresult2.outlier)
#'
#' #-------------------------------------------------------------------
#' #Confirmatory with lavaan
#' model <- 'F1 =~  Remndrs + SntComp + WrdMean
#' F2 =~ MissNum + MxdArit + OddWrds
#' F3 =~ Boots + Gloves + Hatchts'
#'
#' (gCDresult2 <- gCD(holzinger, model, orthogonal=TRUE))
#' (gCDresult2.outlier <- gCD(holzinger.outlier, model, orthogonal=TRUE))
#' plot(gCDresult2)
#' plot(gCDresult2.outlier)
#'
#' # categorical data with mirt
#' library(mirt)
#' data(LSAT7)
#' dat <- expand.table(LSAT7)
#' model <- mirt.model('F = 1-5')
#' result <- gCD(dat, model)
#' plot(result)
#'
#' mod <- mirt(dat, model)
#' res <- mirt::residuals(mod, type = 'exp')
#' cbind(res, gCD=round(result$gCD, 3))
#'
#' }
gCD <- function(data, model, vcov_drop = FALSE, progress = TRUE, ...)
{
    f_numeric <- function(ind, data, model, theta, inv_vcov, vcov_drop){
        tmp1 <- cor(data[-ind,])
        tmp2 <- mlfact(tmp1, model)
        h2 <- tmp2$par
        inv_vcovmat <- if(vcov_drop) tmp2$hessian else inv_vcov
        vcovmat <- solve(inv_vcovmat)
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*%  inv_vcovmat %*% (theta - h2)
        list(dfbeta = DFBETA, gCD = gCD)
    }
    f_sem <- function(ind, data, model, objective, theta, vcov, vcov_drop, ...){
        tmp2 <- sem::sem(model, data=data[-ind, ], objective=objective, ...)
        h2 <- tmp2$coeff
        vcovmat <- if(vcov_drop) tmp2$vcov else vcov
        inv_vcovmat <- solve(vcovmat)
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*%  inv_vcovmat %*% (theta - h2)
        list(dfbeta = DFBETA, gCD = gCD)
    }
    f_lavaan <- function(ind, data, model, theta, vcov, vcov_drop, ...){
        tmp <- lavaan::sem(model, data[-ind, ], ...)
        h2 <- lavaan::coef(tmp)
        vcovmat <- if(vcov_drop) lavaan::vcov(tmp) else vcov
        inv_vcovmat <- solve(vcovmat)
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*% inv_vcovmat %*% (theta - h2)
        list(dfbeta = DFBETA, gCD = gCD)
    }
    f_mirt <- function(ind, data, large, model, theta, sv, vcov, vcov_drop, ...){
        large$Freq[[1L]][ind] <- large$Freq[[1L]][ind] - 1L
        tmp <- mirt::mirt(data, model, large=large, SE=vcov_drop,
                          pars=sv, verbose=FALSE, ...)
        h2 <- mirt::extract.mirt(tmp, 'parvec')
        vcovmat <- if(vcov_drop) mirt::vcov(tmp) else vcov
        inv_vcovmat <- solve(vcovmat)
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*%  inv_vcovmat %*% (theta - h2)
        list(dfbeta = DFBETA, gCD = gCD)
    }

	N <- nrow(data)
    index <- as.list(1:N)
    inv_vcov <- NULL
	if(is.numeric(model)){
	    if(any(is.na(data)))
	        stop('Numeric model requires complete dataset (no NA\'s)')
		mod <- mlfact(cor(data), model)
		theta <- mod$par
		if(!vcov_drop) inv_vcov <- mod$hessian
		tmp <- myLapply(index, FUN=f_numeric, progress=progress,
		                theta=theta, model=model, data=data,
		                inv_vcov=inv_vcov, vcov_drop=vcov_drop)
	} else if(is(model, "semmod")){
	    objective <- if(any(is.na(data))) sem::objectiveFIML else sem::objectiveML
	    mod <- sem::sem(model, data=data, objective=objective, ...)
	    theta <- mod$coeff
	    if(!vcov_drop) vcov <- mod$vcov
	    tmp <- myLapply(index, FUN=f_sem, progress=progress,
	                    theta=theta, model=model, data=data,
                        objective=objective, vcov=vcov,
	                    vcov_drop=vcov_drop, ...)
	} else if(is.character(model)){
	    mod <- lavaan::sem(model, data=data, ...)
	    theta <- lavaan::coef(mod)
	    if(!vcov_drop) vcov <- lavaan::vcov(mod)
        tmp <- myLapply(index, FUN=f_lavaan, progress=progress,
                        theta=theta, model=model, data=data, vcov=vcov,
                        vcov_drop=vcov_drop, ...)
	} else if(is(model, "mirt.model")){
	    large <- mirt::mirt(data=data, model=model, large = 'return')
	    index <- matrix(1L:length(large$Freq[[1L]]))
	    mod <- mirt::mirt(data=data, model=model, large=large, verbose=FALSE,
	                      SE=!vcov_drop, ...)
	    theta <- mirt::extract.mirt(mod, 'parvec')
	    sv <- mirt::mod2values(mod)
	    if(!vcov_drop) vcov <- mirt::vcov(mod)
	    tmp <- myLapply(index, FUN=f_mirt, progress=progress,
	                    theta=theta, model=model, data=data,
	                    large=large, sv=sv,
	                    vcov=vcov, vcov_drop=vcov_drop, ...)
	} else stop('model class not supported')

    gCD <- lapply(tmp, function(x) x$gCD)
    gCD <- do.call(c, gCD)
    dfbetas <- lapply(tmp, function(x) x$dfbeta)
    dfbetas <- do.call(rbind, dfbetas)
    if(!is(model, "mirt.model"))
        names(gCD) <- rownames(data)
    ret <- list(dfbetas = dfbetas, gCD = gCD)
	class(ret) <- 'gCD'
	ret
}

#' @rdname gCD
#' @param x an object of class \code{gCD}
#' @param ncases number of extreme cases to display
#' @param DFBETAS logical; return DFBETA matrix in addition to gCD? If TRUE, a list is returned
#' @param ... additional parameters to be passed
#' @export
print.gCD <- function(x, ncases = 10, DFBETAS = FALSE, ...)
{
    sorted <- x$gCD[order(abs(x$gCD), decreasing = TRUE)]
    ret <- matrix(sorted[1:ncases])
    rownames(ret) <- names(sorted[1:ncases])
    colnames(ret) <- 'gCD'
	if(DFBETAS){
	    dfbetas <- x$dfbetas[names(x$gCD) %in% rownames(ret), ]
	    rownames(dfbetas) <- rownames(ret)
	    ret <- list(gCD=ret, dfbetas=dfbetas)
	}
	return(print(ret))
}

#' @rdname gCD
#' @param y a \code{NULL} value ignored by the plotting function
#' @param main the main title of the plot
#' @param type type of plot to use, default displays points and lines
#' @param ylab the y label of the plot
#' @export
plot.gCD <- function(x, y = NULL, main = 'Generalized Cook Distance',
	type = c('p','h'), ylab = 'gCD', ...)
{
	ID <- 1:length(x$gCD)
	ret <- lattice::xyplot(gCD~ID, x, type = type, main = main, ylab = ylab, ...)
	return(ret)
}

