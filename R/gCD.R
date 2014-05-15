#' Generalized Cook's Distance
#' 
#' Compute generalize Cook's distances (gCD's) for exploratory 
#' and confirmatory FA. Can return DFBETA matrix if requested.
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
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#'   \code{\link{LD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}, \code{\link{setCluster}}
#' @references
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data screening and assumption testing for 
#' exploratory and confirmatory factor analysis. \emph{Frontiers in Psychology, 3}, 1-21. 
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
#' model <- specifyModel()
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
#' }
gCD <- function(data, model, ...)
{	
    f_numeric <- function(ind, data, model, theta){        
        tmp1 <- cor(data[-ind,])
        tmp2 <- mlfact(tmp1, model)	
        vcovmat <- solve(tmp2$hessian)
        h2 <- tmp2$par 			
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
        ret <- list(dfbeta = DFBETA, gCD = gCD)
    }
    f_sem <- function(ind, data, model, objective, theta, ...){
        tmp2 <- sem::sem(model, data=data[-ind, ], objective=objective, ...)
        vcovmat <- tmp2$vcov
        h2 <- tmp2$coeff 			
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
        ret <- list(dfbeta = DFBETA, gCD = gCD)
    }
    f_lavaan <- function(ind, data, model, theta, ...){
        tmp <- lavaan::sem(model, data[-ind, ], ...)
        vcovmat <- lavaan::vcov(tmp)
        h2 <- lavaan::coef(tmp)
        DFBETA <- (theta - h2)/sqrt(diag(vcovmat))
        gCD <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
        ret <- list(dfbeta = DFBETA, gCD = gCD)
    }
    
	N <- nrow(data)
    index <- as.list(1:N)
	if(is.numeric(model)){	
	    if(any(is.na(data)))
	        stop('Numeric model requires complete dataset (no NA\'s)')
		mod <- mlfact(cor(data), model)
		theta <- mod$par
		tmp <- myLapply(index, FUN=f_numeric, theta=theta, model=model, data=data)
	}	
	if(class(model) == "semmod"){
	    objective <- if(any(is.na(data))) sem::objectiveFIML else sem::objectiveML
	    mod <- sem::sem(model, data=data, objective=objective, ...)
	    theta <- mod$coeff
	    tmp <- myLapply(index, FUN=f_sem, theta=theta, model=model, data=data, 
                        objective=objective, ...)    
	}
	if(class(model) == "character"){
	    mod <- lavaan::sem(model, data=data, ...)
	    theta <- lavaan::coef(mod)
        tmp <- myLapply(index, FUN=f_lavaan, theta=theta, model=model, data=data, ...)    
	}
    gCD <- lapply(tmp, function(x) x$gCD)
    gCD <- do.call(c, gCD)
    dfbetas <- lapply(tmp, function(x) x$dfbeta)
    dfbetas <- do.call(rbind, dfbetas)
    ret <- list(dfbetas = dfbetas, gCD = gCD)    
	class(ret) <- 'gCD'
	ret	
}

#' @rdname gCD
#' @param x an object of class \code{gCD}
#' @param head a ratio of how many extreme gCD cases to display
#' @param DFBETAS logical; attach DFBETA matrix attribute to returned result? 
#' @param ... additional parameters to be passed 
#' @export
print.gCD <- function(x, head = .05, DFBETAS = FALSE, ...)
{
	ID <- 1:length(x$gCD)
	ncases <- floor(length(ID)*(1-head))
	sorted <- cbind(ID, x$gCD)
	ranked <- rank(x$gCD)
	ret <- sorted[ranked >= ncases, ]
	colnames(ret) <- c('ID', 'gCD')
	if(DFBETAS){
		attr(ret,'dfbetas') <- x$dfbetas[ret[ ,1], ]
		rownames(attr(ret,'dfbetas')) <- ret[ ,1]	
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

