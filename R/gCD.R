#' Generalized Cook's Distance
#' 
#' Compute generalize Cook's distances (gCD's) for exploratory 
#' and confirmatory FA. Can return DFBETA matrix if requested.
#' 
#'
#' Note that \code{gCD} is not limited to confirmatory factor analysis and 
#' can apply to nearly any model being studied
#' where detection of influential observations is important. 
#'
#' 
#' @aliases gCD
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor analysis. If \code{class(model)} is a sem (semmod), or lavaan (character), 
#' then a confirmatory approach is performed instead
#' @param na.rm logical; remove rows with missing data? Note that this is required for 
#' EFA analysis
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{LD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}
#' @references
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data screening and assumption testing for 
#' exploratory and confirmatory factor analysis. \emph{Frontiers in Psychology, 3}, 1-21. 
#' @keywords cooks
#' @export gCD
#' @examples 
#' 
#' \dontrun{
#' data(holzinger)
#' data(holzinger.outlier)
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
gCD <- function(data, model, na.rm = TRUE, digits = 5, ...)
{	
	if(na.rm) data <- na.omit(data)
	N <- nrow(data)
	if(is.numeric(model)){	
		mod <- mlfact(cor(data), model)
		theta <- mod$par 		
		gCD <- c()	
		DFBETAS <- matrix(0, N, length(theta))
		for(i in 1:N){
			tmp1 <- cor(data[-i,])
			tmp2 <- mlfact(tmp1, model)	
			vcovmat <- solve(tmp2$hessian)
			h2 <- tmp2$par 			
			DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
			gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
		}	
		gCD <- round(gCD,digits)
		DFBETAS <- round(DFBETAS,digits)
		ret <- list(dfbetas = DFBETAS, gCD = gCD)
	}	
	if(class(model) == "semmod"){
	    mod <- sem::sem(model, data=data, ...)
	    theta <- mod$coeff		
	    gCD <- c()	
	    DFBETAS <- matrix(0, N, length(theta))
	    for(i in 1:nrow(data)){
	        tmp1 <- cov(data[-i, ])
	        tmp2 <- sem::sem(model, tmp1, N-1, ...)
	        vcovmat <- tmp2$vcov
	        h2 <- tmp2$coeff 			
	        DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
	        gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
	    }	
	    gCD <- round(gCD,digits)
	    DFBETAS <- round(DFBETAS,digits)
	    ret <- list(dfbetas = DFBETAS, gCD = gCD)    
	}
	if(class(model) == "character"){      
        if(!require(lavaan)) require(lavaan)
	    mod <- lavaan::sem(model, data=data, ...)
	    theta <- coef(mod)
	    gCD <- c()    
	    DFBETAS <- matrix(0, N, length(theta))
	    for(i in 1:nrow(data)){
	        tmp <- lavaan::sem(model, data[-i, ], ...)
	        vcovmat <- vcov(tmp)
	        h2 <- coef(tmp)
	        DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
	        gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
	    }
	    gCD <- round(gCD,digits)
	    DFBETAS <- round(DFBETAS,digits)
	    ret <- list(dfbetas = DFBETAS, gCD = gCD)    
	}
	class(ret) <- 'gCD'
	ret	
}

#' @S3method print gCD
#' @rdname gCD
#' @method print gCD
#' @param x an object of class \code{gCD}
#' @param head a ratio of how many extreme gCD cases to display
#' @param DFBETAS logical; attach DFBETA matrix attribute to returned result? 
#' @param ... additional parameters to be passed 
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

#' @S3method plot gCD
#' @rdname gCD
#' @method plot gCD
#' @param y a \code{NULL} value ignored by the plotting function
#' @param main the main title of the plot
#' @param type type of plot to use, default displays points and lines
#' @param ylab the y label of the plot
plot.gCD <- function(x, y = NULL, main = 'Generalized Cook Distance', 
	type = c('p','h'), ylab = 'gCD', ...)
{
	ID <- 1:length(x$gCD)		
	ret <- lattice::xyplot(gCD~ID, x, type = type, main = main, ylab = ylab, ...)
	return(ret)
}

