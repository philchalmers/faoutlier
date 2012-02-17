#' Generalized Cook's Distance
#' 
#' Compute generalize Cook's distances (gCD's) for exploratory 
#' and confirmatory FA. Can return DFBETA matrix if requested.
#' 
#'
#' Note that \code{gCD} is not limited to confirmatory factor analysis using
#' OpenMx, and can apply to nearly any model being studied
#' where detection of influential observations is important. 
#'
#' 
#' @aliases gCD
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor ansysis. If \code{class(model)} is an OpenMx model then a 
#' confirmatory factor analysis is performed instead
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{LD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}
#' @keywords cooks
#' @export gCD
#' @examples 
#' 
#' \dontrun{
#' data(holzinger)
#' data(holzinger.outlier)
#'
#' ###Exploratory
#' nfact <- 3
#' (gCDresult <- gCD(holzinger, nfact))
#' (gCDresult.outlier <- gCD(holzinger.outlier, nfact))
#' plot(gCDresult)
#' plot(gCDresult.outlier)
#'
#' ###Confirmatory
#' manifests <- colnames(holzinger)
#' latents <- c("G")
#' #specify model, mxData not necessary but useful to check if mxRun works
#' model <- mxModel("One Factor",
#'      type="RAM",
#'       manifestVars = manifests,
#'       latentVars = latents,
#'       mxPath(from=latents, to=manifests),
#'       mxPath(from=manifests, arrows=2),
#'      mxPath(from=latents, arrows=2,
#'             free=FALSE, values=1.0),
#'       mxData(cov(holzinger), type="cov", numObs=nrow(holzinger))
#'	  )
#'	  
#' (gCDresult2 <- gCD(holzinger, model))	  
#' (gCDresult2.outlier <- gCD(holzinger.outlier, model))
#' plot(gCDresult2)
#' plot(gCDresult2.outlier)
#' }
gCD <- function(data, model, na.rm = TRUE, digits = 5)
{	
	if(na.rm) data <- na.omit(data)
	if(is.numeric(model)){	
		mod <- mlfact(cor(data), model)
		theta <- mod$par 
		vcovmat <- mod$hessian
		gCD <- c()	
		DFBETAS <- matrix(0,nrow(data),length(theta))
		for(i in 1:nrow(data)){
			tmp1 <- cor(data[-i,])
			tmp2 <- mlfact(tmp1, model)	  
			h2 <- tmp2$par 			
			DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
			gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
		}	
		gCD <- round(gCD,digits)
		DFBETAS <- round(DFBETAS,digits)
		ret <- list(dfbetas = DFBETAS, gCD = gCD)
	}		
	if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){		
		mxMod <- model
		fullmxData <- mxData(cov(data), type="cov",	numObs = nrow(data))
		fullMod <- mxRun(mxModel(mxMod, fullmxData), silent = TRUE)
		theta <- fullMod@output$estimate
		vcovmat <- solve(fullMod@output$estimatedHessian)
		gCD <- c()	
		DFBETAS <- matrix(0,nrow(data),length(theta))
		for(i in 1:nrow(data)){
			tmpmxData <- mxData(cov(data[-i,]), type="cov",	
				numObs = nrow(data)-1)
			tmpMod <- mxRun(mxModel(mxMod, tmpmxData), silent = TRUE)
			h2 <- tmpMod@output$estimate			
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
	print(ret)
	invisible(ret)
}

#' @S3method plot gCD
#' @rdname gCD
#' @method plot gCD
#' @param y a \code{NULL} value ignored by the plotting function
#' @param main the main title of the plot
#' @param type type of plot to use, default displayes points and lines
#' @param ylab the y label of the plot
plot.gCD <- function(x, y = NULL, main = 'Generalized Cook Distance', 
	type = c('p','h'), ylab = 'gCD', ...)
{
	ID <- 1:length(x$gCD)		
	ret <- xyplot(gCD~ID, x, type = type, main = main, ylab = ylab, ...)
	return(ret)
}

