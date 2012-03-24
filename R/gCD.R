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
#' exploratory factor ansysis. If \code{class(model)} is a sem (or OpenMx model if installed 
#' from github) then a confirmatory approach is performed instead
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
#' #Exploratory
#' nfact <- 3
#' (gCDresult <- gCD(holzinger, nfact))
#' (gCDresult.outlier <- gCD(holzinger.outlier, nfact))
#' plot(gCDresult)
#' plot(gCDresult.outlier)
#'
#' #Confirmatory with sem
#' model <- specifyModel()
#'	  F1 -> V1,    lam11
#' 	  F1 -> V2,    lam21
#' 	  F1 -> V3,    lam31
#' 	  F2 -> V4,    lam41
#' 	  F2 -> V5,    lam52
#' 	  F2 -> V6,    lam62
#' 	  F3 -> V7,    lam73
#'	  F3 -> V8,    lam83
#' 	  F3 -> V9,    lam93
#' 	  F1 <-> F1,   NA,     1
#' 	  F2 <-> F2,   NA,     1
#' 	  F3 <-> F3,   NA,     1
#' 
#' (gCDresult2 <- gCD(holzinger, model))
#' (gCDresult2.outlier <- gCD(holzinger.outlier, model))
#' plot(gCDresult2)
#' plot(gCDresult2.outlier)
#'
#' #Confirmatory using OpenMx (requires github version, see ?faoutlier)
#' manifests <- colnames(holzinger)
#' latents <- c("F1","F2","F3")
#' #specify model, mxData not necessary but useful to check if mxRun works
#' model <- mxModel("Three Factor",
#'       type="RAM",
#'       manifestVars = manifests,
#'       latentVars = latents,
#'       mxPath(from="F1", to=manifests[1:3]),
#' 	     mxPath(from="F2", to=manifests[4:6]),
#' 	     mxPath(from="F3", to=manifests[7:9]),
#'       mxPath(from=manifests, arrows=2),
#'       mxPath(from=latents, arrows=2,
#'             free=FALSE, values=1.0),
#'       mxData(cov(holzinger), type="cov", numObs=nrow(holzinger))
#' 	  )			
#' 	  
#' (gCDresult2 <- gCD(holzinger, model))	  
#' (gCDresult2.outlier <- gCD(holzinger.outlier, model))
#' plot(gCDresult2)
#' plot(gCDresult2.outlier)
#' }
gCD <- function(data, model, na.rm = TRUE, digits = 5)
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
	    mod <- sem(model, cov(data), N)
	    theta <- mod$coeff		
	    gCD <- c()	
	    DFBETAS <- matrix(0, N, length(theta))
	    for(i in 1:nrow(data)){
	        tmp1 <- cov(data[-i, ])
	        tmp2 <- sem(model, tmp1, N-1)
	        vcovmat <- tmp2$vcov
	        h2 <- tmp2$coeff 			
	        DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
	        gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
	    }	
	    gCD <- round(gCD,digits)
	    DFBETAS <- round(DFBETAS,digits)
	    ret <- list(dfbetas = DFBETAS, gCD = gCD)    
	}
	##OPENMX## if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){		
	##OPENMX## 	mxMod <- model
	##OPENMX## 	fullmxData <- mxData(cov(data), type="cov",	numObs = nrow(data))
	##OPENMX## 	fullMod <- mxRun(mxModel(mxMod, fullmxData), silent = TRUE)
	##OPENMX## 	theta <- fullMod@output$estimate		
	##OPENMX## 	gCD <- c()	
	##OPENMX## 	DFBETAS <- matrix(0,nrow(data),length(theta))
	##OPENMX## 	for(i in 1:nrow(data)){
	##OPENMX## 		tmpmxData <- mxData(cov(data[-i,]), type="cov",	
	##OPENMX## 			numObs = nrow(data)-1)
	##OPENMX## 		tmpMod <- mxRun(mxModel(mxMod, tmpmxData), silent = TRUE)
	##OPENMX## 		vcovmat <- solve(tmpMod@output$estimatedHessian)
	##OPENMX## 		h2 <- tmpMod@output$estimate			
	##OPENMX## 		DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
	##OPENMX## 		gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
	##OPENMX## 	}	
	##OPENMX## 	gCD <- round(gCD,digits)
	##OPENMX## 	DFBETAS <- round(DFBETAS,digits)
	##OPENMX## 	ret <- list(dfbetas = DFBETAS, gCD = gCD)
	##OPENMX## }
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

