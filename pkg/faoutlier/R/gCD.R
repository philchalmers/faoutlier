#' Generalized Cook's Distance
#' 
#' Compute generalize Cook's distances (gCD's) for exploratory FA 
#' and SEM. Also returns DFBETA matrix.
#' 
#' 
#' @aliases gCD
#' @param data matrix or data.frame 
#' @param nfact number of factors to extract
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{LD}}, \code{\link{obs.resid}}, \code{link{robustMD}}
#' @keywords cooks
#' @examples 
#' 
#' \dontrun{
#' output <- gCD(data, 2)
#' output
#' }
gCD <- function(data, nfact, na.rm = TRUE, digits = 5)
{
	if(na.rm) data <- na.omit(data)
	if(is.numeric(nfact)){		
		theta <- mlfact(cor(data), nfact)$par 
		gCD <- c()	
		DFBETAS <- matrix(0,nrow(data),length(theta))
		for(i in 1:nrow(data)){
			tmp1 <- cor(data[-i,])
			tmp2 <- mlfact(tmp1, nfact)	  
			h2 <- tmp2$par 
			vcovmat <- solve(tmp2$hessian)
			DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
			gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
		}	
		gCD <- round(gCD,digits)
		DFBETAS <- round(DFBETAS,digits)
		return(list(dfbetas = DFBETAS, gCD = gCD))
	}		
	if(class(nfact) == "MxRAMModel" || class(nfact) == "MxModel" ){		
		mxMod <- nfact
		fullmxData <- OpenMx::mxData(cov(data), type="cov",	numObs = nrow(data))
		fullMod <- OpenMx::mxRun(OpenMx::mxModel(mxMod, fullmxData))
		theta <- fullMod@output$estimate
		gCD <- c()	
		DFBETAS <- matrix(0,nrow(data),length(theta))
		for(i in 1:nrow(data)){
			tmpmxData <- OpenMx::mxData(cov(data[-i,]), type="cov",	numObs = nrow(data)-1)
			tmpMod <- OpenMx::mxRun(OpenMx::mxModel(mxMod, tmpmxData))
			h2 <- tmpMod@output$estimate
			vcovmat <- solve(tmpMod@output$estimatedHessian)
			DFBETAS[i, ] <- (theta - h2)/sqrt(diag(vcovmat))
			gCD[i] <- t(theta - h2) %*%  vcovmat %*% (theta - h2)  
		}	
		gCD <- round(gCD,digits)
		DFBETAS <- round(DFBETAS,digits)
		return(list(dfbetas = DFBETAS, gCD = gCD))
	}
}



