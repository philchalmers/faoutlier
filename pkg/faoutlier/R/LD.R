#' Likelihood Distance
#' 
#' Compute likelihood distances between models when removing the \eqn{i_{th}}
#' case.
#' 
#' 
#' @aliases LD
#' @param data matrix or data.frame 
#' @param nfact number of factors to extract
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{obs.resid}}, \code{link{robustMD}}
#' @keywords cooks
#' @examples 
#' 
#' \dontrun{
#' output <- LD(data, 2)
#' output
#' }
LD <- function(data, nfact, na.rm = TRUE, digits = 5)
{	
	if(na.rm) data <- na.omit(data)
	if(is.numeric(nfact)){		
		MLmod <- factanal(data,nfact)$STATISTIC
		LR <- c()
		for(i in 1:nrow(data)){  
			tmp <- factanal(data[-i, ],nfact)
			LR[i] <- tmp$STATISTIC
		}
	}
	if(class(nfact) == "MxRAMModel" || class(nfact) == "MxModel" ){
		mxMod <- nfact		
		mxData <- OpenMx::mxData(cov(data), type="cov",	numObs = nrow(data))
		fullMod <- OpenMx::mxRun(OpenMx::mxModel(mxMod, mxData))
		MLmod <- fullMod@output$Minus2LogLikelihood - fullMod@output$SaturatedLikelihood 
		LR <- c()
		for(i in 1:nrow(data)){  
			tmpmxData <- OpenMx::mxData(cov(data[-i, ]), type="cov", numObs = nrow(data)-1)
			tmpMod <- OpenMx::mxRun(OpenMx::mxModel(mxMod, tmpmxData))
			LR[i] <- tmpMod@output$Minus2LogLikelihood - tmpMod@output$SaturatedLikelihood
		}	
	}
	deltaX2 <- MLmod - LR	
	deltaX2 <- round(deltaX2, digits)
	deltaX2
}



