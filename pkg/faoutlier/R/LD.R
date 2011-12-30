#' Likelihood Distance
#' 
#' Compute likelihood distances between models when removing the \eqn{i_{th}}
#' case.
#' 
#' 
#' @aliases LD
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor ansysis. If \code{class(model)} is an OpenMx model then a 
#' confirmatory factor analysis is performed instead
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}
#' @keywords cooks
#' @examples 
#' 
#' \dontrun{
#' output <- LD(data, 2)
#' output
#' }
LD <- function(data, model, na.rm = TRUE, digits = 5)
{	
	is.installed('OpenMx')
	require('OpenMx')
	if(na.rm) data <- na.omit(data)
	if(is.numeric(model)){		
		MLmod <- factanal(data,model)$STATISTIC
		LR <- c()
		for(i in 1:nrow(data)){  
			tmp <- factanal(data[-i, ],model)
			LR[i] <- tmp$STATISTIC
		}
	}
	if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){
		mxMod <- model		
		mxData <- OpenMx::mxData(cov(data), type="cov",	numObs = nrow(data))
		fullMod <- OpenMx::mxRun(OpenMx::mxModel(mxMod, mxData))
		MLmod <- fullMod@output$Minus2LogLikelihood - fullMod@output$SaturatedLikelihood 
		LR <- c()
		for(i in 1:nrow(data)){  
			tmpmxData <- OpenMx::mxData(cov(data[-i, ]), type="cov", 
				numObs = nrow(data)-1)
			tmpMod <- OpenMx::mxRun(OpenMx::mxModel(mxMod, tmpmxData))
			LR[i] <- tmpMod@output$Minus2LogLikelihood - tmpMod@output$SaturatedLikelihood
		}	
	}
	deltaX2 <- MLmod - LR	
	deltaX2 <- round(deltaX2, digits)
	class(deltaX2) <- 'LD'
	deltaX2
}

#' @S3method print LD
print.LD <- function(x, ...)
{


}



