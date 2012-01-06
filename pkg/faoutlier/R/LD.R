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
#' @export LD
#' @examples 
#' 
#' \dontrun{
#' output <- LD(data, 2)
#' output
#' }
LD <- function(data, model, na.rm = TRUE, digits = 5)
{	
	is.installed('OpenMx')	
	rownames(data) <- 1:nrow(data)
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
		mxData <- mxData(cov(data), type="cov",	numObs = nrow(data))
		fullMod <- mxRun(mxModel(mxMod, mxData))
		MLmod <- fullMod@output$Minus2LogLikelihood - fullMod@output$SaturatedLikelihood 
		LR <- c()
		for(i in 1:nrow(data)){  
			tmpmxData <- mxData(cov(data[-i, ]), type="cov", 
				numObs = nrow(data)-1)
			tmpMod <- mxRun(mxModel(mxMod, tmpmxData))
			LR[i] <- tmpMod@output$Minus2LogLikelihood - tmpMod@output$SaturatedLikelihood
		}	
	}
	deltaX2 <- MLmod - LR	
	deltaX2 <- round(deltaX2, digits)
	names(deltaX2) <- rownames(data)
	class(deltaX2) <- 'LD'
	deltaX2
}

#' @S3method print LD
print.LD <- function(x, ncases = 10, ...)
{
	sorted <- sort(x)
	if(ncases %% 2 != 0) ncases <- ncases + 1
	sorted <- c(sorted[1:(ncases/2)], 
		sorted[(length(sorted)-(ncases/2 + 1)):length(sorted)])
	ret <- matrix(sorted)
	rownames(ret) <- names(sorted)
	colnames(ret) <- 'deltaX2'
	ret
}

#' @S3method plot LD
plot.LD <- function(x, y = NULL, main = 'LD plot', ...)
{
	x <- abs(as.numeric(x))
	plot(x, type = 'h', main = main, ylab = 'Absolute LD', 
		xlab = 'Case Number', ...)
}


