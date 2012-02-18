#' Likelihood Distance
#' 
#' Compute likelihood distances between models when removing the \eqn{i_{th}}
#' case.
#' 
#' Note that \code{LD} is not limited to confirmatory factor analysis using
#' OpenMx, and can apply to nearly any modelbeing studied
#' where detection of influential observations is important.
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
#' data(holzinger)
#' data(holzinger.outlier)
#'
#' ###Exploratory
#' nfact <- 3
#' (LDresult <- LD(holzinger, nfact))
#' (LDresult.outlier <- LD(holzinger.outlier, nfact))
#' plot(LDresult)
#' plot(LDresult.outlier)
#'
#' ###Confirmatory
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
#' (LDresult <- LD(holzinger, model))	  
#' (LDresult.outlier <- LD(holzinger.outlier, model))
#' plot(LDresult)
#' plot(LDresult.outlier)
#' }
LD <- function(data, model, na.rm = TRUE, digits = 5)
{		
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
		fullMod <- mxRun(mxModel(mxMod, mxData), silent = TRUE)
		MLmod <- fullMod@output$Minus2LogLikelihood - fullMod@output$SaturatedLikelihood 
		LR <- c()
		for(i in 1:nrow(data)){  
			tmpmxData <- mxData(cov(data[-i, ]), type="cov", 
				numObs = nrow(data)-1)
			tmpMod <- mxRun(mxModel(mxMod, tmpmxData), silent = TRUE)
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
#' @rdname LD
#' @method print LD
#' @param x an object of class \code{LD}
#' @param ncases number of extreme cases to display
#' @param ... additional parameters to be passed 
print.LD <- function(x, ncases = 10, ...)
{
	sorted <- sort(x)
	if(ncases %% 2 != 0) ncases <- ncases + 1
	sorted <- c(sorted[1:(ncases/2)], 
		sorted[(length(sorted)-(ncases/2 + 1)):length(sorted)])
	ret <- matrix(sorted)
	rownames(ret) <- names(sorted)
	colnames(ret) <- 'deltaX2'
	print(ret)
	invisible(ret)
}

#' @S3method plot LD
#' @rdname LD
#' @method plot LD
#' @param y a \code{NULL} value ignored by the plotting function
#' @param type type of plot to use, default displayes points and lines
#' @param main the main title of the plot
#' @param ylab the y label of the plot
plot.LD <- function(x, y = NULL, main = 'Likelihood Distance', 
	type = c('p','h'), ylab = 'LD', ...)
{
	LD <- abs(as.numeric(x))
	ID <- 1:length(x)	
	dat <- data.frame(LD,ID)	
	ret <- xyplot(LD~ID, dat, type = type, main = main, ylab = ylab, ...)
	return(ret)
}

