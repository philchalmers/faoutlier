#' Likelihood Distance
#' 
#' Compute likelihood distances between models when removing the \eqn{i_{th}}
#' case.
#' 
#' Note that \code{LD} is not limited to confirmatory factor analysis and 
#' can apply to nearly any model being studied
#' where detection of influential observations is important.
#'
#' @aliases LD
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor ansysis. If \code{class(model)} is a sem (or OpenMx model if installed 
#' from github) then a confirmatory approach is performed instead
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
#' #Exploratory
#' nfact <- 3
#' (LDresult <- LD(holzinger, nfact))
#' (LDresult.outlier <- LD(holzinger.outlier, nfact))
#' plot(LDresult)
#' plot(LDresult.outlier)
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
#' (LDresult <- LD(holzinger, model))	  
#' (LDresult.outlier <- LD(holzinger.outlier, model))
#' plot(LDresult)
#' plot(LDresult.outlier)

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
	if(class(model) == "semmod"){
	    MLmod <- sem(model, cov(data), nrow(data))
	    LR <- c()
	    for(i in 1:nrow(data)){  
	        tmp <- sem(model, cov(data[-i, ]), nrow(data) - 1)            
	        LR[i] <- tmp$criterion * (tmp$N - 1)
	    }	    
	}
	##OPENMX## if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){
	##OPENMX## 	mxMod <- model		
	##OPENMX## 	mxData <- mxData(cov(data), type="cov",	numObs = nrow(data))
	##OPENMX## 	fullMod <- mxRun(mxModel(mxMod, mxData), silent = TRUE)
	##OPENMX## 	MLmod <- fullMod@output$Minus2LogLikelihood - fullMod@output$SaturatedLikelihood 
	##OPENMX## 	LR <- c()
	##OPENMX## 	for(i in 1:nrow(data)){  
	##OPENMX## 		tmpmxData <- mxData(cov(data[-i, ]), type="cov", 
	##OPENMX## 			numObs = nrow(data)-1)
	##OPENMX## 		tmpMod <- mxRun(mxModel(mxMod, tmpmxData), silent = TRUE)
	##OPENMX## 		LR[i] <- tmpMod@output$Minus2LogLikelihood - tmpMod@output$SaturatedLikelihood
	##OPENMX## 	}	
	##OPENMX## }
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

