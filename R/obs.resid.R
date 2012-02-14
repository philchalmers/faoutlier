#' Model predicted residual outliers 
#' 
#' Compute model predicted residuals for each variable using regression
#' estimated factor scores. 
#' 
#' 
#' @aliases obs.resid
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor ansysis. If \code{class(model)} is an OpenMx model then a 
#' confirmatory factor analysis is performed instead
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{LD}}, \code{\link{robustMD}}
#' @keywords covariance
#' @export obs.resid
#' @examples 
#' 
#' \dontrun{
#' data(holzinger)
#'
#' ###Exploratory
#' nfact <- 2
#' (obs.resid.result <- obs.resid(holzinger, nfact))
#' plot(obs.resid.result)
#'
#' ###Confirmatory
#' manifests <- colnames(holzinger)
#' latents <- c("G")
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
#' (obs.resid.result2 <- obs.resid(holzinger, model))	  
#' plot(obs.resid.result2)
#' }
obs.resid <- function(data, model, na.rm = TRUE, digits = 5)
{	
	ret <- list()
	rownames(data) <- 1:nrow(data)
	if(na.rm) data <- na.omit(data)	
	N <- nrow(data)	
	if(is.numeric(model)){		
		R <- cor(data)
		mod <- factanal(data, model, rotation='none', scores = 'regression')
		scores <- mod$scores
		ret$fascores <- scores
		Lambda <- unclass(mod$loadings)
		Theta <- diag(mod$uniquenesses)	
		e <- data - scores %*% t(Lambda)
		VAR <- diag(mod$uniqueness) %*% solve(R) %*% 
			diag(mod$uniqueness) 
		eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(e)) 
		colnames(eji) <- colnames(e) <- colnames(data)		
		ret$res <- e
		ret$std_res <- eji
	}
	if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){		
		mxMod <- model
		fullmxData <- mxData(cov(data), type="cov",	numObs = N)
		fullMod <- mxRun(mxModel(mxMod, fullmxData), silent = TRUE)
		sigHat <- fullMod@objective@info$expCov
		mat <- fullMod@output$matrices
		nfact <- 1:(ncol(mat[[3]]) - sum(mat[[3]]))		
		n <- ncol(data)	
		L <- matrix(mat[[1]][1:n, n+nfact], ncol = length(nfact))
		Phi <- as.matrix(mat[[2]][nfact+n, nfact+n])
		U <- mat[[2]][1:n, 1:n]
		scores <- t(Phi %*% t(L) %*% solve(sigHat) %*% 
			t(data - colMeans(data)))
		e <- data - scores %*% t(L)
		VAR <- U %*% solve(cov(data)) %*%  U
		eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(e)) 
		colnames(eji) <- colnames(e) <- colnames(data)		
		ret$res <- e
		ret$std_res <- eji	
	}
	ret$id <- rownames(data)
	class(ret) <- 'obs.resid'
	ret
}

#' @S3method print obs.resid
#' @rdname obs.resid
#' @method print obs.resid
#' @param x an object of class \code{obs.resid}
#' @param ... additional parameters to be passed 
print.obs.resid <- function(x, ...)
{
	stat <- c()
	for(i in 1:length(x$id))
		stat[i] <- x$std_res[i, ] %*% x$std_res[i, ]	
	ret <- list(ee = stat)
	print(ret)
	invisible(ret)
}

#' @S3method plot obs.resid
#' @rdname obs.resid
#' @method plot obs.resid
#' @param y a \code{NULL} value ignored by the plotting function
#' @param main the main title of the plot
#' @param ylab the y label of the plot
plot.obs.resid <- function(x, y = NULL, main = 'obs.resid plot', 
	ylab = 'Observed residuals', ...)
{
	ID <- x$id		
	stat <- c()
	for(i in 1:length(x$id))
		stat[i] <- x$std_res[i, ] %*% x$std_res[i, ]
	plot(ID, stat, type = 'h', main = main, ylab = ylab, ...)
}


