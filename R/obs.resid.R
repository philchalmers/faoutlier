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
#' output <- obs.resid(data, 2)
#' output
#' }
obs.resid <- function(data, model, na.rm = TRUE, digits = 5)
{	
	ret <- list()
	rownames(data) <- 1:nrow(data)
	if(na.rm) data <- na.omit(data)	
	N <- nrow(data)	
	if(is.numeric(model)){		
		R <- cor(data)
		mod <- fa(R,model,rotate='none')
		scores <- fa(data,model,rotate='none',scores=TRUE)$scores
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
		fullMod <- mxRun(mxModel(mxMod, fullmxData))
		sigHat <- fullMod@objective@expCov
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
plot.obs.resid <- function(x, y = NULL, main = 'obs.resid plot', 
	ylab = 'Observed residuals', ...)
{
	ID <- x$id		
	stat <- c()
	for(i in 1:length(x$id))
		stat[i] <- x$std_res[i, ] %*% x$std_res[i, ]
	plot(ID, stat, type = 'h', main = main, ylab = ylab, ...)
}


