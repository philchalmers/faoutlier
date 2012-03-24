#' Model predicted residual outliers 
#' 
#' Compute model predicted residuals for each variable using regression
#' estimated factor scores. 
#' 
#' 
#' @aliases obs.resid
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor ansysis. If \code{class(model)} is a sem (or OpenMx model if installed 
#' from github) then a confirmatory approach is performed instead
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
#' data(holzinger.outlier)
#'
#' #Exploratory
#' nfact <- 3
#' (ORresult <- obs.resid(holzinger, nfact))
#' (ORresult.outlier <- obs.resid(holzinger.outlier, nfact))
#' plot(ORresult)
#' plot(ORresult.outlier)
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
#' (ORresult <- obs.resid(holzinger, model))	  
#' (ORresult.outlier <- obs.resid(holzinger.outlier, model))
#' plot(ORresult)
#' plot(ORresult.outlier)
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
#' (ORresult <- obs.resid(holzinger, model))	  
#' (ORresult.outlier <- obs.resid(holzinger.outlier, model))
#' plot(ORresult)
#' plot(ORresult.outlier)
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
	if(class(model) == "semmod"){
        C <- cov(data)
        vnames <- colnames(C)
        mod <- sem(model, C, N)
        scores <- fscores(mod, data)        
        ret$fascores <- scores
        lnames <- setdiff(colnames(mod$P), vnames)
        Lambda <- mod$A[1:length(vnames), (length(vnames)+1):ncol(mod$A)]
        Theta <- mod$P[1:length(vnames),1:length(vnames)]
        e <- data - scores %*% t(Lambda)
        VAR <- Theta %*% solve(C) %*% Theta
        eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(e))
        colnames(eji) <- colnames(e) <- colnames(data)
        ret$res <- e
        ret$std_res <- eji        
	}
	##OPENMX## if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){		
	##OPENMX## 	mxMod <- model
	##OPENMX## 	fullmxData <- mxData(cov(data), type="cov",	numObs = N)
	##OPENMX## 	fullMod <- mxRun(mxModel(mxMod, fullmxData), silent = TRUE)
	##OPENMX## 	sigHat <- fullMod@objective@info$expCov
	##OPENMX## 	mat <- fullMod@output$matrices
	##OPENMX## 	nfact <- 1:(ncol(mat[[3]]) - sum(mat[[3]]))		
	##OPENMX## 	n <- ncol(data)	
	##OPENMX## 	L <- matrix(mat[[1]][1:n, n+nfact], ncol = length(nfact))
	##OPENMX## 	Phi <- as.matrix(mat[[2]][nfact+n, nfact+n])
	##OPENMX## 	U <- mat[[2]][1:n, 1:n]
	##OPENMX## 	scores <- t(Phi %*% t(L) %*% solve(sigHat) %*% 
	##OPENMX## 		t(data - colMeans(data)))
	##OPENMX## 	e <- data - scores %*% t(L)
	##OPENMX## 	VAR <- U %*% solve(cov(data)) %*%  U
	##OPENMX## 	eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(e)) 
	##OPENMX## 	colnames(eji) <- colnames(e) <- colnames(data)		
	##OPENMX## 	ret$res <- e
	##OPENMX## 	ret$std_res <- eji	
	##OPENMX## }
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
#' @param type type of plot to use, default displayes points and lines
plot.obs.resid <- function(x, y = NULL, main = 'Observed Residuals', 
	type = c('p','h'), ylab = 'obs.resid', ...)
{
	ID <- as.numeric(x$id)
	stat <- c()
	for(i in 1:length(x$id))
		stat[i] <- x$std_res[i, ] %*% x$std_res[i, ]
	dat <- data.frame(stat,ID)	
	ret <- xyplot(stat~ID, dat, type = type, main = main, ylab = ylab, ...)
	return(ret)
}


