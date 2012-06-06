#' Forward search algorithm for outlier detection
#' 
#' The forward search algorithm begins by selecting a homogeneous subset 
#' of cases based on a maximum likelihood criteria and continues to add individual 
#' cases at each iteration given an acceptance criteria. By default the function
#' add cases that contribute most to the likelihood function and that have 
#' the closest robust Mahalanobis distance, however model implied residuals 
#' may be included as well.
#'
#' Note that \code{forward.search} is not limited to confirmatory factor analysis and 
#' can apply to nearly any model being studied
#' where detection of influential observations is important. 
#' 
#' 
#' @aliases forward.search
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor analysis. If \code{class(model)} is a sem (or OpenMx model if installed 
#' from github) then a confirmatory approach is performed instead
#' @param criteria character strings indicating the forward search method
#' Can contain \code{'LD'} for log-likelihood distance, \code{'mah'} for Mahalanobis
#' distance, or \code{'res'} for model implied residuals 
#' @param n.subsets a scalar indicating how many samples to draw to find 
#' a homogeneous
#' starting base group
#' @param p.base proportion of sample size to use as the base group
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @param print.messages logical; print how many iterations are remaining?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{LD}}, \code{\link{robustMD}}
#' @keywords forward.search
#' @export forward.search
#' @examples 
#' 
#' \dontrun{
#' data(holzinger)
#' data(holzinger.outlier)
#'
#' #Exploratory
#' nfact <- 3
#' (FS <- forward.search(holzinger, nfact))
#' (FS.outlier <- forward.search(holzinger.outlier, nfact))
#' plot(FS)
#' plot(FS.outlier)
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
#' (FS <- forward.search(holzinger, model))	  
#' (FS.outlier <- forward.search(holzinger.outlier, model))
#' plot(FS)
#' plot(FS.outlier)
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
#' (FS <- forward.search(holzinger, model))	  
#' (FS.outlier <- forward.search(holzinger.outlier, model))
#' plot(FS)
#' plot(FS.outlier)
#' }
forward.search <- function(data, model, criteria = c('LD', 'mah'), 
	n.subsets = 1000, p.base= .4, na.rm = TRUE, digits = 5, print.messages = TRUE)
{		
	if(na.rm) data <- na.omit(data)
	N <- nrow(data)
	p <- ncol(data)
	ID <- 1:N
	Samples <- matrix(0, floor(p.base*N), n.subsets)
	for(i in 1:n.subsets)
		Samples[ ,i] <- sample(1:N, floor(p.base*N))	
	if(is.numeric(model)){
		STATISTICS <- rep(NA, n.subsets)
		for(i in 1:n.subsets){
			mod <- factanal(data[Samples[ ,i], ], model)
			STATISTICS[i] <- mod$STATISTIC
		}
		orgbaseID <- baseID <- Samples[ ,(min(STATISTICS) == STATISTICS)]
		nbaseID <- setdiff(ID, baseID)
		basedata <- data[baseID, ]
		basemodels <- list()
		orderentered <- c()
		for (LOOP in 1:length(nbaseID)){			
			tmpcor <- cor(basedata)
			basemodels[[LOOP]] <- mlfact(tmpcor, model)
			basemodels[[LOOP]]$N <- nrow(basedata)
			basemodels[[LOOP]]$R <- tmpcor
 			stat <- c()
			RANK <- rep(0, length(nbaseID))
			if(any(criteria == 'LD')){				
				for(j in 1:length(nbaseID))
					stat[j] <- mlfact(cor(rbind(basedata, data[nbaseID[j], ])), model)$value
				RANK <- RANK + rank(stat)
			}
			if(any(criteria == 'mah')){	
				stat <- mahalanobis(data[nbaseID, ], colMeans(data[baseID, ]), 
					cov(data[baseID, ]))
				RANK <- RANK + rank(stat)	
			}
			if(any(criteria == 'res')){
				stat <- c()
				for(j in 1:length(nbaseID)){					
					tmp <- obs.resid(rbind(basedata, data[nbaseID[j], ]), model)
					stat[j] <- tmp$std_res[nrow(basedata)+1, ] %*% 
						tmp$std_res[nrow(basedata)+1, ]
				}
				RANK <- RANK + rank(stat)
			}
			RANK <- rank(RANK)
			newID <- nbaseID[min(RANK) == RANK]
			if(length(newID) > 1) newID <- newID[1]
			orderentered <- c(orderentered, newID)
			baseID <- c(baseID, newID)
			nbaseID <- setdiff(ID, baseID)
			basedata <- data[baseID, ]
			if(print.messages) {
				cat('Remaining iterations:', length(nbaseID), '\n')	
				flush.console()
			}
		}
		basemodels[[LOOP+1]] <- mlfact(cor(data), model)	
		basemodels[[LOOP+1]]$N <- nrow(data)
		basemodels[[LOOP+1]]$R <- cor(data)
		LRstat <- RMR <- Cooksstat <- c()		
		for(i in 1:(length(basemodels)-1)){
			LRstat[i] <- basemodels[[i]]$value * (length(orgbaseID) + i - 1)
			theta <- basemodels[[i]]$par	
			hess <- basemodels[[i]]$hessian
			theta2 <- basemodels[[i+1]]$par	
			Cooksstat[i] <- (theta-theta2) %*% solve(hess) %*% (theta-theta2)
			Rhat <- basemodels[[i]]$loadings %*% t(basemodels[[i]]$loadings)
			diag(Rhat) <- 1
			RMR[i] <- sqrt(2*sum(((basemodels[[i]]$R - Rhat)^2) /
				(ncol(Rhat)*(ncol(Rhat) + 1))))
		}
		Cooksstat <- c(NA, Cooksstat)
		orderentered <- c(NA, orderentered)
		LRstat[i+1] <- basemodels[[i+1]]$value * N
		Rhat <- basemodels[[i+1]]$loadings %*% t(basemodels[[i+1]]$loadings)
		diag(Rhat) <- 1
		RMR[i+1] <- sqrt(2*sum(((basemodels[[i+1]]$R - Rhat)^2) / 
			(ncol(Rhat)*(ncol(Rhat) + 1))))
		ret <- list(LR=LRstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)		
	}
	if(class(model) == "semmod"){        
	    STATISTICS <- rep(NA, n.subsets)
	    sampleCov <- cov(data)	    
	    for(i in 1:n.subsets){	        
	        samplesemMod <- sem(model, cov(data[-i,]), N-1)
	        STATISTICS[i] <- samplesemMod$criterion * (samplesemMod$N - 1)
	    }        
        orgbaseID <- baseID <- Samples[ ,(min(STATISTICS) == STATISTICS)]
	    nbaseID <- setdiff(ID, baseID)
	    basedata <- data[baseID, ]
	    basemodels <- list()
	    orderentered <- c()
	    for (LOOP in 1:length(nbaseID)){	
	        tmpcov <- cov(basedata)
	        basemodels[[LOOP]] <- sem(model, tmpcov, nrow(basedata))	    	
	        stat <- c()
	        RANK <- rep(0, length(nbaseID))		
			if(any(criteria == 'LD')){	
				for(j in 1:length(nbaseID)){
					tmpcov <- cov(rbind(basedata, data[nbaseID[j], ]))					
					tmpmod <- sem(model, tmpcov, nrow(basedata) + 1)
					stat[j] <- tmpmod$criterion * (tmpmod$N - 1)
				}		
				RANK <- RANK + rank(stat)
			}
			if(any(criteria == 'mah')){	
				stat <- mahalanobis(data[nbaseID, ], colMeans(data[baseID, ]), 
					cov(data[baseID, ]))
				RANK <- RANK + rank(stat)	
			}
			if(any(criteria == 'res')){ 
				stat <- c()
				for(j in 1:length(nbaseID)){					
					tmp <- obs.resid(rbind(basedata, data[nbaseID[j], ]), model)
					stat[j] <- tmp$std_res[nrow(basedata)+1, ] %*% 
						tmp$std_res[nrow(basedata)+1, ]
				}
				RANK <- RANK + rank(stat)
			}
	    	RANK <- rank(RANK)
	    	newID <- nbaseID[min(RANK) == RANK]
	    	if(length(newID) > 1) newID <- newID[1]
	    	orderentered <- c(orderentered, newID)
	    	baseID <- c(baseID, newID)
	    	nbaseID <- setdiff(ID, baseID)
	    	basedata <- data[baseID, ]
	    	if(print.messages){
	    		cat('Remaining iterations:', length(nbaseID), '\n')	
	    		flush.console()		
	    	}
	    }	
	    tmpcov <- cov(data)	    
	    basemodels[[LOOP+1]] <- sem(model, tmpcov, N)
	    LRstat <- RMR <- Cooksstat <- c()		
	    for(i in 1:(length(basemodels)-1)){
	    	LRstat[i] <- basemodels[[i]]$criterion * (basemodels[[i]]$N - 1)			
	    	theta <- basemodels[[i]]$coeff	
	    	vcov <- basemodels[[i]]$vcov
	    	theta2 <- basemodels[[i+1]]$coeff
	    	Cooksstat[i] <- (theta-theta2) %*% vcov %*% (theta-theta2)			
	    	Chat <- basemodels[[i]]$C
	    	C <- basemodels[[i]]$S			
	    	RMR[i] <- sqrt(2*sum(((C - Chat)^2) /
	    		(ncol(C)*(ncol(C) + 1))))
	    }		
	    Cooksstat <- c(NA, Cooksstat)
	    orderentered <- c(NA, orderentered)
	    LRstat[i+1] <- basemodels[[i+1]]$criterion * (basemodels[[i+1]]$N - 1)
	    Chat <- basemodels[[i+1]]$C
	    C <- basemodels[[i+1]]$S
	    RMR[i+1] <- sqrt(2*sum(((C - Chat)^2) /	(ncol(C)*(ncol(C) + 1))))	
	    ret <- list(LR=LRstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)	    
	}
	##OPENMX## if(class(model) == "MxRAMModel" || class(model) == "MxModel" ){	
	##OPENMX## 	STATISTICS <- rep(NA, n.subsets)
	##OPENMX## 	mxMod <- model
	##OPENMX## 	sampleMxData <- mxData(cov(data), type="cov", numObs = nrow(data))
	##OPENMX## 	sampleMxMod <- mxRun(mxModel(mxMod, sampleMxData), silent = TRUE)
	##OPENMX## 	for(i in 1:n.subsets){
	##OPENMX## 		sampleMxData <- mxData(cov(data[Samples[ ,i], ]), type="cov", numObs = nrow(Samples))
	##OPENMX## 		sampleMxMod <- mxRun(mxModel(mxMod, sampleMxData), silent = TRUE)
	##OPENMX## 		STATISTICS[i] <- sampleMxMod@output$Minus2LogLikelihood - 
	##OPENMX## 			sampleMxMod@output$SaturatedLikelihood 			
	##OPENMX## 	}
	##OPENMX## 	orgbaseID <- baseID <- Samples[ ,(min(STATISTICS) == STATISTICS)]
	##OPENMX## 	nbaseID <- setdiff(ID, baseID)
	##OPENMX## 	basedata <- data[baseID, ]
	##OPENMX## 	basemodels <- list()
	##OPENMX## 	orderentered <- c()
	##OPENMX## 	for (LOOP in 1:length(nbaseID)){	
	##OPENMX## 		tmpcov <- cov(basedata)
	##OPENMX## 		Data <- mxData(tmpcov, type="cov", numObs = nrow(Samples))
	##OPENMX## 		basemodels[[LOOP]] <- mxRun(mxModel(mxMod, Data), silent = TRUE)
	##OPENMX## 		stat <- c()
	##OPENMX## 		RANK <- rep(0, length(nbaseID))		
	##OPENMX## 		if(any(criteria == 'LD')){	
	##OPENMX## 			for(j in 1:length(nbaseID)){
	##OPENMX## 				tmpcov <- cov(rbind(basedata, data[nbaseID[j], ]))
	##OPENMX## 				sampleMxData <- mxData(tmpcov, type="cov", numObs = nrow(basedata) + 1)
	##OPENMX## 				sampleMxMod <- mxRun(mxModel(mxMod, sampleMxData), silent = TRUE)
	##OPENMX## 				stat[j] <- sampleMxMod@output$Minus2LogLikelihood - 
	##OPENMX## 					sampleMxMod@output$SaturatedLikelihood	
	##OPENMX## 			}		
	##OPENMX## 			RANK <- RANK + rank(stat)
	##OPENMX## 		}
	##OPENMX## 		if(any(criteria == 'mah')){	
	##OPENMX## 			stat <- mahalanobis(data[nbaseID, ], colMeans(data[baseID, ]), 
	##OPENMX## 				cov(data[baseID, ]))
	##OPENMX## 			RANK <- RANK + rank(stat)	
	##OPENMX## 		}
	##OPENMX## 		if(any(criteria == 'res')){
	##OPENMX## 			stat <- c()
	##OPENMX## 			for(j in 1:length(nbaseID)){					
	##OPENMX## 				tmp <- obs.resid(rbind(basedata, data[nbaseID[j], ]), model)
	##OPENMX## 				stat[j] <- tmp$std_res[nrow(basedata)+1, ] %*% 
	##OPENMX## 					tmp$std_res[nrow(basedata)+1, ]
	##OPENMX## 			}
	##OPENMX## 			RANK <- RANK + rank(stat)
	##OPENMX## 		}
	##OPENMX## 		RANK <- rank(RANK)
	##OPENMX## 		newID <- nbaseID[min(RANK) == RANK]
	##OPENMX## 		if(length(newID) > 1) newID <- newID[1]
	##OPENMX## 		orderentered <- c(orderentered, newID)
	##OPENMX## 		baseID <- c(baseID, newID)
	##OPENMX## 		nbaseID <- setdiff(ID, baseID)
	##OPENMX## 		basedata <- data[baseID, ]
	##OPENMX## 		if(print.messages){
	##OPENMX## 			cat('Remaining iterations:', length(nbaseID), '\n')	
	##OPENMX## 			flush.console()		
	##OPENMX## 		}
	##OPENMX## 	}	
	##OPENMX## 	tmpcov <- cov(data)
	##OPENMX## 	Data <- mxData(tmpcov, type="cov", numObs = nrow(data))
	##OPENMX## 	basemodels[[LOOP+1]] <- mxRun(mxModel(mxMod, Data), silent = TRUE)		
	##OPENMX## 	LRstat <- RMR <- Cooksstat <- c()		
	##OPENMX## 	for(i in 1:(length(basemodels)-1)){
	##OPENMX## 		LRstat[i] <- basemodels[[i]]@output$Minus2LogLikelihood - 
	##OPENMX## 			basemodels[[i]]@output$SaturatedLikelihood			
	##OPENMX## 		theta <- basemodels[[i]]@output$estimate	
	##OPENMX## 		hess <- basemodels[[i]]@output$estimatedHessian		
	##OPENMX## 		theta2 <- basemodels[[i+1]]@output$estimate
	##OPENMX## 		Cooksstat[i] <- (theta-theta2) %*% solve(hess) %*% (theta-theta2)			
	##OPENMX## 		Chat <- basemodels[[i]]@objective@info$expCov			
	##OPENMX## 		C <- basemodels[[i]]@data@observed			
	##OPENMX## 		RMR[i] <- sqrt(2*sum(((C - Chat)^2) /
	##OPENMX## 			(ncol(C)*(ncol(C) + 1))))
	##OPENMX## 	}		
	##OPENMX## 	Cooksstat <- c(NA, Cooksstat)
	##OPENMX## 	orderentered <- c(NA, orderentered)
	##OPENMX## 	LRstat[i+1] <- basemodels[[i+1]]@output$Minus2LogLikelihood - 
	##OPENMX## 			basemodels[[i+1]]@output$SaturatedLikelihood
	##OPENMX## 	Chat <- basemodels[[i+1]]@objective@info$expCov	
	##OPENMX## 	C <- basemodels[[i+1]]@data@observed			
	##OPENMX## 	RMR[i+1] <- sqrt(2*sum(((C - Chat)^2) /
	##OPENMX## 			(ncol(C)*(ncol(C) + 1))))	
	##OPENMX## 	ret <- list(LR=LRstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)
	##OPENMX## }
	class(ret) <- 'forward.search'
	ret
}

#' @S3method print forward.search
#' @rdname forward.search
#' @method print forward.search
#' @param x an object of class \code{forward.search}
#' @param stat type of statistic to use. Could be 'LR', 'RMR', or 'gCD' for 
#' the likelihood ratio, root mean square residual, or generalized Cook's distances,  
#' respectively
#' @param ... additional parameters to be passed
print.forward.search <- function(x, stat = 'LR', ...)
{
	if(stat == 'LR') ret <- x$LR
	if(stat == 'RMR') ret <- x$RMR
	if(stat == 'gCD') ret <- x$gCD
	names(ret) <- x$ord
	return(print(ret))
}

#' @S3method plot forward.search
#' @rdname forward.search
#' @method plot forward.search
#' @param y a \code{null} value ignored by \code{plot}
#' @param main the main title of the plot
#' @param type type of plot to use, default displays points and lines
#' @param ylab the y label of the plot
plot.forward.search <- function(x, y = NULL, stat = 'LR', main = 'Forward Search', 
	type = c('p','h'), ylab = 'obs.resid', ...)
{
	id <- x$ord
	Input <- 1:length(id)
	if(stat == 'LR') stat2 <- x$LR
	if(stat == 'RMR') stat2 <- x$RMR
	if(stat == 'gCD') stat2 <- x$gCD
	dat <- data.frame(stat2,Input,id)	
	ret <- xyplot(stat2~Input, dat, type=type, main=main, ylab=stat, 
		xlab='Input Step', ...)
	return(ret)
}



