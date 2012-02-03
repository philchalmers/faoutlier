#' Forward search algorithm for outlier detection
#' 
#' The forward search algorithm begins by selecting a homogeneous subset 
#' of cases based on a maximum likelihood criteria and continues to add individual 
#' cases at each iteration given an acceptance criteria. By default the function
#' add cases that contribute most to the likelihood function and that have 
#' the closest robust mahalanobis distance, however model implied residuals 
#' may be included as well.
#' 
#' 
#' 
#' @aliases forward.search
#' @param data matrix or data.frame 
#' @param model if a single numeric number declares number of factors to extract in 
#' exploratory factor ansysis. If \code{class(model)} is an OpenMx model then a 
#' confirmatory factor analysis is performed instead
#' @param criteria character strings indicating the forward search method
#' Can contain \code{'LD'} for log-liklihood distance, \code{'mah'} for Mahalanobis
#' distance, or \code{'res'} for model implied residuals 
#' @param n.subsets a scalar indicating how many samples to draw to find 
#' a homogenous
#' starting base group
#' @param p.base proportion of sample size to use as the base group
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{LD}}, \code{\link{robustMD}}
#' @keywords forward.search
#' @export forward.search
#' @examples 
#' 
#' \dontrun{
#' output <- forward.search(data, 2)
#' output
#' }
forward.search <- function(data, model, criteria = c('LD', 'mah'), 
	n.subsets = 1000, p.base= .4, na.rm = TRUE, digits = 5)
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
			cat('Remaining iterations:', length(nbaseID), '\n')	
			flush.console()
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
				(model*(model + 1))))
		}
		Cooksstat <- c(NA, Cooksstat)
		orderentered <- c(NA, orderentered)
		LRstat[i+1] <- basemodels[[i+1]]$value * N
		Rhat <- basemodels[[i+1]]$loadings %*% t(basemodels[[i+1]]$loadings)
		diag(Rhat) <- 1
		RMR[i+1] <- sqrt(2*sum(((basemodels[[i+1]]$R - Rhat)^2) / 
			(model*(model + 1))))
		ret <- list(LR=LRstat, RMR=RMR, gCD=Cooksstat, ord=orderentered)		
	}
	class(ret) <- 'forward.search'
	ret
}

#' @S3method print forward.search
print.forward.search <- function(x, stat = 'LR', ...)
{
	if(stat == 'LR') ret <- x$LR
	if(stat == 'RMR') ret <- x$RMR
	if(stat == 'gCD') ret <- x$gCD
	names(ret) <- x$ord
	print(ret)
}

#' @S3method plot forward.search
plot.forward.search <- function(x, y = NULL, stat = 'LR', ...)
{
	id <- x$ord
	if(stat == 'LR') stat2 <- x$LR
	if(stat == 'RMR') stat2 <- x$RMR
	if(stat == 'gCD') stat2 <- x$gCD
	plot(stat2, type = 'h', main = paste(stat, "plot")) 	
}



