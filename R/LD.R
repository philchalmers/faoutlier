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
#' exploratory factor analysis. If \code{class(model)} is a sem (semmod), lavaan (character), 
#' then a confirmatory approach is performed instead
#' @param na.rm logical; remove rows with missing data? Note that this is required for 
#' EFA analysis
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}
#' @references
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data screening and assumption testing for 
#' exploratory and confirmatory factor analysis. \emph{Frontiers in Psychology, 3}, 1-21.
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
#' #-------------------------------------------------------------------
#' #Confirmatory with sem
#' model <- specifyModel()
#'	  F1 -> Remndrs,    lam11
#' 	  F1 -> SntComp,    lam21
#' 	  F1 -> WrdMean,    lam31
#' 	  F2 -> MissNum,    lam42
#' 	  F2 -> MxdArit,    lam52
#' 	  F2 -> OddWrds,    lam62
#' 	  F3 -> Boots,      lam73
#'	  F3 -> Gloves,     lam83
#' 	  F3 -> Hatchts,    lam93
#' 	  F1 <-> F1,   NA,     1
#' 	  F2 <-> F2,   NA,     1
#' 	  F3 <-> F3,   NA,     1
#' 
#' (LDresult <- LD(holzinger, model))	  
#' (LDresult.outlier <- LD(holzinger.outlier, model))
#' plot(LDresult)
#' plot(LDresult.outlier)
#' 
#' #-------------------------------------------------------------------
#' #Confirmatory with lavaan
#' model <- 'F1 =~  Remndrs + SntComp + WrdMean
#' F2 =~ MissNum + MxdArit + OddWrds
#' F3 =~ Boots + Gloves + Hatchts'
#' 
#' (LDresult <- LD(holzinger, model, orthogonal=TRUE))	  
#' (LDresult.outlier <- LD(holzinger.outlier, model, orthogonal=TRUE))
#' plot(LDresult)
#' plot(LDresult.outlier)
#' 
#' }
LD <- function(data, model, na.rm = TRUE, digits = 5, ...)
{		
	rownames(data) <- 1:nrow(data)
	if(na.rm) data <- na.omit(data)
	LR <- c()
	if(is.numeric(model)){		
		MLmod <- factanal(data,model)$STATISTIC		
		for(i in 1:nrow(data)){  
			tmp <- factanal(data[-i, ],model, ...)
			LR[i] <- tmp$STATISTIC
		}
	}
	if(class(model) == "semmod"){
	    MLmod <- sem::sem(model, data=data)
        MLmod <- MLmod$criterion * MLmod$N	    
	    for(i in 1:nrow(data)){  
	        tmp <- sem::sem(model, cov(data[-i, ]), nrow(data) - 1, ...)            
	        LR[i] <- tmp$criterion * (tmp$N - 1)
	    }	    
	}
	if(class(model) == "character"){        
        MLmod <- lavaan::sem(model, data=data, ...)
        MLmod <- MLmod@Fit@test[[1]]$stat
        for(i in 1:nrow(data)){  
            tmp <- lavaan::sem(model, data[-i, ], ...)            
            LR[i] <- tmp@Fit@test[[1]]$stat
        }
	}
	deltaX2 <- LR - MLmod 	
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
	colnames(ret) <- 'LD'
	return(print(ret))	
}

#' @S3method plot LD
#' @rdname LD
#' @method plot LD
#' @param y a \code{NULL} value ignored by the plotting function
#' @param type type of plot to use, default displays points and lines
#' @param main the main title of the plot
#' @param ylab the y label of the plot
#' @param absolute logical; use absolute values instead of deviations?
plot.LD <- function(x, y = NULL, main = 'Likelihood Distance', 
	type = c('p','h'), ylab = 'LD', absolute = FALSE, ...)
{    
	LD <- if(absolute) abs(as.numeric(x)) else as.numeric(x)
	ID <- 1:length(x)	
	dat <- data.frame(LD,ID)	
	ret <- lattice::xyplot(LD~ID, dat, type = type, main = main, ylab = ylab, 
	                       panel = function(...){
                               panel.xyplot(ID, LD, type='h')
                               panel.xyplot(ID, LD, type='p')
                               panel.abline(0)
                               }, ...)
	return(ret)
}

