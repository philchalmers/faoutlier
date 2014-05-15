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
#'   exploratory factor analysis (requires complete dataset, i.e., no missing). 
#'   If \code{class(model)} is a sem (semmod), or lavaan (character), 
#'   then a confirmatory approach is performed instead
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#'   \code{\link{gCD}}, \code{\link{obs.resid}}, \code{\link{robustMD}}, \code{\link{setCluster}}
#' @references
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data screening and assumption testing for 
#' exploratory and confirmatory factor analysis. \emph{Frontiers in Psychology, 3}, 1-21.
#' @keywords cooks
#' @export LD
#' @examples 
#' 
#' \dontrun{
#'
#' #run all LD functions using multiple cores
#' setCluster()
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
LD <- function(data, model, ...)
{	
    f_numeric <- function(ind, data, model, ...){
        tmp <- factanal(data[-ind, ], model, ...)
        tmp$STATISTIC
    }
    f_sem <- function(ind, data, model, objective, ...){
        tmp <- sem::sem(model, data=data[-ind, ], objective=objective, ...)
        tmp$criterion * (tmp$N - 1)
    }
    f_lavaan <- function(ind, data, model, ...){
        tmp <- lavaan::sem(model, data=data[-ind, ], ...)            
        tmp@Fit@test[[1L]]$stat
    }
    
	rownames(data) <- 1:nrow(data)
    index <- matrix(1L:nrow(data))
	if(is.numeric(model)){		
	    if(any(is.na(data)))
	        stop('Numeric model requires complete dataset (no NA\'s)')
		MLmod <- factanal(data, model, ...)$STATISTIC
		LR <- myApply(index, MARGIN=1L, FUN=f_numeric, data=data, model=model, ...)
	}
	if(class(model) == "semmod"){
        objective <- if(any(is.na(data))) sem::objectiveFIML else sem::objectiveML
	    MLmod <- sem::sem(model, data=data, objective=objective, ...)
        MLmod <- MLmod$criterion * (MLmod$N - 1)
	    LR <- myApply(index, MARGIN=1L, FUN=f_sem, data=data, model=model, objective=objective, ...)
	}
	if(class(model) == "character"){
        MLmod <- lavaan::sem(model, data=data, ...)
        MLmod <- MLmod@Fit@test[[1]]$stat
        LR <- myApply(index, MARGIN=1L, FUN=f_lavaan, data=data, model=model, ...)
	}
	deltaX2 <- LR - MLmod 	
	deltaX2 <- deltaX2
	names(deltaX2) <- rownames(data)
	class(deltaX2) <- 'LD'
	deltaX2
}

#' @rdname LD
#' @param x an object of class \code{LD}
#' @param ncases number of extreme cases to display
#' @param digits number of digits to round in the printed result
#' @param ... additional parameters to be passed 
#' @export
print.LD <- function(x, ncases = 10, digits = 5, ...)
{
	sorted <- sort(x)
	if(ncases %% 2 != 0) ncases <- ncases + 1
	sorted <- c(sorted[1:(ncases/2)], 
		sorted[(length(sorted)-(ncases/2 + 1)):length(sorted)])
	ret <- matrix(sorted)
	rownames(ret) <- names(sorted)
	colnames(ret) <- 'LD'
	return(print(round(ret, digits)))	
}

#' @rdname LD
#' @param y a \code{NULL} value ignored by the plotting function
#' @param type type of plot to use, default displays points and lines
#' @param main the main title of the plot
#' @param ylab the y label of the plot
#' @param absolute logical; use absolute values instead of deviations?
#' @export
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

