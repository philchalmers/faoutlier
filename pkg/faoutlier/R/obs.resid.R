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
#' @examples 
#' 
#' \dontrun{
#' output <- obs.resid(data, 2)
#' output
#' }
obs.resid <- function(data, model, na.rm = TRUE, digits = 5)
{
	is.installed('OpenMx')	
	ret <- list()
	if(na.rm) data <- na.omit(data)	
	N <- nrow(data)
	R <- cor(data)
	mod <- fa(R,model,rotate='none')
	scores <- fa(data,model,rotate='none',scores=TRUE)$scores
	ret$fascores <- scores
	Lambda <- unclass(mod$loadings)
	Theta <- diag(mod$uniquenesses)	
	e <- data - scores %*% t(Lambda)
	
	VAR <- diag(mod$uniqueness) %*% solve(R) %*%  diag(mod$uniqueness) 
	eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(e)) 
	colnames(eji) <- colnames(e) <- colnames(data)		
	ret$res <- e
	ret$std_res <- eji
	class(ret) <- 'obs.resid'
	ret
}

#' @S3method print obs.resid
print.obs.resid <- function(x, ...)
{


}



