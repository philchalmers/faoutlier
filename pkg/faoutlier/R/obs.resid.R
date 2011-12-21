#' Model predicted residual outliers 
#' 
#' Compute model predicted residuals for each variable using regression
#' estimated factor scores. 
#' 
#' 
#' @aliases obs.resid
#' @param data matrix or data.frame 
#' @param nfact number of factors to extract
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{LD}}, \code{link{robustMD}}
#' @keywords covariance
#' @examples 
#' 
#' \dontrun{
#' output <- obs.resid(data, 2)
#' output
#' }
obs.resid <- function(data, nfact, na.rm = TRUE, digits = 5)
{
	ret <- list()
	if(na.rm) data <- na.omit(data)	
	N <- nrow(data)
	R <- cor(data)
	mod <- fa(R,nfact,rotate='none')
	scores <- fa(data,nfact,rotate='none',scores=TRUE)$scores
	ret$fascores <- scores
	Lambda <- unclass(mod$loadings)
	Theta <- diag(mod$uniquenesses)	
	e <- data - scores %*% t(Lambda)
	
	VAR <- diag(mod$uniqueness) %*% solve(R) %*%  diag(mod$uniqueness) 
	eji <- t(solve(diag(sqrt(diag(VAR)))) %*% t(e)) 
	colnames(eji) <- colnames(e) <- colnames(data)		
	ret$res <- e
	ret$std_res <- eji
	ret
}



