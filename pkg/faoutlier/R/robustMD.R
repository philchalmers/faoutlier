#' Robust Mahalanobis 
#' 
#' Compute Mahalanobis distances using the robust 
#' computing methods found in the \code{MASS} package.
#' 
#' 
#' @aliases robustMD 
#' @param data matrix or data.frame 
#' @param method type of estimation for robust means and covariance
#' (see \code{\link{cov.rob}}
#' @param na.rm logical; remove cases with missing data?
#' @param digits number of digits to round in the final result
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso
#' \code{\link{gCD}}, \code{\link{obs.resid}}, \code{\link{LD}}
#' @keywords covariance
#' @examples 
#' 
#' \dontrun{
#' output <- robustMD(data)
#' output
#' }
robustMD <- function(data, method = 'mve', na.rm = TRUE, digits = 5)
{
	ret <- list()
	if(na.rm) data <- na.omit(data)	
	id <- 1:nrow(data)
	rob <- cov.rob(data, method = method)	
	ret$ID <- id
	ret$mah <- mahalanobis(data, rob$center, rob$cov)
	ret$mah <- round(ret$mah, digits)
	ret$mah_p <- round(pchisq(ret$mah, ncol(data), lower.tail = FALSE), digits)
	ret$normmah <- mahalanobis(data, colMeans(data), cov(data))
	ret$normmah_p <- round(pchisq(ret$normmah, ncol(data), 
		lower.tail = FALSE), digits)
	class(ret) <- 'robmah'
	ret
}

#' @S3method print robmah
print.robmah <- function(x, ...)
{



}


