#' Robust Mahalanobis 
#' 
#' Compute Mahalanobis distances using a trimmed estimate of the covariance matrix.
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
#' \code{\link{gCD}}, \code{\link{obs.resid}}, \code{link{LD}}
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
	ret$mah <- mahalanobis(data, rob$center, rob$cov)
	ret$mah <- round(ret$mah, digits)
	ret$mah_p <- round(pchisq(ret$mah, ncol(data), lower.tail = FALSE), digits)
	ret
}


