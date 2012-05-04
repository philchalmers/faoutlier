#' Influential case detection methods for factor analysis and SEM
#' 
#' Implements robust Mahalanobis methods, generalized Cook's distances,
#' likelihood ratio tests, model implied residuals, and various 
#' graphical methods to help detect and summarize influential 
#' cases that can affect exploratory and confirmatory factor analyses. The 
#' package can also use \code{OpenMx} for evaluating confirmatory models, 
#' although this version is not yet available through CRAN directly. To obtain an 
#' \code{OpenMx} compatible version do the following:
#' 
#' \describe{ 
#' 	\item{devtools}{If not installed, obtain the \code{devtools} package with 
#' 		\code{install.packages('devtools')}}
#'	\item{install}{Install the github maintained \code{faoutlier} package 
#' 		using \code{devools::install_github('faoutlier', username='philchalmers', branch='omxversion'}}
#'	\item{load}{Now load the new faoutlier normally with \code{library(faoutlier)}}
#' }
#'  
#' 
#'  
#' 
#' @name faoutlier
#' @docType package
#' @title Influential case detection methods for FA and SEM
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @import MASS
# @import OpenMx
#' @import sem
#' @import lattice
#' @keywords package
NULL

#' Description of holzinger data
#' 
#' A sample of 100 simulated cases from the infamous Holzinger dataset
#' using 9 variables. 
#' 
#' 
#' @name holzinger
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL

#' Description of holzinger data with 1 outlier
#' 
#' A sample of 100 simulated cases from the infamous Holzinger dataset
#' using 9 variables, but with 1 outlier added to the dataset. The first row was replaced 
#' by adding 2 to five of the observed variables (odd-numbered items) and subtracting 2 from 
#' the other four observed variables (even-numbered items). 
#' 
#' 
#' @name holzinger.outlier
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL
