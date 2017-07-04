#' Define a parallel cluster object to be used in internal functions
#'
#' This function defines a object that is placed in a relevant internal environment defined in faoutlier.
#' Internal functions will utilize this object automatically to capitalize on parallel
#' processing architecture. The object defined is a call from \code{parallel::makeCluster()}. Note that
#' if you are defining other parallel objects (for simulation designs, for example) it is not recommended
#' to define a cluster.
#'
#' @aliases setCluster
#' @param spec input that is passed to \code{parallel::makeCluster()}. If no input is given the
#'   maximum number of available local cores will be used
#' @param ... additional arguments to pass to \code{parallel::makeCluster}
#' @param remove logical; remove previously defined cluster object?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#'
#' Chalmers, R. P. & Flora, D. B. (2015). faoutlier: An R Package for Detecting 
#'   Influential Cases in Exploratory and Confirmatory Factor Analysis.
#'   \emph{Applied Psychological Measurement, 39}, 573-574. \doi{10.1177/0146621615597894}
#'
#' Flora, D. B., LaBrish, C. & Chalmers, R. P. (2012). Old and new ideas for data 
#' screening and assumption testing for exploratory and confirmatory factor analysis. 
#'  \emph{Frontiers in Psychology, 3}, 1-21. \doi{10.3389/fpsyg.2012.00055}
#' @keywords parallel
#' @export setCluster
#' @examples
#'
#' \dontrun{
#'
#' #make 4 cores available for parallel computing
#' setCluster(4)
#'
#' #' #stop and remove cores
#' setCluster(remove = TRUE)
#'
#' #use all available cores
#' setCluster()
#'
#' }
setCluster <- function(spec, ..., remove = FALSE){
    if(missing(spec))
        spec <- parallel::detectCores()
    if(remove){
        if(is.null(.faoutlierClusterEnv$CLUSTER)){
            message('There is no visible CLUSTER() definition')
            return(invisible())
        }
        parallel::stopCluster(.faoutlierClusterEnv$CLUSTER)
        .faoutlierClusterEnv$CLUSTER <- NULL
        .faoutlierClusterEnv$ncores <- 1L
        return(invisible())
    }
    if(!is.null(.faoutlierClusterEnv$CLUSTER)){
        message('CLUSTER() has already been defined')
        return(invisible())
    }
    .faoutlierClusterEnv$CLUSTER <- parallel::makeCluster(spec)
    .faoutlierClusterEnv$ncores <- length(.faoutlierClusterEnv$CLUSTER)
    parSapply(.faoutlierClusterEnv$CLUSTER, 1L:.faoutlierClusterEnv$ncores*2L,
              function(x) invisible())
    return(invisible())
}
