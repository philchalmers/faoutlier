#common utility functions

faoutlierClusterEnv <- new.env()
faoutlierClusterEnv$ncores <- 1L

myApply <- function(X, MARGIN, FUN, ...){
    if(!is.null(faoutlierClusterEnv$CLUSTER)){
        return(t(parallel::parApply(cl=faoutlierClusterEnv$faoutlierCluster, X=X, 
                                    MARGIN=MARGIN, FUN=FUN, ...)))
    } else {
        return(t(apply(X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    }
}

myLapply <- function(X, FUN, ...){
    if(!is.null(faoutlierClusterEnv$CLUSTER)){
        return(parallel::parLapply(cl=faoutlierClusterEnv$faoutlierCluster, X=X, 
                                   fun=FUN, ...))
    } else {
        return(lapply(X=X, FUN=FUN, ...))
    }
}