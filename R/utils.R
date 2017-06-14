#common utility functions

.faoutlierClusterEnv <- new.env(parent=emptyenv())
.faoutlierClusterEnv$ncores <- 1L

myApply <- function(X, MARGIN, FUN, progress = FALSE, ...){
    usepar <- if(MARGIN == 1L){
            ifelse(length(X[,1]) > 2*.faoutlierClusterEnv$ncores, TRUE, FALSE)
        } else {
            ifelse(length(X[1,]) > 2*.faoutlierClusterEnv$ncores, TRUE, FALSE)
        }
    if(!is.null(.faoutlierClusterEnv$CLUSTER) && usepar){
        if(progress){
            return(pbapply::pbapply(cl=.faoutlierClusterEnv$CLUSTER, X=X,
                                     MARGIN=MARGIN, FUN=FUN, ...))
        } else return(t(parallel::parApply(cl=.faoutlierClusterEnv$CLUSTER, X=X,
                                    MARGIN=MARGIN, FUN=FUN, ...)))
    } else {
        if(progress){
            return(pbapply::pbapply(X=X, MARGIN=MARGIN, FUN=FUN, ...))
        } else return(t(apply(X=X, MARGIN=MARGIN, FUN=FUN, ...)))
    }
}

myLapply <- function(X, FUN, progress = FALSE, ...){
    usepar <- ifelse(length(X) > 2*.faoutlierClusterEnv$ncores, TRUE, FALSE)

    if(!is.null(.faoutlierClusterEnv$CLUSTER) && usepar){
        if(progress){
            return(pbapply::pblapply(cl=.faoutlierClusterEnv$CLUSTER, X=X,
                                     fun=FUN, ...))
        } else return(parallel::parLapply(cl=.faoutlierClusterEnv$CLUSTER, X=X,
                                   fun=FUN, ...))
    } else {
        if(progress){
            return(pbapply::pblapply(X=X, FUN=FUN, ...))
        } else return(lapply(X=X, FUN=FUN, ...))
    }
}