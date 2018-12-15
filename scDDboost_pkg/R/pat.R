
#' generating partition patterns
#'
#' @param K number of elements
#' @return all possible partition of K elements



pat<-function(K){
    if(K %% 1 != 0){
        stop("input must be integer")
    }
    .Call('pat', K=K)
}
