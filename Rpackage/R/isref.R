
#' check refinement relation between two clusters
#'
#' @param x a cluster
#' @param y a cluster
#' @return whether x refines y




isref <- function(x,y) {
    if(length(x) != length(y)){
        stop("x, y must have equal lenth")
    }
    if(sum(x %% 1) != 0 || sum(y %% 1) != 0){
        stop("x, y must be integer vectors")
    }
    .Call("isref", X=x,Y=y)
}
