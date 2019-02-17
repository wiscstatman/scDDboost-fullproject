
#' modal-cluster for poisson gamma model
#'
#' @param x an integer vector
#' @param y mass parameter
#' @param z shape and scale parameters for gamma distribution
#' @return cluster of x
#' @examples
#' MCP(1:10)

MCP <- function(x,y = 1,z = c(1,1)) {
    if(!is.vector(x)){
        stop("cluster object must be vector")
    }
    if(sum(x %% 1) != 0){
        stop("x must be vector of integers")
    }
    if(length(z) != 2){
        stop("incorrect number of hyper parameters")
    }
    fit<-.Call("MCP",X = x,MASS = y, PARAM=z)
    return(fit)
}


