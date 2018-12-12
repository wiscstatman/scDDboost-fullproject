#' likelihood function for hyperparameters estimation
#'
#' @param D distance matrix of cells
#' @return return hyperparameteres a.
#' @export



LL = function(param, x){
    a0 = 1   #shape for prior
    d0 = param[1]   #rate for prior
    a1 = param[2]  #shape for sampling model, consistent with dividing
    
    n = length(x)
    
    nc = ncol(x)
    
    I = matrix(1,nc,nc)
    
    I = I - diag(nc)
    
    C = d0 + a1 * x
    
    res = (n - nc) * (lgamma(a0 + a1) - lgamma(a0) - lgamma(a1) + a0 * log(d0) + a1 * log(a1)) + sum((a1 - 1) * log(x + diag(nc))) - sum((a0 + a1) * log(C) * I)
    
    return(-res)
}
