
#' MLE for random weighting parameter
#'
#' @param D, distance matrix of cells
#' @param reltol, tolerance of convergence
#' @return MLE of random weighting parameter


rwMLE = function(D,reltol){
    ctrl = list()
    ctrl$rel.tol = reltol
    
    invD = 1 / D
    n = ncol(D)
    for(i in 1:n){
        invD[i,i] = 0
    }
    mm = sum(invD) / (n * (n - 1))
    vv = sum(invD * invD) / (n * (n - 1)) - mm^2
    
    fit3 <- suppressWarnings(nlminb( start=c(2,2), objective=LL, x=D, d0 = mm / vv, lower=c(0,0), control = ctrl))
    
    a0 = fit3$par[1]
    a1 = fit3$par[2]
    a = a0 + a1
    
    return(a)
}
