
#' add random noise to distance matrix and generate new EBS distance
#'
#' @param data transcripts
#' @param K number of clusters
#' @param gcl gene cluster
#' @param sz size factor
#' @param D cell distance matrix
#' @param hp hyper parameters
#' @param Posp partition patterns
#' @param iter max number of iterations
#' @lambda parameters for noise exponential distribution
#' @return new distance matrix

EBS_D_random = function(data, K, gcl, sz, D, hp, Posp, iter, lambda){
    E=rexp(ncol(D),rate = lambda)
    n = ncol(D)
    d_mean = mean(D) * n / (n-1)
    weights=d_mean*E/2
    R_D = sapply(1:ncol(D), function(i) sapply(1:ncol(D), function(j) {if(i != j){return(D[i,j]+weights[i]+weights[j])}else{return(0)}}))
    ccl = pam(R_D, k = K, diss = T)$clustering
    res = EBS(data,ccl,gcl,sz,iter,hp,Posp)
    DE = res$DEpattern
    Y = apply(DE, 1, which.max)
    Ebs_D = EBS_D(ccl, Posp, Y)
    return(Ebs_D)
}
