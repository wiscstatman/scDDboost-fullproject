
#' determine number of clusters
#'
#' @param data transcripts
#' @param D distance matrix
#' @param hp hyper parameter
#' @param K number of clusters
#' @param ncores number of cores for parallel computing
#' @return new distance matrix

det_K = function(data, D, hp, K, ncores, Posp, iter, lambda){
    gcl = 1:nrow(data)
    sz = rep(1, ncol(data))
    bp <- BiocParallel::MulticoreParam(ncores)
    result = bplapply(1:100, function(i) {EBS_D_random(data, K, gcl, sz, D, hp, Posp, iter, lambda)}, BPPARAM = bp)
    
    boot_D = matrix(0,nrow=ncol(D),ncol = ncol(D))
    for(i in 1:100){
        boot_D = boot_D + result[[i]]
    }
    boot_D = boot_D / 100
    
    return(boot_D)
}


