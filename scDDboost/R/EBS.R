
#' accelerated empirical bayesian
#'
#' @param data single cell expression matrix, row as genes column as cells
#' @param conditions partition of cells
#' @param gclus partition of genes
#' @param sf size factors
#' @param iter maximum iteration step of EM
#' @param hyper hyper parameters for beta distributions
#' @param different pattern of partitions
#' @return posterior probability of mean expression pattern

EBS <- function(data,conditions,gclus,sf,iter = 10,hyper,PP) {
    if(!is.matrix(data)){
        stop("data must be a numerical matrix")
    }
    if(length(conditions) != ncol(data)){
        stop("incorrect length of conditions")
    }
    if(length(gclus) != nrow(data)){
        stop("incorrect length of gene cluster")
    }
    if(length(sf) != ncol(data)){
        stop("incorrect length of size factors")
    }
    if(length(hyper) != nrow(data) + 1){
        stop("incorrect length of hyper parameters")
    }
    if(!is.matrix(PP)){
        stop("partition pattern must be a numerical matrix")
    }
    
    .Call('EBS', X=data,Y=conditions,Z=gclus,W=sf, iter=iter,hyper=hyper,part=PP )
}



