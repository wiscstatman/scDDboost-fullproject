
#' calculate PDD when add random noise in distance matrix
#'
#' @param data normalized preprocessed transcripts
#' @cd index of conditions of cells
#' @param ncores number cores for parallel computing
#' @param K number of subgroups
#' @param D distance matrix of cells
#' @param hp hyper parameters for EBSeq
#' @param Posp parition patterns
#' @param iter max number of iterations for EM
#' @param lambda parameter for random noise
#' @return posterior probabilities under random distance matrix


PDD_random = function(data, cd, K, D, hp, Posp, iter, lambda, seed){
    set.seed(seed)
    E=rexp(ncol(D),rate = lambda)
    n = ncol(D)
    d_mean = sd(D)
    weights=d_mean*E/2
    PD = rep(0, nrow(data))
    for(i in 1:iter){
    R_D = sapply(1:ncol(D), function(i) sapply(1:ncol(D), function(j) {if(i != j){return(D[i,j]+weights[i]+weights[j])}else{return(0)}}))
    ccl = pam(R_D, k = K, diss = T)$clustering
    n1 = table(cd)[1]
    z1<-c(1:K)
    z2<-c(1:K)
    for(i in 1:K){
        ##current index
        cur<-which(ccl==i)
        z1[i]<-length(which(cur<=n1))
        z2[i]<-length(which(cur>n1))
    }
    post = MDD(z1, z2, Posp)
    np = nrow(Posp)
    modified_p = sapply(1:np,function(i) sum(post[which(ref[[K]][,i] == 1)]))
    res = EBS(data,ccl,gcl,sz,iter,hp,Posp)
    DE = res$DEpattern
    PED = DE%*%modified_p
    PDD = 1 - PED
    PD = PD + PDD
    }
    PD = PD/iter
    return(PD)
}
