
#' calculate PDD when add random noise in distance matrix
#'
#' @param data normalized preprocessed transcripts
#' @param cd condition label
#' @param K number of subgroups
#' @param D distance matrix of cells
#' @param hp hyper parameters for EBSeq
#' @param Posp parition patterns
#' @param iter max number of iterations for EM in EBSeq
#' @param lambda parameter for random noise
#' @return posterior probabilities under random distance matrix
#' @export

PDD_random = function(data, cd, K, D, sz, hp, Posp, iter, lambda, seed){
    set.seed(seed)
    ## E=rexp(ncol(D),rate = lambda)
    n = ncol(D)
    d_mean = mean(D)
    p_ = 1 / K
    noise = matrix(rbinom(n^2, prob=1-p_,size=1),nrow=n)
    ## to make noise be a symmetric matrix
    up_ = 1 * upper.tri(noise, diag = FALSE)
    noise = noise * up_
    noise = noise + t(noise)
    
    ## weights = d_mean * E
    PD = rep(0, nrow(data))
    #R_D = D_c + weights %o% rep(1,n) + rep(1,n) %o% weights
    R_D = R_D + noise * 0.0005
//    ccl = pam(R_D, k = K, diss = T)$clustering
    R_D = as.dist(R_D)
    hc = hclust(R_D)
    ccl = cutree(hc,k = K)
    gcl = 1:nrow(data)
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
    #modified_p = sapply(1:np,function(i) sum(post[which(ref[[K]][,i] == 1)]))
    modified_p = t(ref[[K]]) %*% post
    res = EBS(data,ccl,gcl,sz,iter,hp,Posp)
    DE = res$DEpattern
    PED = DE%*%modified_p
    #PDD = (1 - DE[,1]) * post[1]
    #PDD = PDD / (PED + PDD)
    PDD = 1 - PED
    return(PDD)
}
