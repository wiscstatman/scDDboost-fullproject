
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
    n = ncol(D)
    fit3 <- suppressMessages(nlminb( start=c(0.1,0.1), objective=LL, x=D, lower=c(0,0) , upper = c(2,Inf)))
    a0 = fit3$par[1]
    d0 = fit3$par[2]
    a = a0 + 1
    e <- rgamma(n,shape=(a + 1), rate=(a) )
    bar = D/outer(e,e,"+")
    dst.star <- as.dist(bar)
    cstar = cutree(hclust(dst.star), k = K)
    
    gcl = 1:nrow(data)
    n1 = table(cd)[1]
    z1<-c(1:K)
    z2<-c(1:K)
    for(i in 1:K){
        ##current index
        cur<-which(cstar==i)
        z1[i]<-length(which(cur<=n1))
        z2[i]<-length(which(cur>n1))
    }
    post = MDD(z1, z2, Posp)
    np = nrow(Posp)
    #modified_p = sapply(1:np,function(i) sum(post[which(ref[[K]][,i] == 1)]))
    modified_p = t(ref[[K]]) %*% post
    res = EBS(data,cstar,gcl,sz,iter,hp,Posp)
    DE = res$DEpattern
    PED = DE%*%modified_p
    #PDD = (1 - DE[,1]) * post[1]
    #PDD = PDD / (PED + PDD)
    PDD = 1 - PED
    return(PDD)
}
