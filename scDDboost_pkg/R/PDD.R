
#' calculate posterior probabilities of a gene to be differential distributed
#'
#' @param data normalized preprocessed transcripts
#' @param cd conditions label
#' @param ncores number of cores for parallel computing
#' @param K number of subgroups
#' @param D distance matrix of cells or cluster of cells or a given clustering
#' @param hp hyper parameters for EBSeq
#' @param Posp parition patterns
#' @param iter max number of iterations for EM
#' @param lambda parameter for random noise
#' @param random boolean indicator of whether randomzation has been been implemented on distance matrix
#' @param nrandom number of bagging times
#' @return posterior probabilities of a gene to be differential distributed

#' @export


PDD = function(data, cd, ncores, K, D,sz, hp, Posp, iter, random, UP, nrandom){
    #data(ref.RData)
    gcl = 1:nrow(data)

    
    #if(hp == 0){
    #    hp = rep(1, 1 + nrow(data))
    #}
    
    #if(Posp == 0){
    #   Posp = pat(K)[[1]]
    #}
    
    #if(D == 0){
    #    gc = g_cl(data, ncores)
    #    D = cal_D(gc)
    #}
    
    if(!random){
        if(is.matrix(D)){
            ccl = pam(D, k = K, diss = T)$clustering
        }
        else{
            ccl = D
        }
        res = EBS(data,ccl,gcl,sz,iter,hp,Posp)
        DE = res$DEpattern
        n1 = table(cd)[1]
        z1<-rep(0, K)
        z2<-rep(0, K)
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
        PED = DE%*%modified_p
        
        #PDD = (1 - DE[,1]) * post[1]
        #PDD = PDD / (PED + PDD)
        PDD = 1 - PED
        return(PDD)
    }
    else{
        
        fit3 <- suppressMessages(nlminb( start=c(0.1,0.1), objective=LL, x=D, lower=c(0,0) , upper = c(Inf,UP)))
        d0 = fit3$par[1]
        a1 = fit3$par[2]
        a = a1
        
        
        bp <- BiocParallel::MulticoreParam(ncores)
        result = bplapply(1:nrandom, function(i) {PDD_random(data, cd, K, D, a, sz, hp, Posp, iter, i)}, BPPARAM = bp)
        
        
        boot = matrix(0,nrow=length(result[[1]]),ncol = nrandom)
        for(i in 1:nrandom){
            boot[,i]=result[[i]]
        }
         return (rowSums(boot) / nrandom)
    }
    
}
