
#' calculate posterior probabilities of a gene to be differential distributed
#'
#' @param data normalized preprocessed transcripts
#' @cd index of conditions of cells
#' @param ncores number cores for parallel computing
#' @param K number of subgroups
#' @param D distance matrix of cells or cluster of cells
#' @param hp hyper parameters for EBSeq
#' @param Posp parition patterns
#' @param iter max number of iterations for EM
#' @param lambda parameter for random noise
#' @return posterior probabilities of a gene to be differential distributed



PDD = function(data, cd, ncores, K, D, hp, Posp, iter, random, lambda, nrandom){
    #data(ref.RData)
    gcl = 1:nrow(data)
    sz = rep(1, ncol(data))
    
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
        modified_p = sapply(1:np,function(i) sum(post[which(ref[[K]][,i] == 1)]))
        PED = DE%*%modified_p
        PDD = 1 - PED
        return(PDD)
    }
    else{
        bp <- BiocParallel::MulticoreParam(ncores)
        result = bplapply(1:nrandom, function(i) {PDD_random(data, cd, K, D, hp, Posp, iter, lambda, i)}, BPPARAM = bp)
        
        
        boot = matrix(0,nrow=length(result[[1]]),ncol = nrandom)
        for(i in 1:nrandom){
            boot[,i]=result[[i]]
        }
         return (rowSums(boot) / nrandom)
    }
    
}
