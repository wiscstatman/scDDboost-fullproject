
#' calculate posterior probabilities of a gene to be differential distributed
#'
#' @param data normalized preprocessed transcripts
#' @param cd conditions label
#' @param ncores number of cores for parallel computing
#' @param D distance matrix of cells or cluster of cells or a given clustering
#' @param epi tol for change of validity score in determining number of clusters
#' @param random boolean indicator of whether randomzation has been been implemented on distance matrix
#' @param Upper bound for hyper parameters optimization
#' @param nrandom number of random generated distance matrix
#' @param iter max number of iterations for EM
#' @return posterior probabilities of a gene to be differential distributed

#' @export


PDD = function(data, cd, ncores,D, random = T, epi = 0.15, Upper = 1000, nrandom = 30, iter = 20){
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
    sz = rep(1, ncol(D))
    
    hp = rep(1, 1 + nrow(data))
    
    
    if(!random){
        if(is.matrix(D)){
            K = detK(D,epi)
            message(paste0("estimated number of subtypes: ",K))
            ccl = pam(D, k = K, diss = T)$clustering
        }
        else{
            if(is.numeric(D))
            {
                ccl = D
                K = max(ccl)
            }
            else
            {
               stop("input ccl must either be a dissimilarity measure or a partition for cells")
            }
            
        }
        
        
        Posp = pat(K)[[1]]
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
        alpha1 = rep(1,K)
        alpha2 = rep(1,K)
        post = MDD(z1, z2, Posp, alpha1, alpha2)
        np = nrow(Posp)
        #modified_p = sapply(1:np,function(i) sum(post[which(ref[[K]][,i] == 1)]))
        REF = g_ref(Posp)
        modified_p = t(REF) %*% post
        PED = DE%*%modified_p
        
        #PDD = (1 - DE[,1]) * post[1]
        #PDD = PDD / (PED + PDD)
        PDD = 1 - PED
        return(PDD)
    }
    else{
        K = detK(D,epi)
        message(paste0("estimated number of subtypes: ",K))
        Posp = pat(K)[[1]]
        REF = g_ref(Posp)
        
        fit3 <- suppressWarnings(nlminb( start=c(0.1,2), objective=LL, x=D, lower=c(0,2) , upper = c(Upper,Upper)))
        a0 = 1
        d0 = fit3$par[1]
        a1 = fit3$par[2]
        a = a1 + a0
        b = a1
        #message(paste0("param of weights: ", a1))
        
        bp <- BiocParallel::MulticoreParam(ncores)
        result = bplapply(1:nrandom, function(i) {PDD_random(data, cd, K, D, a,b, sz, hp, Posp, iter, REF, i)}, BPPARAM = bp)
        
        
        boot = matrix(0,nrow=length(result[[1]]),ncol = nrandom)
        for(i in 1:nrandom){
            boot[,i]=result[[i]]
        }
         return (rowSums(boot) / nrandom)
    }
    
}
