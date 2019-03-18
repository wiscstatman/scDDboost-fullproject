
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
#' @param reltol relative tolerance for optim on weighting paramters
#' @return posterior probabilities of a gene to be differential distributed

#' @export


PDD = function(data, cd, ncores,D, random = T, epi = 1, Upper = 1000, nrandom = 30, iter = 20,reltol = 1e-3){
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
    sz = tryCatch({MedianNorm(data)},error = function(e){
        rep(1, ncol(D))
    })
    alpha = 0.4
    beta = 2
    
    hp = rep(beta, 1 + nrow(data))
    hp[1] = alpha
    
    
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
        if(K >= 2){
            res = EBS(data,ccl,gcl,sz,iter,hp,Posp)
            DE = res$DEpattern
        }
        else if(K == 1){
            message("There is only one cluster, no postive")
            return(rep(0,nrow(data)))
        }
        #        else
        #        {
        #            message("small number of clusters, using exact EBSeq")
        #            rowmeans = apply(data_counts,1,mean)
        #           tmp = which(rowmeans > 0)
        #            data_tmp = data[tmp,]
        #            res = EBTest(Data = data_tmp, Conditions = ccl, sizeFactors = sz, maxround = 3)
        #            tmpDE = rep(0, nrow(data))
        #            tmpDE[tmp] = res$PPDE
        #            DE = cbind(1 - tmpDE, tmpDE)
        #        }
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
        
        if(K == 1){
            message("There is only one cluster, no postive")
            return(rep(0,nrow(data)))
        }
        
        Posp = pat(K)[[1]]
        REF = g_ref(Posp)
        
        ctrl = list()
        ctrl$rel.tol = reltol
        
        fit3 <- suppressWarnings(nlminb( start=c(2,2), objective=LL, x=D, d0 = 1, lower=c(0,0), control = ctrl))
        a0 = fit3$par[1]
        a1 = fit3$par[2]
        a = a0 + a1
        #b = a1
        #message(paste0("param of weights: ", a1))
        
        bp <- BiocParallel::MulticoreParam(ncores)
        result = bplapply(1:nrandom, function(i) {PDD_random(data, cd, K, D, a, sz, hp, Posp, iter, REF, i)}, BPPARAM = bp)
        
        
        boot = matrix(0,nrow=length(result[[1]]),ncol = nrandom)
        for(i in 1:nrandom){
            boot[,i]=result[[i]]
        }
         return (rowSums(boot) / nrandom)
    }
    
}
