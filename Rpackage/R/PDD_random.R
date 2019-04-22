
#' calculate PDD when add random noise in distance matrix
#'
#' @param data normalized preprocessed transcripts
#' @param cd condition label
#' @param K number of subgroups
#' @param D distance matrix of cells
#' @param a shape param for weights
#' @param b rate param for weights
#' @param sz size factors 
#' @param hp hyper parameters for EBSeq
#' @param Posp parition patterns
#' @param iter max number of iterations for EM in EBSeq
#' @param REF refinement relation matrix
#' @param random seed
#' @return posterior probabilities under random distance matrix
#' @export

PDD_random = function(data, cd, K, D, a, sz, hp, Posp, iter, REF, stp1, stp2,seed){
    
    
    cstar = genRClus(D,a,K,seed)
    
    
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
    alpha1 = rep(1,K)
    alpha2 = rep(1,K)
    post = MDD(z1, z2, Posp, alpha1, alpha2)
    np = nrow(Posp)
    #modified_p = sapply(1:np,function(i) sum(post[which(ref[[K]][,i] == 1)]))
    modified_p = t(REF) %*% post
    
    if(K >= 2){
        res = EBS(data,cstar,gcl,sz,iter,hp,Posp,stp1,stp2)
        DE = res$DEpattern
    }
    #    else{
    #        message("small number of clusters, using exact EBSeq")
    #       rowmeans = apply(data_counts,1,mean)
    #        tmp = which(rowmeans > 0)
    #        data_tmp = data[tmp,]
    #        res = EBTest(Data = data_tmp, Conditions = cstar, sizeFactors = sz, maxround = 2)
    #        tmpDE = rep(0, nrow(data))
    #        tmpDE[tmp] = res$PPDE
    #        DE = cbind(1 - tmpDE, tmpDE)
    #    }
    PED = DE%*%modified_p
    #PDD = (1 - DE[,1]) * post[1]
    #PDD = PDD / (PED + PDD)
    PDD = 1 - PED
    return(PDD)
}
