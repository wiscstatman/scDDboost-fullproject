library(SingleCellExperiment)
library(scDDboost)
library(ggplot2)

load("Fucci.RData")

K = 6

#ccl = pam(D_c,K,diss = T)$clustering
sz = rep(1,ncol(data_counts))
iter = 10
hp = rep(2, 1 + nrow(data_counts))
hp[1] = 0.4
Posp = pat(K)[[1]]
stp1 = 1e-6
stp2 = 1e-2
gcl = 1:nrow(data_counts)

alphaList = list()

betaList = list()

rlist = list()

clist = list()



for(i in 1:10){
    cstar = genRClus(D_c, a, K, i)
    clist[[i]] = cstar
    res = EBS(data_counts, cstar, gcl, sz, iter, hp, Posp, stp1,
                stp2)
    
    alphaList[[i]] = res$Alpha
    betaList[[i]] = res$Beta
    rlist[[i]] = res$r[,1]
}


fitVal = function(x,clus,r,K,alpha,Beta){
    mn = tapply(x,clus,sum)
    cn = tapply(rep(1,length(x)),clus,sum)
    #p = mn / (r  + mn)
    
    post_alpha = cn * r + alpha
    post_beta = mn + Beta
    
    
    fit = c()
    for(i in 1:K){
        for(j in 1:cn[i]){
            local_q  = rbeta(1,post_alpha[i],post_beta[i])
            local_x = rnbinom(1,r,prob = local_q)
            fit = c(fit,local_x)
        }
    }
    return(fit)
}


########### plot ###########


## random sample 6 non zero genes

set.seed(1600)

rM = rowMeans(data_counts)

#NZ = sample(which(rM > 0),6)

callout = lsz(newPDD,0.05)

names(callout) = genenames[callout]

NZ = sample(callout,6)

## random sample 3 non-constant shape genes
names(trouble) = genenames[trouble]
NC = sample(trouble,3)


total = c(NZ,NC)

den = list()
den_o = list()
for(i in 1:9){
    x = data_counts[total[i],]
    den_o[[i]] = ecdf(x)
    local_den = list()
    for(j in 1:10){
        tmp = fitVal(x,clist[[j]],rlist[[j]][total[i]],K,alphaList[[j]],betaList[[j]][total[i]])
        local_den[[j]] = ecdf(tmp)
    }
    den[[i]] = local_den
}


par(mfrow = c(3,3),mai = c(0.4, 0.3, 0.1, 0.1))

for(i in 1:9){
    
    x = data_counts[total[i],]
    
    im = values_to_impute(x)
    
    if(i %in% c(1,4,7)){
        plot(den_o[[i]],lwd = 2, main = "", col = "green",do.points=F, ylab = "" ,xlab = "",xlim = c(0,max(x)))
    }else{
        plot(den_o[[i]],lwd = 2, main = "", col = "green",do.points=F, ylab = "",yaxt='n' ,xlab = "",xlim = c(0,max(x)))
    }
    legend("bottomright",genenames[total[i]])
    
    
    
    
    start = 10
    for(j in 1:10){
        lines(den[[i]][[j]], lwd = 0.5, col = alpha(t_col("pink",percent = start),0.5),do.points=F)
        
        if(j == 1){
            meanCDF = den[[i]][[j]](im)
        }else{
            meanCDF = meanCDF + den[[i]][[j]](im)
        }
        
    }
    
    meanCDF = meanCDF / 10

    lines(x = im, y = meanCDF,lwd = 2, col = "pink", type = "l")
    
}
