library(ggplot2)


load("marginal.RData")

## data_counts : transcripts matrix
## clist: list of randomized clusters (in total 10 randomizations)
## rlist: list of r estimated by EBSeq
## alphaList: list of hyperparameter alpha estimated by EBSeq
## betaListt: list of hyperparameter beta estimated by EBSeq

## trouble: non-constant shape genes
## total: selected genes




## samples from mixture NB by estimation

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

## for averging ecdf
values_to_impute = function(x, impute_resolution = 100){
    return(seq(
    min(x)
    , max(x)
    , length.out = impute_resolution
    ))
}


## make transparent color

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}


## random sample 6 non zero genes

set.seed(1200)

rM = rowMeans(data_counts)

NZ = sample(which(rM > 0),6)


## random sample 3 non-constant shape genes
names(trouble) = genenames[trouble]
NC = sample(trouble,3)


total = c(NZ,NC)

## 3 violate constant shape
intersect(total,trouble)

## 6 subtypes
K = 6

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


pdf("ecdf.pdf")

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
dev.off()
