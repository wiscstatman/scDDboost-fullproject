library(ggplot2)


load("marginal.RData")

## samples from mixture NB by estimation
fitVal = function(x,clus,r,K){
  mn = tapply(x,clus,mean)
  cn = tapply(rep(1,length(x)),clus,sum)
  p = mn / (r  + mn) 
  fit = c()
  for(i in 1:K){
    fit = c(fit,rnbinom(cn[i],r,prob = 1-  p[i]))
  }
  return(fit)
}

values_to_impute = function(x, impute_resolution = 100){
    return(seq(
    min(x)
    , max(x)
    , length.out = impute_resolution
    ))
}


K = 6

den = list()
den_o = list()
for(i in 1:9){
  x = data_counts[NZ[i],]
  den_o[[i]] = ecdf(x)
  local_den = list()
  for(j in 1:10){
    tmp = fitVal(x,clist[[j]],rlist[[j]][NZ[i]],K)
    local_den[[j]] = ecdf(tmp)
  }
  den[[i]] = local_den
}



pdf("ecdf.pdf")

par(mfrow = c(3,3),mai = c(1, 0.3, 0.1, 0.1))

for(i in 1:9){
  
  x = data_counts[NZ[i],]
  
  im = values_to_impute(x)
  
  if(i %in% c(1,4,7)){
    plot(den_o[[i]],lwd = 2, main = "", col = "green",do.points=F, ylab = "" ,xlab = "",xlim = c(0,max(x)))
  }else{
    plot(den_o[[i]],lwd = 2, main = "", col = "green",do.points=F, ylab = "",yaxt='n' ,xlab = "",xlim = c(0,max(x)))
  }
  legend("bottomright",genenames[NZ[i]])
  
  
  
  
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

