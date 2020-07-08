## test validity of NB fit 




################################################################## detail of computing, can be skipped

library(scDDboost)
library(MASS)
library(DESeq2)

### log likelihood of NB
OBJ = function(r,p,x){
  n = length(x)
  obj = sum(x) * log(p) + n * r * log(1 - p) + sum(lgamma(r + x)) - n * lgamma(r) - sum(lgamma(x + 1))
  
  return(-obj)
}

OBJwrapper = function(param,x){
  r = param[[1]]
  p = param[[2]]
  
  return(OBJ(r,p,x))
}


### goodness of fit test for NB
gof = function(xx){
  x = sample(xx,length(xx),T)
  uniq = length(unique(x))
  s = c(0,0.7,0.85,1)
  ns = length(s) - 1
  if(uniq > 4){
    df = ns
  }else{
    return(NA)
  }
  
  prop = length(which(x < 0.01)) / length(x)
  
  if(prop > 0.7)
  {
    return(NA)
  }    
  
  
  
  tryCatch(
    {
      fit = suppressWarnings(nlminb(start = c(1,0.5),objective = OBJwrapper, x = x, lower = c(0,0)))
      sz = fit$par[1]
      p = fit$par[2]
      Pt = pnbinom(xx,size = sz, prob = 1 - p)
      Ot = rep(0,ns)
      for(i in 1:ns)
      {
        Ot[i] = sum((Pt >= s[i]) * (Pt < s[i + 1]))
      }
      n = length(x)
      Exp = n * diff(s)
      ts = sum((Ot - Exp) ^ 2 / Exp)
      res = 1 - pchisq(ts,df)
      return(res)
    },error = function(cond){
      message(cond)
      return(NA)
    })
  
}


K = max(ccl)

dataList = list()

for(i in 1:K){
  dataList[[i]] = data_counts[,which(ccl == i)]
}

pval = list()
G = nrow(data_counts)
for(i in 1:K){
  pv = rep(0,G)
  for(j in 1:G){
    tmp = dataList[[i]][j,]
    #boot = sample(tmp,length(tmp),T)
    pv[j] = gof(tmp)
  }
  pval[[i]] = pv
}

pvalAll = c()
for(i in 1:K){
  pvalAll = c(pvalAll,pval[[i]])
}

padjAll = p.adjust(pvalAll,"fdr")

padj = list()
for(i in 1:K){
  start = (i - 1) * G + 1
  end = i * G
  padj[[i]] = padjAll[start:end]
}

##  t-test p values

tpval = rep(NA,nrow(data_counts))
for(i in 1:nrow(data_counts)){
  x = data_counts[i,which(cd == 1)]
  y = data_counts[i,which(cd == 2)]
  fit = t.test(x,y)
  tpval[i] = fit$p.value
}

t_pval = p.adjust(tpval,"fdr")


###############################################################################


load("NBfit1.RData")


# pval: list of K vectors, each vector representing unadjusted p values of goodness of fit of NB at each cluster
# padj: list of K vectors, each vector representing fdr adjusted p values of goodness of fit of NB at each cluster
###  p values are adjusted by pooling all cluster specific p values together and do the BH procedure. 

# p_scDDboost: posterior of being DD 
# p_deseq2,p_mast,p_scdd, adjusted p value for deseq2, mast and scDD

# data_counts, transcripts matrix
# cd, condition label
# ccl, clustering 

K = max(ccl)

## cdf plot

den = list()
K = 5
for(i in 1:K){
  den[[i]] = ecdf(padj[[i]])
}
c1 <- rainbow(4)


plot(den[[1]],main = "Goodness of Fit", xlab = "BH adjusted pvalue", ylab = "cumulative mass")
for(i in 2:K){
  lines(den[[i]],col = c1[i - 1])
}


## position may against NB assumption will be labelled 1
nG = nrow(data_counts)
WHICH = rep(0,nG)
for(i in 1:nG){
  flag = 0
  for(j in 1:K){
    tmp = padj[[j]][i]
    if(!is.na(tmp)){
      if(padj[[j]][i] < 0.05){
        flag = 1
      }
    }
  }
  if(flag == 1){
    WHICH[i] = 1
  }
}


p1 = which(p_scDDboost[which(WHICH == 1)] > 0.95)
p2 = which(tpval.adj[which(WHICH == 1)] < 0.05)
length(p1)
length(p2)
length(intersect(p1,p2))







