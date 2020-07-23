## test validity of NB fit 




################################################################## detail of computing, can be skipped

library(scDDboost)
library(MASS)
library(DESeq2)


### data (EMTAB2805)

load("NBfit1.RData")

## ks-test p values
kpval = rep(0,nrow(data_counts))
for(i in 1:nrow(data_counts)){
    x = data_counts[i,which(cd == 1)]
    y = data_counts[i,which(cd == 2)]
    fit = ks.test(x,y)
    kpval[i] = fit$p.value
}

## fdr correction
k_pval = p.adjust(kpval,"fdr")










