## test validity of NB fit 




################################################################## detail of computing, can be skipped

library(scDDboost)
library(MASS)
library(DESeq2)

lsz = function (pDD, FDR = 0.01)
{
    ee <- 1 - pDD
    oe <- sort(ee)
    or = order(ee)
    ff <- cumsum(oe)/(1:length(oe))
    return(or[which(ff < FDR)])
}


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



### data (Fucci)

load("Fucci-Data.RData")

## data_counts: the expression matrix
## cd: condition label

## genenames: gene names
## cell_cycle_genes: name for cell cycle related genes

## p_scDDboost: probabilities of DD by scDDboost
## p_scDD, p_mast, p_deseq2, p_t, p_ks: fdr adjusted p values by scDD, MAST, DESeq2, t-test, ks-test.

## get DD genes with 0.05 threshold
scDDb_dd = lsz(p_scDDboost,0.05)

## get the corresponding gene names
scDDb_genes = genenames[scDDb_dd]


## similar for MAST, scDD, DESeq2, t-test and ks-test
mst_dd = which(p_mast < 0.05)
scdd_dd = which(p_scDD < 0.05)
des_dd = which(p_deseq2 < 0.05)
t_dd = which(p_t < 0.05)
ks_dd = which(p_ks < 0.05)

mst_genes = genenames[mst_dd]
scdd_genes = genenames[scdd_dd]
des_genes = genenames[des_dd]
t_genes = genenames[t_dd]
ks_genes = genenames[ks_dd]

## union of DD genes by MAST,scDD,DESeq2 and t-test
U_genes = union(union(union(mst_genes,scdd_genes),des_genes),t_genes)

## uniquely by scDDboost
uni_genes = setdiff(scDDb_genes, U_genes)

## intersection with cell cycle related, 137 uniquely found by scDDboost
sum(scDDb_genes %in% cell_cycle_genes)

## only 2 genes found by k_genes also in those uni_genes by scDDboost

sum(k_genes %in% intersect(uni_genes,cell_cycle_genes))

## which are

k_genes[which(k_genes %in% intersect(uni_genes,cell_cycle_genes) == T)]

## they are not the 3 genes we found support for DD in G1 vs G2/M


## 237 genes find by ks also cell-cycle related

sum(k_genes %in% cell_cycle_genes)

## while scDDboost found 409 cell-cycle related genes
sum(scDDb_genes %in% cell_cycle_genes)

## excluding those union found by MAST,DESeq2,scDD and t-test
## ks test find 2 cell-cycle related genes
sum(setdiff(k_genes,U_genes) %in% cell_cycle_genes)


