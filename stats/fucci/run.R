

load("results.RData")


library(scDDboost)
library(SingleCellExperiment)

source("scDDboost/stats/simulation/eval/eval_deseq2.R")
source("scDDboost/stats/simulation/eval/eval_MAST.R")
source("scDDboost/stats/simulation/eval/eval_scDD.R")


## number of cores to use
ncores = 4

## run scDDboost

D_c = cal_D(data_counts,ncores)
newPDD = PDD(data = data_counts, cd = cd, ncores = ncores, D = D_c, norm = T)

## run MAST

res_mast = eval_MAST(data_counts,cd)

## run DESeq2

res_des = eval_DESeq2(data_counts,cd,T,ncores)

## run scDD, set number of permutations as 100

res_scdd = eval_scDD(data_counts,cd,ncores,100)


