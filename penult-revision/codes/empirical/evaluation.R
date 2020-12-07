library(SingleCellExperiment)
library(scDDboost)
library(ggplot2)

source("codesToRunOtherMethods/eval_MAST.R")
source("codesToRunOtherMethods/eval_deseq2.R")
source("codesToRunOtherMethods/eval_scDD.R")

eval_scDDboost <- function(data_counts, cd, ncores,epi = 1,nrd = 50){
  #distance matrix
  D_c = cal_D(data_counts, ncores)
  #refinement relation
  #ref = list()
  #ref[[K]] = g_ref(pat(K)[[1]])
  #default hyper parameters
  #s = 1
  #hp = c(s, rep(1,nrow(data_counts)))
  
  #default size factor
  sz = rep(1, ncol(data_counts))
  ##posterior based on random dist
  pDD = PDD(data = data_counts, cd = cd, ncores = ncores, D = D_c, nrandom = nrd, epi = epi)
  ##posterior based on non random dist
  #pDD_nr = PDD(data = data_counts, cd = cd, ncores = ncores, D = D_c, epi = 0.17,random = F)
  
  res = list()
  res$RPDD = pDD
  #res$PDD = pDD_nr
  res$D = D_c
  return(res)
  
  
}

##
readDir = "user-specified"
load(readDir)

ncores = 4
res_PDD = eval_scDDboost(data_counts,cd,ncores)
res_mast = eval_MAST(data_counts,cd)
res_des = eval_DESeq2(data_counts,cd,T,ncores)
res_scdd = eval_scDD(data_counts,cd,ncores,permu = 100)

threshold = 0.05

## number of DD genes
length(which(res_mast[,3,3] < threshold))
length(which(res_des$padj < threshold))
length(which(res_scdd$combined.pvalue.adj < threshold))

listsize(res_PDD$RPDD,threshold)
