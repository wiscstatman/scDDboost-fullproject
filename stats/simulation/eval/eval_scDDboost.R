##run scDDboost
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 
#!param ncores number of cores for parallel computation of dist
#!param K number of clusters
#!param nrd number of how many random dist are generated.

suppressPackageStartupMessages(library(scDDboost))

eval_scDDboost <- function(data_counts, cd, ncores){
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
  pDD = PDD(data = data_counts, cd = cd, ncores = ncores, D = D_c, epi = 0.17)
  ##posterior based on non random dist
  pDD_nr = PDD(data = data_counts, cd = cd, ncores = ncores, D = D_c, epi = 0.17)
  
  res = list()
  res$RPDD = pDD
  res$PDD = pDD_nr
  res$D = D_c
  return(res)
  
  
}
