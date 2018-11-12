##run scDDboost
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 
#!param ncores number of cores for parallel computation of dist
#!param K number of clusters
#!param nrd number of how many random dist are generated.

suppressPackageStartupMessages(library(scDDboost))

eval_scDDboost <- function(data_counts, cd, ncores, K, nrd = 30){
  #distance matrix
  D_c = cal_D(data_counts, ncores)
  #refinement relation
  ref = list()
  ref[[K]] = g_ref(pat(K)[[1]])
  #default hyper parameters
  s = 1
  hp = c(s, rep(1,nrow(data_counts)))
  
  #default size factor
  sz = rep(1, ncol(data_counts))
  ##posterior based on random dist
  pDD = PDD(data = data_counts, cd = cd, ncores = ncores, K = K, D = D_c,
            sz = sz, hp, pat(K)[[1]], 1, random = T, 
            lambda = mean(D_c), nrandom = nrd)
  ##posterior based on non random dist
  pDD_nr = PDD(data = data_counts, cd = cd, ncores = 1, K = K, D = D_c,
               sz = sz, hp, pat(K)[[1]], 1, random = F,
               lambda = 1, nrandom = 0)
  
  res = list()
  res$RPDD = pDD
  res$PDD = pDD_nr
  return(res)
  
  
}
