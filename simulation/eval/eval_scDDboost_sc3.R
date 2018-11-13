##run scDDboost based on sc3 clustering method
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 
#!param ncores number of cores for parallel computation of dist
#!param nrd number of how many random dist are generated.

suppressPackageStartupMessages(library(scDDboost))
suppressPackageStartupMessages(library(SC3))
eval_scDDboost_sc3 <- function(data_counts, cd, ncores, nrd = 30){
  #distance matrix
  X <- SingleCellExperiment(assays = list(normcounts = data_counts), 
                            colData = colnames(data_counts))
  counts(X) <- normcounts(X)
  logcounts(X) <- log2(normcounts(X) + 1)
  rowData(X)$feature_symbol = rownames(data_counts)
  X <- sc3_prepare(X)
  sce <- sc3_estimate_k(X)
  K = metadata(sce)$sc3$k_estimation
  dst.sc3 = sc3_calc_dists(X)
  tran.D = sc3_calc_transfs(dst.sc3)
  pre_output = sc3_kmeans(tran.D,K)
  consen.D = sc3_calc_consens(pre_output)
  print(names(metadata(consen.D)$sc3))
  consen_matrix = metadata(consen.D)$sc3$consensus[[1]]$consensus
  
  #refinement relation
  ref = list()
  ref[[K]] = g_ref(pat(K)[[1]])
  #default hyper parameters
  s = 1
  hp = c(s, rep(1,nrow(data_counts)))
  
  #default size factor
  sz = rep(1, ncol(data_counts))
  ##posterior based on random dist
  pDD_sc3 = PDD(data = data_counts, cd = cd, ncores = ncores, K = K, D = consen_matrix,
            sz = sz, hp, pat(K)[[1]], 1, random = T, 
            lambda = mean(D_c), nrandom = nrd)
  ##posterior based on non random dist
  pDD_nr = PDD(data = data_counts, cd = cd, ncores = 1, K = K, D = consen_matrix,
               sz = sz, hp, pat(K)[[1]], 1, random = F,
               lambda = 1, nrandom = 0)
  
  res = list()
  res$RPDD = pDD_sc3
  res$PDD = pDD_nr
  return(res)
  
}
