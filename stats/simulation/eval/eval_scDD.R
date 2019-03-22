#####run scDD
#!param data_counts is a matrix of single cell counts
#!param cd is condition label

suppressPackageStartupMessages(library(scDD))
suppressPackageStartupMessages(library(SingleCellExperiment))

eval_scDD = function(data_counts, cd, ncores, npermu){
  
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  condition = factor(cd)
  X =  SingleCellExperiment(assays = list(normcounts = data_counts), 
                            colData = data.frame(condition))
  
  
  X_scDD <- scDD(X, prior_param=prior_param, testZeroes=T, param = BiocParallel::MulticoreParam(workers = ncores), permutations=npermu)
  RES = scDD::results(X_scDD)
  # EDD_sc = which(RES$nonzero.pvalue.adj < 0.05)
  return(RES)
  # return(RES$nonzero.pvalue.adj)
}
