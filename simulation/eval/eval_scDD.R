

suppressPackageStartupMessages(library(scDD))

eval_scDD = function(data_counts, cd){
  
  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  condition = factor(cd)
  X =  SingleCellExperiment(assays = list(normcounts = data_counts), 
                            colData = data.frame(condition))
  
  
  X_scDD <- scDD(X, prior_param=prior_param, testZeroes=T)
  RES = scDD::results(X_scDD)
  # EDD_sc = which(RES$nonzero.pvalue.adj < 0.05)
  return(RES)
  # return(RES$nonzero.pvalue.adj)
}