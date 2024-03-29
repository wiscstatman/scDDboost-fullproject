#####run DESeq2
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 

suppressPackageStartupMessages(library(DESeq2))

eval_DESeq2 = function(data_counts, cd, parallel, ncores){
  cd = factor(cd)
  dds <- DESeqDataSetFromMatrix(countData = round(data_counts), 
                                colData = data.frame(condition = cd), 
                                design = ~condition)
  dds = estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq(dds, parallel = parallel, BPPARAM = MulticoreParam(workers=ncores))
  res <- DESeq2::results(dds, contrast = c("condition", levels(factor(cd))[1], 
                                           levels(factor(cd))[2]), alpha = 0.05)
  
  return(res)
}
