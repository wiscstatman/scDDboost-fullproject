#####run DESeq2 no beta prior
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 

suppressPackageStartupMessages(library(DESeq2))

eval_DESeq2 = function(data_counts, cd){
  cd = factor(cd)
  dds <- DESeqDataSetFromMatrix(countData = round(data_counts), 
                                colData = data.frame(condition = cd), 
                                design = ~condition)
  dds <- DESeq(dds,betaPrior = False)
  res <- DESeq2::results(dds, contrast = c("condition", levels(factor(cd))[1], 
                                           levels(factor(cd))[2]), alpha = 0.05)
  
  return(res)
}
