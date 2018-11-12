#####run MAST
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 

suppressPackageStartupMessages(library(MAST))


eval_MAST = function(data_counts, cd){
  freq_expressed <- 0.2
  FCTHRESHOLD <- log2(1.5)
  ngeneon=apply(data_counts,2,function(x) length(which(x>0)))
  cd_=cd
  cdata=data.frame(conditions=cd_,ngeneon=1:length(cd))
  fdata=data.frame(gene=factor(1:nrow(data_counts)))
  MNZ10=FromMatrix(log(data_counts + 1),cdata,fdata)
  cond<-factor(colData(MNZ10)$conditions)
  colData(MNZ10)$condition<-cond
  cdr2 <-colSums(assay(MNZ10)>0)
  colData(MNZ10)$cngeneson <- scale(cdr2)
  zlmCond <- zlm(~condition + cngeneson, MNZ10)
  # summaryCond <- summary(zlmCond, doLRT='condition2') 
  mast <- lrTest(zlmCond, "condition")
  return(mast)
}
