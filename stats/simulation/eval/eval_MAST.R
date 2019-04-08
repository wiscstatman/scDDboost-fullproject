#####run MAST
#!param data_counts is a matrix of single cell counts
#!param cd is condition label 

suppressPackageStartupMessages(library(MAST))


eval_MAST = function(data_counts, cd){
    #freq_expressed <- 0.2
    #FCTHRESHOLD <- log2(1.5)
  #ngeneon=apply(data_counts,2,function(x) length(which(x>0)))
  
  cdata=data.frame(conditions=cd,ngeneon=1:length(cd))
  fdata=data.frame(gene=factor(1:nrow(data_counts)))
  
  
  DATA=FromMatrix(log(data_counts + 1),cdata,fdata)
  
  cond<-factor(colData(DATA)$conditions)
  
  colData(DATA)$condition<-cond
  ## cellular detection rate
  cdr2 <-colSums(assay(DATA)>0)
  colData(DATA)$cngeneson <- scale(cdr2)
  
  zlmCond <- zlm(~condition + cngeneson, DATA)
  mast <- lrTest(zlmCond, "condition")
  return(mast)
}
