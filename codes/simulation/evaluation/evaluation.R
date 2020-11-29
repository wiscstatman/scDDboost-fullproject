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

run_de = function(loc,scale,nrep, ncores, DIR, epi = 1){
  for(II in 1:nrep){
    if(loc < 0){
      save_vec = c("sim","_","neg",abs(loc),'_',scale,'_',II,".rds")
      filename = paste(save_vec,collapse = '')
    }
    else{
      save_vec = c("sim","_",abs(loc),'_',scale,'_',II,".rds")
      filename = paste(save_vec,collapse = '')
    }
    
    readDIR = paste0("user-specified",filename)
    
    load(readDIR)
    
    mast = eval_MAST(data_counts, cd)
    
    des = eval_DESeq2(data_counts, cd,T,ncores)
    
    scddres = eval_scDD(data_counts, cd, ncores,0)
     
    #U = 3
     
    scddboost_res = eval_scDDboost(data_counts, cd, ncores, epi = epi)
    
    #scddboost_sc3_res = eval_scDDboost_sc3(data_counts, cd, ncores, U, niter = 20, nrd = 30, K = 7)

    if(loc < 0){
      save_vec = c("res","_","neg",abs(loc),'_',scale,'_',II,".RData")
      savename = paste(save_vec,collapse = '')
    }
    else{
      save_vec = c("res","_",abs(loc),'_',scale,'_',II,".RData")
      savename = paste(save_vec,collapse = '')
    }
    
    save(file = paste0(DIR,savename),
         DD,ED,mast,des,scddres,scddboost_res)
  }
}

