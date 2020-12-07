

##function to split data to 2 conditions used for FDR analysis

split_data = function(data_counts, cd, saveName){
  NL = which(cd == 1)
  n = length(NL)
  tmp = data_counts[,NL]
  G1 = sample(1:n,round(n / 2))
  G2 = setdiff(1:n,G1)
  data1 = tmp[,G1]
  data2 = tmp[,G2]
  data_null = cbind(data1, data2)
  cd_null = c(rep(1, length(G1)), rep(2, length(G2)))
  save(file = saveName, data_null, cd_null)
}

##data sets used to consider
loadName = c("GSE79102.RData", "GSE48968.RData", "TASIC.RData", "GSE94383.RData",
             "DEC_EC.RData", "Fucci.RData", "GSE60749.RData","GSE74596.RData",
             "EMTAB2805.RData","GSE63818.RData")

##
for(XX in loadName){
  load(paste0("empirical/",XX))
  for(j in 1:5)
  {
    if(j == 1){
      saveName = paste0("FDR/", XX)
    }
    else{
    saveName = paste0("FDR", (j - 1))
    saveName = paste0(saveName,"/")
    saveName = paste0(saveName, XX)
    }
    split_data(data_counts, cd, saveName)
  }
}


##after we have the data
#import function for evaluation
source("~/Desktop/scDDboost/stats/simulation/eval/eval_deseq2.R")
source("~/Desktop/scDDboost/stats/simulation/eval/eval_MAST.R")
source("~/Desktop/scDDboost/stats/simulation/eval/eval_scDD.R")
source("~/Desktop/scDDboost/stats/simulation/eval/eval_scDDboost.R")

for(XX in loadName){
  for(i in 1:5){
    if(i > 1)
      TMP = paste0("FDR", (i -1))
    else
      TMP = "FDR"
    TMP = paste0(TMP,"/")
    load(paste0(TMP,XX))
    res_mast = eval_MAST(data_null, cd_null)
    res_des = eval_DESeq2(data_null,cd_null)
    res_scdd = eval_scDD(data_null,cd_null, 2)
    res_scddboost = eval_scDDboost(data_counts = data_null, cd = cd_null, ncores = 2)
    #res_scddboost_sc3 = eval_scDDboost_sc3(data_counts = data_null, cd = cd_null,U = 100, ncores = 2, K = 5, niter = 20,nrd = 30)
    sD = paste0(TMP,"res_")
    save(file = paste0(sD,XX),res_scddboost, res_scdd, res_des, res_mast)
  }
}


