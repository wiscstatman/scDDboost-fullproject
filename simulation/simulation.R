library(SingleCellExperiment)

filename = paste0("sim",i)
readDIR_ = paste0("/ua/xiuyu/simulation_data/",filename)
readDIR = paste0(readDIR_,".rds")

dat = readRDS(readDIR)

data_counts_all = assays(dat)$counts

group = colData(dat)$Group

gp = sort(unique(group))



names(rowData(dat))


rowData(dat)$GeneMean[1:10]

length(rowData(dat)$BaseGeneMean)

rowData(dat)$BaseGeneMean[1:10]


###get index of DD genes 
tmp = c()
for(i in 1:K){
  a = paste0("rowData(dat)$DEFacGroup",i)
  tmp = union(tmp,which(eval(parse(text = a)) != 1))
}

DD = tmp
ED = setdiff(1:nrow(data_counts),DD)

###split into two conditions
g1 = sample(1:ncol(data_counts),200)
g2 = setdiff(1:ncol(data_counts),400)

data_counts = data_counts_all[,c(g1,g2)]
cd = c(rep(1,200),rep(2,200))


savefile_ = paste0("simdata",i)
savefile = paste0(savefile_,".rds")

saveRDS(dat,data_counts,cd,g1,g2,group,DD,ED, file = savefile)



for(I in 2:20){
  
  K = 7
  filename = paste0("sim",I)
  readDIR_ = paste0("/ua/xiuyu/simulation_data/",filename)
  readDIR = paste0(readDIR_,".rds")
  
  dat = readRDS(readDIR)
  
  
  data_counts_all = assays(dat)$counts
  
  group = colData(dat)$Group
  
  
  tmp = c()
  for(i in 1:K){
    a = paste0("rowData(dat)$DEFacGroup",i)
    tmp = union(tmp,which(eval(parse(text = a)) != 1))
  }
  
  DD = tmp
  ED = setdiff(1:nrow(data_counts_all),DD)
  
  
  g1 = sample(1:ncol(data_counts_all),200)
  g2 = setdiff(1:ncol(data_counts_all),g1)
  
  data_counts = data_counts_all[,c(g1,g2)]
  cd = c(rep(1,200),rep(2,200))
  
  
  savefile_ = paste0("simdata",I)
  savefile = paste0(savefile_,".RData")
  
  save(dat,data_counts,cd,g1,g2,group,DD,ED, file = savefile)
  
}


##apply mast to simulation

load("simdata1.RData")


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

library(doParallel)
cl <- makeCluster(10)
registerDoParallel(cl)
nsim = 20


nsim = 2
results <- foreach(i=3:nsim) %dopar% {
  library(MAST)
  filename_ = paste0("simdata",i)
  filename = paste0(filename_,".RData")
  load(filename)
  eval_MAST(data_counts, cd)
}



summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='condition2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='condition2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 
stopCluster(cl)

for(i in 3:nsim){
  filename_ = paste0("mast_sim",i)
  filename = paste0(filename_,".RData")
  mast = results[[i - 2]]
  save(mast,file = filename)
}

length(mast[, "hurdle", "Pr(>Chisq)"])
###for real using edger get 

grp <- cd
names(grp) = 1:400
dge <- DGEList(counts = data_counts)
dge <- edgeR::calcNormFactors(dge)
cpms <- cpm(dge)
sca <- FromMatrix(exprsArray = log2(cpms + 1), 
                  cData = data.frame(wellKey = names(grp), 
                                     grp = grp))
zlmdata <- zlm.SingleCellAssay(~grp, sca)
mast <- lrTest(zlmdata, "grp")
})

hist(mast[, "hurdle", "Pr(>Chisq)"], 50)






#####run deseq2
eval_DESeq2 = function(data_counts, cd){
  dds <- DESeqDataSetFromMatrix(countData = round(data_counts), 
                                colData = data.frame(condition = cd), 
                                design = ~condition)
  dds <- DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("condition", levels(factor(cd))[1], 
                                   levels(factor(cd))[2]), alpha = 0.05)

  return(res)
}


###run scDD
scDatList <- list()
(groups <- unique(cd))
for (i in 1:length(groups)) {
  scDatList[[paste0("G", i)]] <- as.matrix(data_counts[, which(cd == groups[i])])
}
datNorm.scran <- scDD::preprocess(scDatList, 
                                  ConditionNames = names(scDatList),
                                  zero.thresh = 1, median_norm = TRUE)
condition <- L$condt[colnames(datNorm.scran)]
condition <- as.numeric(as.factor(condition))
names(condition) <- colnames(datNorm.scran)

SDSumExp <- SummarizedExperiment(assays = list("NormCounts" = datNorm.scran),
                                 colData = data.frame(condition))
prior_param <- list(alpha = 0.01, mu0 = 0, s0 = 0.01, a0 = 0.01, b0 = 0.01)
scd <- scDD(SDSumExp, prior_param = prior_param, testZeroes = FALSE,
            param = BiocParallel::MulticoreParam(workers = 1), 
            condition = "condition", min.size = 3, min.nonzero = NULL)
res <- results(scd)



eval_scDD = function(data_counts, cd){

  prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
  condition = factor(cd)
  X =  SingleCellExperiment(assays = list(normcounts = data_counts), 
                            colData = data.frame(condition))
  
  
  X_scDD <- scDD(X, prior_param=prior_param, testZeroes=T)
  RES = scDD::results(X_scDD)
  EDD_sc = which(RES$nonzero.pvalue.adj < 0.05)

  return(RES$nonzero.pvalue.adj)
}





