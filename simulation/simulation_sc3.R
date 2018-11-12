library(SingleCellExperiment)
library(scDDboost)
library(SC3)
library(scDD)
library(MAST)
library(DESeq2)



filename = paste0("ano_sim",i)
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




loc = -2
scale = 2
split_data(-2,2,2)
split_data = function(loc,scale,n){
  for(II in 1:n){
    
    K = 7
    if(loc < 0){
      name_vec = c("sim","_","neg",abs(loc),'_',scale,'_',II,".rds")
      filename = paste(name_vec,collapse = '')
    }
    else{
      name_vec = c("sim","_",abs(loc),'_',scale,'_',II,".rds")
      filename = paste(name_vec,collapse = '')
    }
    
    # filename = paste0("sim",I)
    readDIR = paste0("/ua/xiuyu/simulation_data/",filename)
    
    dat = readRDS(readDIR)
    
    
    
    data_counts_all = assays(dat)$counts
    
    group = colData(dat)$Group
    
    ###make group1 belong to condition 1 and group 7 belong to condition2
    g11 = which(group == 1)
    
    g22 = which(group == 7)
    
    rest = setdiff(1:400,c(g11,g22))
    
    
    gp1 = sample(rest, 100 - length(g11))
    
    gp2 = setdiff(rest, gp1)
    
    g1 = c(g11,gp1)
    
    g2 = c(g22,gp2)
    
    label1 = c(group[g11],group[gp1])
    
    label2 = c(group[g22],group[gp2])
    
    tmp = c()
    for(i in 1:K){
      a = paste0("rowData(dat)$DEFacGroup",i)
      tmp = union(tmp,which(eval(parse(text = a)) != 1))
    }
    
    DD = tmp
    ED = setdiff(1:nrow(data_counts_all),DD)
    
    
    # g1 = sample(1:ncol(data_counts_all),200)
    # g2 = setdiff(1:ncol(data_counts_all),g1)
    
    data_counts = data_counts_all[,c(g1,g2)]
    cd = c(rep(1,200),rep(2,200))
    
    if(loc < 0){
      save_vec = c("simdata","_","neg",abs(loc),'_',scale,'_',II,".RData")
      savename = paste(save_vec,collapse = '')
    }
    else{
      save_vec = c("simdata","_",abs(loc),'_',scale,'_',II,".RData")
      savename = paste(save_vec,collapse = '')
    }
    # savefile_ = paste0("simdata",I)
    # savefile = paste0(savefile_,".RData")
    save_dir = paste0("/ua/xiuyu/simulation_data/",savename)
    save(dat,data_counts,cd,g1,g2,group,DD,ED,label1,label2, file = save_dir)
    
  }
}

###test for significance
D_c = cal_D(data_counts, 2)
K = 7

sz = MedianNorm(dat)
hp = c(1, rep(1,nrow(dat))) ##get hyper parameter
pDD7 = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
           sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 0.5, nrandom = 50)

EDDb = which(pDD7 > 0.95)
length(EDDb)



run_de(-0.1,0.3,2)

run_de = function(loc,scale,n){
  for(II in 1:n){
    if(loc < 0){
      save_vec = c("simdata","_","neg",abs(loc),'_',scale,'_',II,".RData")
      filename = paste(save_vec,collapse = '')
    }
    else{
      save_vec = c("simdata","_",abs(loc),'_',scale,'_',II,".RData")
      filename = paste(save_vec,collapse = '')
    }
    
    readDIR = paste0("/ua/xiuyu/simulation_data/",filename)
    
    load(readDIR)
    
    mast = eval_MAST(data_counts, cd)
    
    des = eval_DESeq2(data_counts, cd)
    
    scddres = eval_scDD(data_counts, cd)
     
     
    D_c = cal_D(data_counts, 10)
    K = 7
    s = 1
    ref = list()
    ref[[K]] = g_ref(pat(K)[[1]])
    hp = c(s, rep(1,nrow(data_counts)))
    sz = MedianNorm(data_counts)
    sz = rep(1, ncol(data_counts))
    pDD = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
              sz = sz, hp, pat(K)[[1]], 1, random = T, lambda = 5, nrandom = 30)
              
    pDD_nr = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
              sz = sz, hp, pat(K)[[1]], 1, random = F, lambda = 1, nrandom = 0)
    
    ###sc3 random dist
    X <- SingleCellExperiment(assays = list(normcounts = data_counts), 
                              colData = colnames(data_counts))
    counts(X) <- normcounts(X)
    logcounts(X) <- log2(normcounts(X) + 1)
    ##prepare for further sc3 reladted cal
    rowData(X)$feature_symbol = rownames(data_counts)
    X <- sc3_prepare(X)
    sce <- sc3_estimate_k(X)
    K = metadata(sce)$sc3$k_estimation
    K = 7
    dst.sc3 = sc3_calc_dists(X)
    tran.D = sc3_calc_transfs(dst.sc3)
    pre_output = sc3_kmeans(tran.D,K)
    consen.D = sc3_calc_consens(pre_output)
    print(names(metadata(consen.D)$sc3))
    consen_matrix = metadata(consen.D)$sc3$consensus[[1]]$consensus
  
    pDD_sc3 = PDD(data_counts,cd, ncores = 2, K, consen_matrix, sz, hp, pat(K)[[1]],iter = 1, random = T, lambda = 0.2 * mean(consen_matrix), nrandom = 30)
    
    if(loc < 0){
      save_vec = c("res","_","neg",abs(loc),'_',scale,'_',II,".RData")
      savename = paste(save_vec,collapse = '')
    }
    else{
      save_vec = c("res","_",abs(loc),'_',scale,'_',II,".RData")
      savename = paste(save_vec,collapse = '')
    }
    
    save(file = paste0("/ua/xiuyu/simulation_data/",savename), pDD,pDD_sc3,pDD_nr,DD,ED,mast,des,scddres)
    
  }
  
  
}



####
X <- SingleCellExperiment(assays = list(normcounts = data_counts),
colData = colnames(data_counts))
counts(X) <- normcounts(X)
logcounts(X) <- log2(normcounts(X) + 1)
##prepare for further sc3 reladted cal
rowData(X)$feature_symbol = rownames(data_counts)
X <- sc3_prepare(X)
sce <- sc3_estimate_k(X)
K = metadata(sce)$sc3$k_estimation
K = 7
dst.sc3 = sc3_calc_dists(X)
tran.D = sc3_calc_transfs(dst.sc3)
pre_output = sc3_kmeans(tran.D,K)
consen.D = sc3_calc_consens(pre_output)
print(names(metadata(consen.D)$sc3))
consen_matrix = metadata(consen.D)$sc3$consensus[[1]]$consensus
###
pDD_sc3_nr = PDD(data_counts,cd, ncores = 4, K, consen_matrix, sz, hp, pat(K)[[1]],iter = 1, random = F, lambda = mean(consen_matrix), nrandom = 1)




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

pDDM = mast[, "hurdle", "Pr(>Chisq)"]
EDDM = which(pDDM < 0.05)

library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
nsim = 20


nsim = 2
results <- foreach(i=1:nsim) %dopar% {
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

for(i in 1:nsim){
  filename_ = paste0("mast_sim",i)
  filename = paste0(filename_,".RData")
  mast = results[[i]]
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
  cd = factor(cd)
  dds <- DESeqDataSetFromMatrix(countData = round(data_counts), 
                                colData = data.frame(condition = cd), 
                                design = ~condition)
  dds <- DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("condition", levels(factor(cd))[1], 
                                   levels(factor(cd))[2]), alpha = 0.05)

  return(res)
}


results <- foreach(i=2:nsim) %dopar% {
  library(DESeq2)
  library(SingleCellExperiment)
  filename_ = paste0("simdata",i)
  filename = paste0(filename_,".RData")
  load(filename)
  eval_DESeq2(data_counts, cd)
}

for(i in 2:nsim){
  filename_ = paste0("des_sim",i)
  filename = paste0(filename_,".RData")
  mast = results[[i]]
  save(mast,file = filename)
}
save(dds,file = filename)
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
  # EDD_sc = which(RES$nonzero.pvalue.adj < 0.05)
  return(RES)
  # return(RES$nonzero.pvalue.adj)
}

nsim = 20
results <- foreach(i=1:nsim) %dopar% {
  library(SingleCellExperiment)
  library(scDD)
  filename_ = paste0("simdata",i)
  filename = paste0(filename_,".RData")
  load(filename)
  eval_scDD(data_counts, cd)
}



for(i in 1:nsim) {
  library(SingleCellExperiment)
  library(scDD)
  filename_ = paste0("simdata",i)
  filename = paste0(filename_,".RData")
  load(filename)
  scd = eval_scDD(data_counts, cd)
  filename_ = paste0("scDD_sim",i)
  filename = paste0(filename_,".RData")
  save(scd,file = filename)
}

for(i in 1:nsim){
  filename_ = paste0("scDD_sim",i)
  filename = paste0(filename_,".RData")
  scd = results[[i]]
  save(scd,file = filename)
}


ref = list()
for(K in 2:8){
  ref[[K]] = g_ref(pat(K)[[1]])
}
for(i in 2:20){
  filename_ = paste0("simdata",i)
  filename = paste0(filename_,".RData")
  load(filename)
  D_c = cal_D(data_counts, 10)
  K = 7
  s = 1
  hp = c(s, rep(1,nrow(data_counts)))
  sz = MedianNorm(data_counts)
  pDD = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
             sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 1, nrandom = 20)
  filename_ = paste0("scDDboost",i)
  filename = paste0(filename_,".RData")
  save(file = filename,pDD)
}




load("sctrim.RData")

library(SummarizedExperiment)
library(scDDboost)
library(SCnorm)
library(SingleCellExperiment)
library(data.table)

data_all = round(edat)
data_all = as.matrix(data_all)
##generate distance matrix
D_all = cal_D(data_all, 10)

##data according to different phase
G1 = data_all[,which(cells == 'G1')]
G2 = data_all[,which(cells == 'G2')]
S = data_all[,which(cells== 'S')]



##8 subtypes
K = 8
Posp = pat(K)[[1]] #possible partition for 8 groups

##generate refinement relation between partitions, K range from 2 to 9
ref = list()
for(i in 2:9){
    ref[[i]] = g_ref(pat(i)[[1]])
}




label_ = c(which(cells == 'G1'),which(cells == 'G2'))
cd = c(rep(1,length(which(cells == 'G1'))),rep(2,length(which(cells == 'G2'))))
data_counts = cbind(G1,G2)
D_c = D_all[label_,label_]
ccl = pam(D_c,K, T)$clustering
##can view the different composition of subtypes between G1 and G2
table(ccl[which(cells == 'G1')])
table(ccl[which(cells == 'G2')])

dat_mean = rowSums(data_counts)
rm_gene = which(dat_mean == 0)

selected = which(dat_mean > 0)
sel = data_counts[selected,]


sz = MedianNorm(data_counts)
hp = c(1, rep(1,nrow(data_counts))) ##get hyper parameter

pDD = PDD(data = sel, cd = cd, ncores = 10, K = K, D = D_c,
sz = sz, hp, pat(K)[[1]], 1, random = T, lambda = 0.01, nrandom = 50)


pDD = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 0.5, nrandom = 20)

pDD = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
sz = sz, hp, pat(K)[[1]], 10, random = F, lambda = 0.0000, nrandom = 1)

EDDb = which(pDD > 0.95)
length(EDDb)

K = 2


hp = c(0.378, rep(1.41,nrow(sel)))
res2 = EBS(sel,ccl,1:nrow(sel),sz,iter = 0,hp,pat(K)[[1]])
DE = res2$DEpattern
uni = DE[,2]


res = EBTest(data_counts, NgVector = NULL, cd, sz, 1,
Pool = F, NumBin = 1000, ApproxVal = 10^-10, Alpha = NULL,
Beta = NULL, PInput = NULL, RInput = NULL,
PoolLower = .25, PoolUpper = .75, Print = T, Qtrm = 1,QtrmCut=0)

res1 = EBTest(sel, NgVector = NULL, cd, sz, 100,
Pool = F, NumBin = 1000, ApproxVal = 10^-10, Alpha = 0.378,
Beta = 1.41, PInput = NULL, RInput = NULL,
PoolLower = .25, PoolUpper = .75, Print = T, Qtrm = 1,QtrmCut=0)

PP = res1$PPDE
cor(PP[which(!is.na(PP))], pDDM[which(!is.na(PP))])


cor(PP[intersect(which(!is.na(PP)), which(!is.na(EDD_sc)))], EDD_sc[intersect(which(!is.na(PP)), which(!is.na(EDD_sc)))])


loc = c(-3,-2,-0.3,-0.1)
scale = c(2,2,2,0.3)
for(i in 1:4){
    update_scDDboost(loc[i], scale[i], 2)
}

update_scDDboost(-0.1,0.3,2)
update_scDDboost = function(loc,scale,n){
    for(II in 1:n){
        if(loc < 0){
            save_vec = c("simdata","_","neg",abs(loc),'_',scale,'_',II,".RData")
            filename = paste(save_vec,collapse = '')
        }
        else{
            save_vec = c("simdata","_",abs(loc),'_',scale,'_',II,".RData")
            filename = paste(save_vec,collapse = '')
        }
        
        readDIR = paste0("/ua/xiuyu/simulation_data/",filename)
        
        load(readDIR)
        if(loc < 0){
            save_vec = c("res","_","neg",abs(loc),'_',scale,'_',II,".RData")
            savename = paste(save_vec,collapse = '')
        }
        else{
            save_vec = c("res","_",abs(loc),'_',scale,'_',II,".RData")
            savename = paste(save_vec,collapse = '')
        }
        
        load(file = paste0("/ua/xiuyu/simulation_data/",savename))
        
        D_c = cal_D(data_counts, 10)
        K = 7
        s = 1
        ref = list()
        ref[[K]] = g_ref(pat(K)[[1]])
        hp = c(s, rep(1,nrow(data_counts)))
        sz = MedianNorm(data_counts)
        sz = rep(1, ncol(data_counts))
        pDD = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
        sz = sz, hp, pat(K)[[1]], 1, random = T, lambda = 0.01, nrandom = 30)
        
        pDD_nr = PDD(data = data_counts, cd = cd, ncores = 10, K = K, D = D_c,
        sz = sz, hp, pat(K)[[1]], 1, random = F, lambda = 1, nrandom = 0)
        
        ###sc3 random dist
        X <- SingleCellExperiment(assays = list(normcounts = data_counts),
        colData = colnames(data_counts))
        counts(X) <- normcounts(X)
        logcounts(X) <- log2(normcounts(X) + 1)
        ##prepare for further sc3 reladted cal
        rowData(X)$feature_symbol = rownames(data_counts)
        X <- sc3_prepare(X)
        sce <- sc3_estimate_k(X)
        K = metadata(sce)$sc3$k_estimation
        K = 7
        dst.sc3 = sc3_calc_dists(X)
        tran.D = sc3_calc_transfs(dst.sc3)
        pre_output = sc3_kmeans(tran.D,K)
        consen.D = sc3_calc_consens(pre_output)
        print(names(metadata(consen.D)$sc3))
        consen_matrix = metadata(consen.D)$sc3$consensus[[1]]$consensus
        
        pDD_sc3 = PDD(data_counts,cd, ncores = 10, K, consen_matrix, sz, hp, pat(K)[[1]],iter = 1, random = T, lambda = 0.01, nrandom = 30)
        
        
        
        
        save(file = paste0("/ua/xiuyu/simulation_data/",savename), pDD,pDD_sc3,pDD_nr,DD,ED,mast,des,scddres)
        
    }
    
    
}












