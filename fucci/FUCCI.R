load("sctrim.RData")

library(SummarizedExperiment)
library(scDDboost)
library(SCnorm)
library(SingleCellExperiment)
##first normalize the data
cd = rep(1,247) ##treat all cells to be in one condition for normalization
DataNorm<-SCnorm(Data=data_count,Conditions=cd,PrintProgressPlots=F,FilterCellNum=10,NCores = 4)
data_all = assays(DataNorm)$Counts

##generate distance matrix
D_all = cal_D(data_all, 10)
##8 subtypes
K = 8
Posp = pat(K)[[1]] #possible partition for 7 groups


##first only consider DD genes between G1 vs G2
ccl = pam(D_all,K, T)$clustering[1:167]  ##get clustering
data_count = data_all[,1:167]
D_c = D_all[1:167,1:167]
ccl ##can view the different composition of subtypes between G1 and G2
sz = MedianNorm(data_count)
hp = c(1, rep(1,nrow(data_count))) ##get hyper parameter
pDD7 = PDD(data = data_count, cd = cd, ncores = 1, K = K, D = D_c,
           sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 1, nrandom = 100)
EDDb = which(pDD7 > 0.95)
length(EDDb) ##having 252 DD genes

pDD7_no_random = PDD(data = data_count, cd = cd, ncores = 1, K = K, D = D_c,
                     sz = sz, hp, pat(K)[[1]], 10, random = F, lambda = 1, nrandom = 0)
length(which(pDD7_no_random > 0.95))
##having 1352 DD genes without randomization by set



##MAST
library(MAST)
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)
ngeneon=apply(data_count,2,function(x) length(which(x>0)))
cd_=cd
cdata=data.frame(conditions=cd_,ngeneon=1:length(cd))
fdata=data.frame(gene=factor(1:nrow(data_count)))
MNZ10=FromMatrix(log(data_count + 1),cdata,fdata)

cond<-factor(colData(MNZ10)$conditions)
colData(MNZ10)$condition<-cond
cdr2 <-colSums(assay(MNZ10)>0)
colData(MNZ10)$cngeneson <- scale(cdr2)
zlmCond <- zlm(~condition + cngeneson, MNZ10)
summaryCond <- summary(zlmCond, doLRT='condition2') 


summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[contrast=='condition2' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryDt[contrast=='condition2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') 

fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
fcHurdleSig <- merge(fcHurdle[fdr<.05 & abs(coef)> 0], as.data.table(mcols(MNZ10)), by='primerid')

EDDM = as.numeric(fcHurdleSig$gene)

##MAST has 504 DD genes
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
condition = factor(cd)
X =  SingleCellExperiment(assays = list(normcounts = data_count), 
                          colData = data.frame(condition))

X_scDD <- scDD(X, prior_param=prior_param, testZeroes=T)
RES = scDD::results(X_scDD)
EDD_dz = which(RES$nonzero.pvalue.adj < 0.05)
#DDc = RES$DDcategory
#table(DDc)
##scDD has 1575 DD genes. 
