load("sctrim.RData")

library(SummarizedExperiment)
library(scDDboost)
library(SCnorm)
library(SingleCellExperiment)
##first normalize the data
cd = rep(1,247) ##treat all cells to be in one condition for normalization
DataNorm<-SCnorm(Data=data_count,Conditions=cd,PrintProgressPlots=F,FilterCellNum=10,NCores = 4)
data_all = assays(DataNorm)$Counts
data_all = round(edat)
data_all = as.matrix(data_all)
##generate distance matrix
D_all = cal_D(data_all, 10)
##7 subtypes
K = 9
Posp = pat(K)[[1]] #possible partition for 7 groups

cd = c(rep(1,91),rep(2,80))
##generate refinement relation between partitions
ref[[K]] = g_ref(Posp)
ref = list()
for(i in 2:8){
  ref[[i]] = g_ref(pat(i)[[1]])
}

##first only consider DD genes between G1 vs G2
ccl = pam(D_all,K, T)$clustering[c(1:91,168:247)]  ##get clustering
#ccl = pam(D_c,K,T)$clustering
G1 = data_all[,1:91]
G2 = data_all[,92:167]
S = data_all[,168:247]
data_count = cbind(G1,S)
D_c = D_all[c(1:91,168:247),c(1:91,168:247)]
ccl ##can view the different composition of subtypes between G1 and G2
sz = MedianNorm(data_count)
#sz = rep(1,ncol(data_count))
hp = c(1, rep(1,nrow(data_count))) ##get hyper parameter
pDD7= PDD(data = data_count, cd = cd, ncores = 1, K = K, D = D_c,
           sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 1, nrandom = 10)
EDDb = which(pDD7 > 0.95)
length(EDDb)

pDD9_no_random = PDD(data = data_count, cd = cd, ncores = 1, K = K, D = D_c,
                     sz = sz, hp, pat(K)[[1]], 10, random = F, lambda = 1, nrandom = 0)
length(which(pDD9_no_random > 0.99))
EDD9 = which(pDD9_no_random > 0.99)
gcl = 1:nrow(data_count)
##having 1352 DD genes without randomization by set
ebres = EBS(data_count,ccl,gcl,sz,10,hp,Posp)
DE = ebres$DEpattern

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
fcHurdleSig <- merge(fcHurdle[fdr<.01 & abs(coef)> 0], as.data.table(mcols(MNZ10)), by='primerid')

EDDM = as.numeric(fcHurdleSig$gene)
length(EDDM)
##MAST has 504 DD genes
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
condition = factor(cd)
X =  SingleCellExperiment(assays = list(normcounts = data_count), 
                          colData = data.frame(condition))

X_scDD <- scDD(X, prior_param=prior_param, testZeroes=T, permutations = F)
RES = scDD::results(X_scDD)
EDD_dz = which(RES$nonzero.pvalue.adj < 0.01)
EDD_sc = which(RES$combined.pvalue < 0.01)
#DDc = RES$DDcategory
#table(DDc)
##scDD has 1575 DD genes. 
un1 = union(EDDM,EDD_sc)
scb_uni = setdiff(EDD9, un1)

tran = list()
J = 6
tmpp = sample(EDD9,J)
tmp = scb_uni[1:J]
for(i in 1:J){
  tran[[i]] = data_count[tmp,]
}


cur_rn = rn[tmp]

cond_ind = c(rep("1", length(which(cd ==1))),
             rep("2", length(which(cd==2))))

cur_rn = factor(cur_rn, levels = cur_rn)
df = data.frame(x = rep(cond_ind, J), y = log(do.call(c, tran) + 1), z = rep(cur_rn, each = length(cd)))


#pdf("density_G48_dd.pdf")

pp<-ggplot(df,aes(factor(x),y))+ geom_violin(aes(colour = factor(x))) + geom_point(size = 0.1,position = position_jitter(w = 0.05, h = 0)) +
  xlab("conditions") + 
  ylab("gene expressions")+
  theme(panel.background = element_rect(
    fill = 'white', colour = 'black'),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    legend.title=element_blank(),
    panel.grid.minor.x = element_line(size = 0.5),
    panel.grid.minor.y = element_line(size = 0.5),
    panel.grid.major.x = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.5),
    panel.grid.major = element_line(colour = "grey"))
pp + facet_wrap( ~ z, ncol = 2)
#dev.off()






