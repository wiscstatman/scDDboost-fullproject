load("sctrim.RData")

library(SummarizedExperiment)
library(scDDboost)
library(SCnorm)
library(SingleCellExperiment)
library(data.table)
##no need to normalize data
##first normalize the data
##cd = rep(1,247) ##treat all cells to be in one condition for normalization
##DataNorm<-SCnorm(Data=data_count,Conditions=cd,PrintProgressPlots=F,FilterCellNum=10,NCores = 4)
##data_all = assays(DataNorm)$Counts

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

##first only consider DD genes between G1 vs G2
label_ = c(which(cells == 'G1'),which(cells == 'G2'))
cd = c(rep(1,length(which(cells == 'G1'))),rep(2,length(which(cells == 'G2'))))
data_count = cbind(G1,G2)
D_c = D_all[label_,label_]
ccl = pam(D_c,K, T)$clustering  
##can view the different composition of subtypes between G1 and G2
table(ccl[which(cells == 'G1')])
table(ccl[which(cells == 'G2')])

##size factor for EBSeq
sz = MedianNorm(data_count)
hp = c(1, rep(1,nrow(data_count))) ##get hyper parameter

pDD8 = PDD(data = data_count, cd = cd, ncores = 1, K = K, D = D_c,
           sz = sz, hp, pat(K)[[1]], 10, random = T, lambda = 0.5, nrandom = 20)
EDDb = which(pDD8 > 0.95)
length(EDDb)
#scDD has 6805 DD genes


#pDD8_no_random = PDD(data = data_count, cd = cd, ncores = 1, K = K, D = D_c,
#                     sz = sz, hp, pat(K)[[1]], 10, random = F, lambda = 1, nrandom = 0)
#length(which(pDD8_no_random > 0.95))
#EDD8 = which(pDD8_no_random > 0.99)
#gcl = 1:nrow(data_count)
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
length(EDDM)
##MAST has 483 DD genes


prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
condition = factor(cd)
X =  SingleCellExperiment(assays = list(normcounts = data_count), 
                          colData = data.frame(condition))

X_scDD <- scDD(X, prior_param=prior_param, testZeroes=T, permutations = F)
RES = scDD::results(X_scDD)
EDD_dz = which(RES$nonzero.pvalue.adj < 0.05)
EDD_sc = which(RES$combined.pvalue < 0.05)
#DDc = RES$DDcategory
#table(DDc)
length(EDD_sc)
##scDD has 5023 DD genes. 



un1 = union(EDDM,EDD_sc)
scb_uni = setdiff(EDDb, un1)
un2 = union(EDDM, EDDb)
sc_uni = setdiff(EDD_sc,un2)
un3 = union(EDDb,EDD_sc)
m_uni = setdiff(EDDM,un3)



tran = list()
J = 6
tmpp = sample(scb_uni,J)
#tmp = scb_uni[tmpp]
for(i in 1:J){
  tran[[i]] = data_count[tmpp[i],]
}


cur_rn = rn[tmpp]

cond_ind = c(rep("1", length(which(cd ==1))),
             rep("2", length(which(cd==2))))

cur_rn = factor(cur_rn, levels = cur_rn)
df = data.frame(x = rep(cond_ind, J), y = log(do.call(c, tran) + 1), z = rep(cur_rn, each = length(cd)))


#pdf("density_G48_dd.pdf")
pdf("DD_by_scb.pdf")
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
dev.off()

###plot of pDD under different number of subtypes
pdf("sub7_vs_sub8.pdf")
plot(pDD7,pDD8)
dev.off()




