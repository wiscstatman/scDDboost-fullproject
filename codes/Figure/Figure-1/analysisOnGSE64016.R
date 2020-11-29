library(SingleCellExperiment)
library(scDDboost)
library(ggplot2)
source("codesToRunOtherMethods/eval_MAST.R")
source("codesToRunOtherMethods/eval_deseq2.R")
source("codesToRunOtherMethods/eval_scDD.R")

## Fucci data set
### data_counts: counts matrix
### cd: condition label (G1 vs G2/M)
load("GSE64016.RData")


## number of cores to use
ncores = 4

## run scDDboost
# calculating distance matrix
D_c = cal_D(data_counts,ncores)

newPDD = PDD(data = data_counts, cd = cd, ncores = ncores, D = D_c, norm = T)

## run MAST
## and adjusted p value from MASTT
res_mast = eval_MAST(data_counts,cd)
mst = res_mast[,"hurdle","Pr(>Chisq)"]

## run DESeq2
## and adjusted p value from DESeq2
res_des = eval_DESeq2(data_counts,cd,T,ncores)
des = res_des$padj

## run scDD, set number of permutations as 100
res_scdd = eval_scDD(data_counts,cd,ncores,100)
scdd = res_scdd$combined.pvalue.adj

## run t-test and wilcox test
tdat <- log2( data_counts + 1/2 )

X <- tdat[,cd==1]
Y <- tdat[,cd==2]

pv.t <- rep(1,nrow(X))
pv.w <- rep(1,nrow(X))

for( i in 1:nrow(X) )
{
    if( var( c(X[i,], Y[i,] ) ) > 0 )
    {
        pv.t[i]  <- ( t.test( X[i,], Y[i,] ) )$p.value
        pv.w[i]  <- ( wilcox.test( X[i,], Y[i,] ) )$p.value
    }
}

qv.t <- p.adjust( pv.t, method="BH" )
qv.w <- p.adjust( pv.w, method="BH" )

threshold  = 0.05

genenames = rownames(data_counts)

des_dd = which(res_des$padj < threshold)
scdd_dd = which(res_scdd$combined.pvalue.adj < threshold)

mast_dd = which(mst < threshold)

scddboost_dd = lsz(newPDD, threshold)

t_dd = which(qv.t < threshold)
w_dd = which(qv.w < threshold)

## DD genes by other methods
others = union(des_dd,union(mast_dd,scdd_dd))
others = union(others,union(t_dd,w_dd))

## uniquely DD by scDDboost
uni = setdiff(scddboost_dd, others)

## DD genes by all methods
all = union(scddboost_dd,others)

## ED genes by all methods
ED = setdiff(1:nrow(data_counts),all)

## others identified but not identified by scDDboost
others_id = setdiff(others,scddboost_dd)

## match with GO annotation
GO_pool = read.csv("GO_pool.csv")
smb = as.character(unique(GO_pool$Symbol))
## convert genenames to upper case for matching genenames from Fucci
smb = sapply(smb,function(x) toupper(x))
scb_genes = genenames[scddboost_dd]
unique_genes = genenames[uni]

## uniquely identified cell cycle annotated genes
uni_cyc = intersect(smb,uni_gene)


############## plot
## we have 3 cell cycle genes related to DD in G1 vs G2/M
tran = list()
cycleGenes = c("BIRC5", "CKAP2", "HMMR")
#cycleGenes = sample(uni_cyc,20)
J = length(cycleGenes)
p_cur = mapToIndex(cycleGenes,genenames)
for(i in 1:J){
  tran[[i]] = log(data_counts[p_cur[i],] + 1)
}
cur_rn = genenames[p_cur]

cond_ind = c(rep("G1", length(which(cd ==1))),
             rep("G2/M", length(which(cd==2))))

cur_rn = factor(cur_rn, levels = cur_rn)
TRAN = do.call(c,tran)
df = data.frame(x = rep(cond_ind, J), y = TRAN, z = rep(cur_rn, each = length(cd)))


pp<-ggplot(df,aes(factor(x),y))+ geom_violin(aes(colour = factor(x)),lwd = 2) + geom_point(size = 0.1,position = position_jitter(w = 0.05, h = 0)) +
  xlab("conditions") +
  ylab("gene expressions")+
  theme(panel.background = element_rect(
    fill = 'white', colour = 'black'),
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.title.y=element_blank(),
    legend.title=element_blank(),
    legend.text = element_text(size = 12, face = "bold"),
    panel.grid.minor.x = element_line(size = 0.5),
    panel.grid.minor.y = element_line(size = 0.5),
    panel.grid.major.x = element_line(size = 0.5),
    panel.grid.major.y = element_line(size = 0.5),
    panel.grid.major = element_line(colour = "grey"))
pp + facet_wrap( ~ z, ncol = 1,scales = "free") + theme(strip.text = element_text(size = 12, face = "bold"))


#### heatmap for overall uniquely identified genes with GO cell cycle annotation
library(RColorBrewer)

winsor = function (x, fraction=.05)
{
  if(length(fraction) != 1 || fraction < 0 ||
     fraction > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}

mapToIndex = function (subGeneNames, genenames)
{
    res = rep(0, length(subGeneNames))
    for (i in 1:length(subGeneNames)) res[i] = which(genenames ==
        subGeneNames[i])
    return(res)
}

p_cur = mapToIndex(uni_cyc,genenames)
ccl = pam(D_c,6,T,cluster.only = T)
ccl_cd1 = ccl[1:91]
ccl_cd2 = ccl[92:167]
cd1_ord = order(ccl_cd1)
cd2_ord = order(ccl_cd2)
D1 = data_counts[p_cur,1:91]
D2 = data_counts[p_cur,92:167]
plotData = cbind(D1[,cd1_ord],D2[,cd2_ord])

my_group=as.numeric(as.factor(substr(c(ccl_cd1[cd1_ord],ccl_cd2[cd2_ord]), 1 , 1)))
my_col=brewer.pal(9, "Set1")[my_group]
PL = 1:length(p_cur)
#PL = sample(1:length(p_cur))
tmp = plotData[PL,]

for(i in 1:length(PL)){
  tmpp = tmp[i,]
  tmp[i,] = winsor(tmpp)
}
mean_tmp = apply(tmp,1,mean)
low_expr = which(mean_tmp < 1)
tmp = tmp[-low_expr,]

i1 = which(rownames(tmp) == "BIRC5")

i2 = which(rownames(tmp) == "HMMR")

i3 = which(rownames(tmp) == "CKAP2")



makeRects <- function(cells){
  nx = 167
  ny = 137
  coords = expand.grid(ny:1, 1:nx)[cells,]
  xl=coords[,2]-0.49
  xl = min(xl)
  
  yb=coords[,1]-0.49
  yb = unique(yb)

  xr=coords[,2]+0.49
  xr = max(xr)
    
  yt=coords[,1]+0.49
  yt = unique(yt)

  rect(xl,yb,xr,yt,border="black",lwd=1)
  abline(v = 91.5, lwd = 3)
}



selection = matrix(F,nrow = nrow(tmp), ncol = ncol(tmp))
selection[c(i1,i2,i3),] = T

par(mar=c(1,1,1,1))
p2 = heatmap.2(tmp, scale="row", col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),Colv = NA,labCol = "",
trace = "none", add.expr = makeRects(selection), dendrogram = "none" ,
          ColSideColors = my_col,symm = F, labRow = "",key = T,
         #lmat=rbind(c(4,1),c(3,2)),
         #lhei=c(0.05,0.85),
         #lwid = c(0.15,1),
         )



