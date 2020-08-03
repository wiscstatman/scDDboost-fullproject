### fucci data analysis 
library(scDDboost)
load("Fucci.RData")

### content in the data ########

## data_counts: the Fucci data set G1 vs G2/M
## cd: condition label, 1 for G1, 2 for G2/M
## D_c: distance matrix by default clustering
## genenames: gene names for Fucci data
## newPDD: estimated PDD by scDDboost
## res_des: results by DESeq2
## res_mast: results by MAST
## res_scdd: results by scDD
## GO_pool: gene ontology annotation for cell cycle related genes (GO:0007049)
## WHICH: index for which genes against the constant shape assumption


#############

### get PDD by scDDboost
## get distance matrix, using 10 cores
D_c = cal_D(data_counts,10)

## using 4 cores
pdd = PDD(data_counts, cd, 10, D_c)

## should output the same results as newPDD


### get adjusted p values of other DE methods

p_des = res_des$padj

p_mast = res_mast[,3,3]

p_scdd = res_scdd$nonzero.pvalue.adj


## t-test
tpval = rep(0,nrow(data_counts))
for(i in 1:nrow(data_counts)){
  x = data_counts[i,which(cd == 1)]
  y = data_counts[i,which(cd == 2)]
  fit = t.test(x,y)
  tpval[i] = fit$p.value
}

t_pval = p.adjust(tpval,"fdr")


### get uniquely identified genes against constant shape 
t_ = which(t_pval < 0.05) 

des_ = which(p_des < 0.05)
mast_ = which(p_mast < 0.05)
scdd_ = which(p_scdd < 0.05)

UNION = union(scdd_,union(mast_,union(des_,t_)))

trouble = intersect(which(newPDD > 0.95), WHICH)

## have 468 genes called DD but against constant shape
length(trouble)

## 323 of them also called by other methods
length(intersect(UNION,trouble))

uniq = setdiff(trouble,UNION)

## most of them are not cell cycle related genes
tolower(genenames[uniq]) %in% sapply(GO_pool$Symbol,tolower)

##### let us look at some uniquely called genes but against constant shape

longVec = c()
L = 5
sam = sample(uniq,5)
for(i in 1:L){
  longVec = c(longVec,data_counts[sam[i],which(cd == 1)],data_counts[sam[i],which(cd == 2)])
}

dfm = data.frame(value = longVec)
dfm$condition = as.factor(rep(cd,L))
dfm$gene = rep(1:L,each = length(cd))

ggplot(dfm, aes(x=condition,y = value,color = condition)) + 
  geom_violin() + facet_wrap(.~gene,nrow = 2,scales = "free") + theme_classic() + labs(y = "log expression", x = " ") 





