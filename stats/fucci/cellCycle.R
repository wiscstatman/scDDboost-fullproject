load("results.RData")
library(scDDboost)

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

qv.t <- p.adjust( pv.t, method="BH" )  ## looks like a lot
qv.w <- p.adjust( pv.w, method="BH" )  ## l

#sum( qv.t <= 0.05 )  ## I find 3361...not sure how this compares to other methods
#sum( qv.w <= 0.05 )  ## I find 3496 ...not sure how this compares to other methods

genenames = rownames(data_counts)

des_dd = which(res_des$padj < 0.05)
scdd_dd = which(res_scdd$combined.pvalue.adj < 0.05)


mst = p.adjust(res_mast[,"hurdle","Pr(>Chisq)"],"fdr")
mast_dd = which(mst < 0.05)

scddb_dd = lsz(newPDD, 0.05)

t_dd = which(qv.t < 0.05)
w_dd = which(qv.w < 0.05)

## DD genes by other methods
others = union(des_dd,union(mast_dd,scdd_dd))
others = union(others,union(t_dd,w_dd))

## uniquely DD by scDDboost
uni = setdiff(scddb_dd, others)

## DD genes by all methods
all = union(scddb_dd,others)

## ED genes by all methods
ED = setdiff(1:nrow(data_counts),all)

## others identified but not identified by scDDboost
others_id = setdiff(others,scddb_dd)


## commonly identified DD
common = intersect(others,scddb_dd)

uni_gene = genenames[uni]

GO_pool = read.csv("/stats/fucci/GO_pool.csv")
smb = as.character(unique(GO_pool$Symbol))
## convert genenames to upper case for matching genenames from Fucci
smb = sapply(smb,function(x) toupper(x))

scb_genes = genenames[scddb_dd]
mst_genes = genenames[mast_dd]
sc_genes = genenames[scdd_dd]
des_genes = genenames[des_dd]
t_genes = genenames[t_dd]
w_genes = genenames[w_dd]

## uniquely identified cell cycle related genes by scDDboost
uni_cyc = intersect(smb,uni_gene)

## cell cycle genes by scDDboost
scb_cyc = intersect(scb_genes,smb)

## cell cycle genes by MAST
mst_cyc = intersect(mst_genes,smb)

## cell cycle genes by scDD
sc_cyc = intersect(sc_genes,smb)

## cell cycle genes by DESeq2
des_cyc = intersect(des_genes,smb)

## cell cycle genes by t test
t_cyc = intersect(t_genes,smb)

## cell cycle genes by wilcox test
w_cyc = intersect(w_genes,smb)



length(scb_cyc)
length(mst_cyc)
length(sc_cyc)
length(des_cyc)
length(t_cyc)
length(w_cyc)


