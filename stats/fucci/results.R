load("results.RData")
library(scDDboost)

## containing
## data_counts for the normalized counts of FUCCI
## cd condition label, 1 for G1 and 2 for G2
## newPDD as the PDD given by scDDboost

# would gave the number of DD genes under 0.05 threshold
listsize(newPDD,0.05)

# index for those DD genes under 0.05 threshold
scDDb_dd = lsz(newPDD,0.05)

## res_des is the result after applying DESeq2
# I used the pvalue adjusted by "BH"
des_dd = which(res_scdd$padj < 0.05)

## res_mast is the result for MAST

# apply BH adjusted to p value from MAST
mst = p.adjust(res_mast[,"hurdle","Pr(>Chisq)"],"fdr")

mast_dd = which(mst < 0.05)

## res_scdd is the result for scDD

# I used the combined p-value of scDD adjusted by "BH"
scdd_dd = which(res_scdd$combined.pvalue.adj < 0.05)
