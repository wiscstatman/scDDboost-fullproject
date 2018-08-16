#randomized with sc3

#using FUCCI data
load("sctrim.RData")
library(SingleCellExperiment)
library(SC3)
data_count = as.matrix(edat)
X <- SingleCellExperiment(assays = list(normcounts = data_count), 
                          colData = colnames(edat))
counts(X) <- normcounts(X)
logcounts(X) <- log2(normcounts(X) + 1)
##prepare for further sc3 reladted cal
rowData(X)$feature_symbol = rownames(edat)
X <- sc3_prepare(X)
dst.sc3 = sc3_calc_dists(X)
tran.D = sc3_calc_transfs(dst.sc3)
pre_output = sc3_kmeans(tran.D,3)
consen.D = sc3_calc_consens(pre_output)
print(names(metadata(consen.D)$sc3))
consen_matrix = metadata(consen.D)$sc3$consensus[[1]]$consensus
dim(consen_matrix)

