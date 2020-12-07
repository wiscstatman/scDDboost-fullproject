## scDDboost with sc3 clustering


library(scDDboost)
library(SC3)
library(SingleCellExperiment)



load("######### simulation data ########")
### simulation data is put in
## /z/Comp/newtongroup/xiuyu/reviewSimu/
## due to the large size



data_counts = as.matrix(data_counts)
X <- SingleCellExperiment(assays = list(normcounts = data_counts),
                          colData = colnames(data_counts))
counts(X) <- normcounts(X)
logcounts(X) <- log2(normcounts(X) + 1)
##prepare for further sc3 reladted cal
rowData(X)$feature_symbol = rownames(data_counts)
X <- sc3_prepare(X)
sce <- sc3_estimate_k(X)
K = metadata(sce)$sc3$k_estimation
dst.sc3 = sc3_calc_dists(X)
tran.D = sc3_calc_transfs(dst.sc3)
pre_output = sc3_kmeans(tran.D,K)
consen.D = sc3_calc_consens(pre_output)
print(names(metadata(consen.D)$sc3))
consen_matrix = metadata(consen.D)$sc3$consensus[[1]]$consensus
dim(consen_matrix)


pdd_sc3 = PDD(data_counts,cd, ncores = 10, D = consen_matrix)

## FDR
1 - length(intersect(lsz(pdd_sc3,0.05),DD)) / length(lsz(pdd_sc3,0.05))
