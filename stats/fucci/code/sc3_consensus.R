#randomized with sc3

#using FUCCI data
load("sctrim.RData")
library(SingleCellExperiment)
library(SC3)
load("E2805.RData") #K = 6
load("G48968.RData") #K = 6
load("G45719.RData") # K = 6
load("G74596.RData") # K = 5
load("G63818.RData") # K = 5
load("Tasic_1.RData")# K = 7
load("G60749.RData") # K = 6
load("DEC_EC.RData") # K = 4
load("DEC-NPC.RData") # K = 7
load("H12.RData") # K = 0
data_count = as.matrix(edat)
X <- SingleCellExperiment(assays = list(normcounts = data_count), 
                          colData = colnames(data_count))
counts(X) <- normcounts(X)
logcounts(X) <- log2(normcounts(X) + 1)
##prepare for further sc3 reladted cal
rowData(X)$feature_symbol = rownames(data_count)
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
####

pdd_sc3 = PDD(data_count,cd, ncores = 10, K, consen_matrix, sz, hp, pat(K)[[1]],iter = 10, random = T, lambda = 10, nrandom = 10)
length(which(pdd_sc3 > 0.95))













