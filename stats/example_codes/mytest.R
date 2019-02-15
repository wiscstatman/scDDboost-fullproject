
# a new example

#source("https://bioconductor.org/biocLite.R")
#biocLite("Oscope")
#biocLite("BiocParallel")

#install("/Users/michaelnewton/wip/scDDboost/scDDboost/scDDboost_pkg")

library(scDDboost)

###first to calculate 
D_c = cal_D(data_count, 10) ##calculate distance matrix using 10 cores
hp = c(1, rep(1,nrow(data_count))) ## set default hyper parameter of EBSeq to be 1
gcl = 1:nrow(data_count) ##gene cluster are set to be each gene form one cluster

K = 5 ##number of subtypes 
Posp = pat(K)[[1]] ##possible partitions
sz = MedianNorm(data_count) ##sizefactor required by EBSeq


#calculating Posterior probability of differential distributed. 
###INPUT:
##data_count: transcript matrix
##cd: conditon label
##ncores: number of cores 
##K: number of subtypes
##D_c: distance matrix
##sz: sizefactor of cells in the EBSeq part
##hp: hyper parameter in EBSeq
##iter: max number of iteration for EM, typically converge in 5 times
##random: bool indicating whether using random distance
##lambda: parameter for exponential noise
##nrandom: number of random times

pdd5 = PDD(data_count,cd, ncores = 10,K, D_c, sz, hp, pat(K)[[1]],iter = 10, random = T, lambda = 1, nrandom = 100)
EDDb = which(pdd5 > 0.95)