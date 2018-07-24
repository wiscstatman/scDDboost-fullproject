load("results.RData")


##data_count is the transcripts for G1 and G2
##cd is the condition label
##rn genes names
##ccl refers to cell subgroup label under distance matrix D_c
##table(ccl[which(cd == 1)])
##table(ccl[which(cd == 2)]) ##different of subtypes proportions
##pDD8 refers to probability of DD given by scDDboost under 8 subtypes
##pDD7 refers to probibility of DD given by scDDboost under 7 subtypes
##plot(pDD7,pDD8)
##p_scDD : probability of DD given by scDD
##p_MAST: probability of DD given by MAST
##length(which(pDD8 > 0.95))
##length(which(p_scDD > 0.95))
##length(which(p_MAST > 0.95))
