load("spla1.RData")

###in total 7 subtypes
###proportion of subtypes: p and q in two conditions

##index of Equivalent distributed gene: ED
##index of Differential distributed gene: DD
##rn: gene name

##EDDM DD genes identified by MAST, probability of a gene being :
##EDD_sc DD genes idenfitied by scDD, probability of a gene being : p_sc 
##De_dd DE genes identified by DESeq2, probability of a gene being : fcHurdle model has column fdr
##EDDb DD genes identified by scDDboost, probability of a gene being DD: pdd7


##ROC plot
##scDDboost ROC
TPR = rep(0, 1001)
FPR = rep(0, 1001)

for(i in 1:1001){
  tmp = which(pdd7 > (i-1)/1000)
  TPR[i] = length(intersect(tmp, DD))/length(DD)
  FPR[i] = (length(tmp) - length(intersect(tmp, DD)))/nED
}
TPR_1 = TPR
FPR_1 = FPR



##scDD_roc
for(i in 1:1001){
  tmp = which(p_sc < (i-1)/1000)
  TPR[i] = length(intersect(tmp, DD))/length(DD)
  FPR[i] = (length(tmp) - length(intersect(tmp, DD)))/nED
}

TPR_sc = TPR
FPR_sc = FPR

##MAST roc
for(i in 1:1001){
  fcHurdleSig <- merge(fcHurdle[fdr<(i-1)/1000], as.data.table(mcols(MNZ10)), by='primerid')
  tmp = as.numeric(fcHurdleSig$gene)
  TPR[i] = length(intersect(tmp, DD))/length(DD)
  FPR[i] = (length(tmp) - length(intersect(tmp, DD)))/nED
}

TPR_m = TPR
FPR_m = FPR



##DEseq2 
for(i in 1:1001){
  tmp = which(res$pvalue < (i-1)/1000)
  TPR[i] = length(intersect(tmp, DD))/length(DD)
  FPR[i] = (length(tmp) - length(intersect(tmp, DD)))/nED
}

TPR_DE = TPR
FPR_DE = FPR


df = data.frame(x = c(FPR_1,FPR_sc,FPR_m,FPR_DE),
                y = c(TPR_1,TPR_sc,TPR_m,TPR_DE),
                method = rep(c("scDDboost","scDD","MAST","DESeq2"),each = 1001) )
clr = c("red","green","blue","pink")
ltp = factor(1:4)
pp = ggplot(df, aes(x,y,group = method)) + geom_line(aes(colour = factor(method),linetype = factor(method)), size = 2)
pdf("roc_sim.pdf")
pp + geom_abline(intercept = 0, slope = 1, size = 2) + theme(
  axis.text.x = element_text(face="bold", color="#993333", 
                             size=14),
  axis.text.y = element_text(face="bold", color="#993333", 
                             size=14),
  panel.background = element_rect(
    fill = 'white', colour = 'black'),
  panel.grid.minor.x = element_line(size = 0.5),
  panel.grid.minor.y = element_line(size = 0.5),
  panel.grid.major.x = element_line(size = 0.5),
  panel.grid.major.y = element_line(size = 0.5),
  panel.grid.major = element_line(colour = "grey")) + guides(fill = F)+ ylab("true positive rate") + xlab("false positive rate")
dev.off()

plot(FPR_1, TPR_1, type = "l", lwd = 8, col = "magenta", xlab = 'n', ylab = 'n')
lines(FPR_sc, TPR_sc, lwd = 8, col = "deepskyblue")
lines(FPR_m, TPR_m, lwd = 8, col = "khaki1")
lines(FPR_DE,TPR_DE, lwd = 8, col = "lawngreen")
abline(0,1, lwd = 8)






###pca plot
pc = prcomp(t(log(data_count + 1)))
library(ggfortify)
pdf("pc.pdf")
autoplot(pc, colour = as.factor(gp_label))+
theme(panel.background = element_rect(
fill = 'white', colour = 'black'),
panel.grid.major = element_line(colour = "grey"),
axis.text.x = element_blank(),
axis.text.y = element_blank())+ guides(fill = F)
dev.off()







