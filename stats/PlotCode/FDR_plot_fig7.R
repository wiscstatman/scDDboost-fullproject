library(ggplot2)


loadName = c("GSE79102.RData", "GSE48968.RData", "TASIC.RData", "GSE94383.RData",
             "DEC_EC.RData", "Fucci.RData", "GSE60749.RData","GSE74596.RData",
             "EMTAB2805.RData","GSE63818.RData")

nSet = length(loadName)

##in total we have 10 datasets
##each one has 5 replicates

mst_fdr = rep(0,nSet * 5)
des_fdr = mst_fdr
scdd_fdr = mst_fdr
scddboost_fdr = mst_fdr
scddboost_sc3_fdr = mst_fdr



for(II in 1:nSet){
  for(J in 1:5){
    if(J > 1)
      TMP = paste0("FDR", (J -1))
    else
      TMP = "FDR"
    TMP = paste0(TMP,"/res_")
    x = loadName[II]
    load(paste0(TMP,x))
    
    mst = res_mast[,"hurdle","Pr(>Chisq)"]
    print(x)
    nr = nrow(data_null)
    print("MAST")
    print(length(which(mst < 0.05))/ nr)
    mst_fdr[nSet * J - nSet + II] = length(which(mst < 0.05))/ nr
    print("DESeq2")
    print(length(which(res_des$padj < 0.05))/nr)
    des_fdr[nSet * J - nSet + II] = length(which(res_des$padj < 0.05))/nr
    print("scDD")
    print(length(which(res_scdd$combined.pvalue.adj < 0.05))/nr)
    scdd_fdr[nSet * J - nSet + II] = length(which(res_scdd$combined.pvalue.adj < 0.05))/nr
    print("scDDboost")
    print(length(which(res_scddboost$RPDD > 0.95))/nr)
    scddboost_fdr[nSet * J - nSet + II] = length(which(res_scddboost$RPDD > 0.95))/nr
  }
}


##FDR plot
df_fdr = data.frame(methods = rep(c("DESeq2", "MAST", "scDD", "scDDboost"),each = length(mst_fdr)),
                    y = c(des_fdr, mst_fdr, scdd_fdr, scddboost_fdr))

p = ggplot(df_fdr,aes(methods,y, color = methods)) + geom_boxplot(alpha = 0.7,
                                                                  outlier.colour = "#1F3552", outlier.shape = 20) + geom_jitter()
pdf("fdr.pdf", height = 5, width = 8)
p + theme(panel.background = element_rect(
  fill = 'white', colour = 'black'),
  axis.text.x = element_text( color="black", 
                              size= 14),
  axis.text.y = element_text(face="bold", color="#993333", 
                             size=14),
  legend.text=element_text(size = 14),
  legend.title = element_blank(),
  axis.title.x=element_blank(),
  axis.title=element_text(size=14,face="bold"),
  panel.grid.minor.x = element_line(size = 0.5),
  panel.grid.minor.y = element_line(size = 0.5),
  panel.grid.major.x = element_line(size = 0.5),
  panel.grid.major.y = element_line(size = 0.5),
  panel.grid.major = element_line(colour = "grey"))+ ylab("False Positive Rate") + geom_hline(yintercept=0.05, linetype="dashed", color = "red")

dev.off()
