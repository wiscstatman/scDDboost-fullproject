library(ggplot2)
library(ggpubr)
library(scDDboost)

### need to remove the # of first line in those txt file,(result from runD3E

PATH = "./D3E_res/"
filenames = list.files(PATH)
d3eRes = list()
i = 1
for(FILE in filenames)
{
    d3eRes[[i]] = read.table(paste0(PATH,FILE), sep = "\t", header = T)
    i = i + 1
    
}

load("./DEC_EC.RData")
genenames = rownames(data_counts)

thre = 0.05

EDDb = lsz(newPDD,thre)

EDD_des = which(res_des$padj <  thre)
EDDM = which(res_mast[,"hurdle","Pr(>Chisq)"] < thre)
EDD_sc = which(res_scdd$combined.pvalue.adj < thre)

for(I in 1:length(d3eRes))
{
    IND = genenames %in% d3eRes[[I]]$GeneID
    
    #IND = d3eRes[[I]]$GeneID %in% GeneB
    tmp = which(IND)
    d3eRes[[I]]$pB = newPDD[tmp]
    d3eRes[[I]]$pD = res_des$padj[tmp]
    d3eRes[[I]]$pM = p.adjust(res_mast[tmp,"hurdle","Pr(>Chisq)"],"fdr")
    d3eRes[[I]]$pS = res_scdd$combined.pvalue.adj[tmp]
    #tmp_control = which(!IND)
#     d3eB[[I]] = d3eRes[[I]][tmp,]
#     d3eB_C[[I]] = d3eRes[[I]][tmp_control,]
}
# d3eB = do.call(rbind,d3eB)
# d3eB_C = do.call(rbind,d3eB_C)
d3e = do.call(rbind,d3eRes)


## boxplot for comparison
WHICH1 = which(d3e$pB > 0.95)

others = union(which(d3e$pD < 0.01),union(which(d3e$pS < 0.01), which(d3e$pM < 0.01)))

uni = setdiff(WHICH1,others)
#uni = order(-d3e$pB)[1:500]
#uni = WHICH1
#common = intersect(WHICH1,rest)


a1 = d3e$a1
a2 = d3e$a2
b1 = d3e$b1
b2 = d3e$b2
g1 = d3e$g1
g2 = d3e$g2

Type = rep(c(rep("others",length(others)), rep("uni", length(uni))), 3)

variety = c(rep("mean of activation rate", length(others) + length(uni)),
            rep("var of activation rate", length(others) + length(uni)),
           rep("rate of transcription", length(others) + length(uni)))

variety = c(rep("promoter activation", length(others) + length(uni)),
            rep("promoter inactivation", length(others) + length(uni)),
           rep("rate of transcription", length(others) + length(uni)))

LFC = c()

#tmp1 = a1 / (a1 + b1)
#tmp2 = a2 / (a2 + b2)
tmp1 = a1
tmp2 = a2
tmp = abs(log(tmp1 / tmp2))

tmp = tmp

LFC = c(LFC, c(tmp[others],tmp[uni]))


#tmp1 = a1 * b1 / (a1 + b1)^2 / (a1 + b1 + 1)
#tmp2 = a2 * b2 / (a2 + b2)^2 / (a2 + b2 + 1)
tmp1 = b1
tmp2 = b2
tmp = abs(log(tmp1 / tmp2))

tmp = tmp

LFC = c(LFC, c(tmp[others],tmp[uni]))

tmp1 = g1
tmp2 = g2
tmp = abs(log(tmp1 / (tmp2 + 1e-3) + 1e-3))

tmp = tmp

LFC = c(LFC, c(tmp[others],tmp[uni]))

data = data.frame(variety,Type,LFC)

new_order = with(data,reorder(variety,LFC,mean,na.rm = T))

par(mar=c(3,4,3,1))
#pdf("D3E_box.pdf")
myplot=boxplot(LFC ~Type * new_order , data=data  , boxwex=0.4 , ylab="change of bursting parameters",
        main="" , col=c("blue" , "red") ,  xaxt="n", outline = F)
 

my_names=sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names=my_names[seq(1 , length(my_names) , 2)]
axis(1, at = seq(1.5 , 7 , 2), labels = my_names , tick=FALSE , cex=0.3)
for(i in seq(0.5 , 20 , 2)){ abline(v=i,lty=1, col="grey")}
 
# Add a legend
legend("topleft", legend = c("others", "unique"), col=c("blue" , "red"),
       pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))
                
#dev.off()
