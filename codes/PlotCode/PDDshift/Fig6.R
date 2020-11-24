library(ggplot2)
library(ggpubr)
library(scDDboost)


## boxplot for comparison






variety = c(rep("GSE79102NULL", length(pdditer) + length(pdd0)),
rep("NULL", length(pdditer) + length(pdd0)),
rep("1NULL", length(pdditer) + length(pdd0)))

Type = c()

variety = c()

pdd = c()

load("Desktop/scDDboost/stats/PlotCode/PDDshift/EMTAB2805NULL.RData")

pdd = c(pdd, c(pdditer,pdd0))

Type = c(Type,c(rep("pooled",length(pdditer)), rep("solely", length(pdd0))) )

variety = c(variety, rep("EMTAB2805", length(pdditer) + length(pdd0)))

load("Desktop/scDDboost/stats/PlotCode/PDDshift/TASICNULL.RData")

pdd = c(pdd, c(pdditer,pdd0))

Type = c(Type,c(rep("pooled",length(pdditer)), rep("solely", length(pdd0))) )

variety = c(variety, rep("TASIC", length(pdditer) + length(pdd0)))

load("Desktop/scDDboost/stats/PlotCode/PDDshift/FUCCINULL.RData")

pdd = c(pdd, c(pdditer,pdd0))

Type = c(Type,c(rep("pooled",length(pdditer)), rep("solely", length(pdd0))) )

variety = c(variety, rep("FUCCI", length(pdditer) + length(pdd0)))

data = data.frame(variety,Type,pdd)

new_order = with(data,reorder(variety,pdd,mean,na.rm = T))

par(mar=c(3,4,3,1))
pdf("nullpdd_box.pdf")
myplot=boxplot(pdd ~Type * new_order , data=data  , boxwex=0.4 , ylab="PDD under null case",
main="" , col=c("slateblue1" , "tomato") ,  xaxt="n", outline = F)


my_names=sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names=my_names[seq(1 , length(my_names) , 2)]
axis(1, at = seq(1.5 , 7 , 2), labels = my_names , tick=FALSE , cex=0.3)
for(i in seq(0.5 , 20 , 2)){ abline(v=i,lty=1, col="grey")}

# Add a legend
legend("topleft", legend = c("pooled", "unpooled"), col=c("slateblue1" , "tomato"),
pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))

dev.off()
