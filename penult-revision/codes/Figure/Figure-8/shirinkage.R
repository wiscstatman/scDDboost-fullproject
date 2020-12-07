Type = c()

variety = c()

clussz = c()



load("TASIC_clussz.RData")


##tmp is the size of the cluster each gene belongs to.
## others are index for those DD genes identified by other methods at 0.01 threshold
## uni are index for those uniquely identified DD genes


clussz = c(clussz, c(tmp[others],tmp[uni]))

Type = c(Type, c(rep("others",length(others)), rep("uni", length(uni))))

variety = c(variety,rep("TASIC", length(others) + length(uni)) )

load("EMTAB_clussz.RData")

clussz = c(clussz, c(tmp[others],tmp[uni]))


Type = c(Type, c(rep("others",length(others)), rep("uni", length(uni))))

variety = c(variety,rep("EMTAB2805", length(others) + length(uni)) )


load("FUCCI_clussz.RData")

clussz = c(clussz, c(tmp[others],tmp[uni]))

Type = c(Type, c(rep("others",length(others)), rep("uni", length(uni))))

variety = c(variety,rep("FUCCI", length(others) + length(uni)) )

data = data.frame(variety,Type,clussz)

new_order = with(data,reorder(variety,clussz,mean,na.rm = T))

par(mar=c(3,4,3,1))
pdf("clussz_box.pdf")
myplot=boxplot(clussz ~ Type * new_order , data=data  , boxwex=0.4 , ylab="cluster size on log scale",
main="right shifts of DE pattern distribution" , col=c("slateblue1" , "tomato") ,  xaxt="n", outline = F,
ylim = c(2,10))


my_names=sapply(strsplit(myplot$names , '\\.') , function(x) x[[2]] )
my_names=my_names[seq(1 , length(my_names) , 2)]
axis(1, at = seq(1.5 , 7 , 2), labels = my_names , tick=FALSE , cex=0.3)
for(i in seq(0.5 , 20 , 2)){ abline(v=i,lty=1, col="grey")}

# Add a legend
legend("topleft", legend = c("others", "unique"), col=c("slateblue1" , "tomato"),
pch = 15, bty = "n", pt.cex = 3, cex = 1.2,  horiz = F, inset = c(0.1, 0.1))

dev.off()
