
## codes for heat map
load("heatPlot.RData")
## containing 4 things: 
## 1. cycleAll matrix for total cellcycle genes identified by scDDboost at 0.05 threshold
## 2. cycleUni matrix for cellcycle genes uniquely identified by scDDboost at 0.05 threshold
## 3. cd condition label: first 91 columns represent G1, last 76 columns represent G2M
## 4. my_group subtypes label, ordered from small to large of label values for each condition. 

library(RColorBrewer)
library(gplots)
##winsorize function
winsor = function (x, fraction=.05)
{
  if(length(fraction) != 1 || fraction < 0 ||
     fraction > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(fraction, 1-fraction))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  return(x)
}

## color for the color bar on top of the heat map
my_col=brewer.pal(9, "Set1")[my_group]

## apply winsorization and remove low expressed genes.
ApplyWin = function(x){
  n = nrow(x)
  for(i in 1:n){
    tmp = x[i,]
    x[i,] = winsor(tmp)
  }
  mean_win = apply(x,1,mean)
  low_expr = which(mean_win < 1)
  x = x[-low_expr,]
  return(x)
}

my_palette = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))

# heat map for total 
winAll = ApplyWin(cycleAll)

heatmap(winAll, scale="row", col = my_palette,Colv = NA,labCol = "", 
        add.expr = abline(v=91.5), Rowv = NA,ColSideColors = my_col,symm = F, labRow = "")

par(mar=c(1,1,1,1))
heatmap.2(winAll, scale="row", col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),Colv = NA,labCol = "",
trace = "none", add.expr = abline(v=91.5), dendrogram = "none" ,
ColSideColors = my_col,symm = F, labRow = "",key = F,
lmat=rbind(c(4,1),c(3,2)),
lhei=c(0.05,0.9),
lwid = c(0.15,1)
)

# heat map for unique
winUni = ApplyWin(cycleUni)

heatmap(winUni, scale="row", col = my_palette,Colv = NA,labCol = "", 
        add.expr = abline(v=91.5), Rowv = NA,ColSideColors = my_col,symm = F, labRow = "")


par(mar=c(1,1,1,1))
heatmap.2(winUni, scale="row", col = rev(colorRampPalette(brewer.pal(10, "RdBu"))(256)),Colv = NA,labCol = "",
trace = "none", add.expr = abline(v=91.5), dendrogram = "none" ,
ColSideColors = my_col,symm = F, labRow = "",key = F,
lmat=rbind(c(4,1),c(3,2)),
lhei=c(0.05,0.9),
lwid = c(0.15,1)
)
