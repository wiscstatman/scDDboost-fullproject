
load("results.RData")


names(pDD8) <- rn ## gene symbols on posterior probs
names(P_scDD) <- rn
names(P_mast) <- rn

yboost <- 1*(pDD8 >= 0.99)
yscDD <- 1*(P_scDD <= 0.01 )
ymast <- 1*(P_mast <= 0.01 )

nc <- 5
library(allez)
library(org.Hs.eg.db)

se <- as.list(org.Hs.egSYMBOL2EG)

gs <- names(yboost)
ee <- se[gs]
ok.name <- sapply(ee, function(x) !is.null(x))

yboost2 <- yboost[ok.name]
yscDD2 <- yscDD[ok.name]
ymast2 <- ymast[ok.name]
ff <- ee[ok.name]   ## at least one Entrez ID


## pick one
eg1 <- sapply(ff, function(x) x[1] )
names(yboost2) <- eg1

ag <- allez(yboost2, lib="org.Hs.eg" ) 

u <- ag$aux$set.data
v1 <- (1:nrow(u))[u$go_id=="GO:0007049"]   ## cell cycle
##v1 <- (1:nrow(u))[u$go_id=="GO:0051726"]   ## regulation of cell cycle

## gene symbols in two  GO terms
w1 <- u[v1,3] # gene symbols  

### let's march down prob lists and count GO genes, instead of thresholding


pDD82 <- pDD8[ok.name]
P_scDD2 <- P_scDD[ok.name]
P_mast2 <- P_mast[ok.name]

rDD83 <- rank( -pDD82 )
y1 <- rep(0, length(rDD83)); names(y1) <- names(rDD83)
y1[w1] <- 1
plot( cumsum(y1) , pch="." )
abline( 0, sum(y1)/length(y1), col="red" )

r_scDD3 <- rank( P_scDD2 )
y2 <- rep(0, length(r_scDD3)); names(y2) <- names(r_scDD3)
y2[w1] <- 1
plot( cumsum(y2) , pch="." )
abline( 0, sum(y2)/length(y2), col="red" )



r_mast3 <- rank(P_mast2)



## no sign that any of these DE schemes are ranking cell cycle better
## than other genes....inconclusive

