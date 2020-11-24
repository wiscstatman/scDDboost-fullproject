load("simplot.RData")
library(ggplot2)
library(ggpubr)

K = 6
DF.K = data.frame(X = pDD[[K]],Y = pDD[[K + 1]])


basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(colour = NA, aes(
  fill = cut(..count.., c(0, 30, 70,
                          100, Inf))),
  size = 3, bins = 100, alpha = 0.75,
) + coord_equal() +
  labs(fill = NULL) 


# 
# basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(aes(x = X, y = Y,fill = NULL,
#                              colour = cut(..count.., c(0, 10, 50, 100,
#                                                       500,  Inf))),
#                          geom = "point", size = 3, bins = 30, alpha = 0.75,
#                          data = DF.K) 


p1 = basic.hexmap + scale_colour_brewer(palette = "OrRd",
                                        labels = c("<50", "50-100",
                                                   "100-500", "500-1000",
                                                   ">1000")) + theme(
                                                     axis.text.x = element_text(face="bold", color="#993333", 
                                                                                size=6),
                                                     axis.text.y = element_text(face="bold", color="#993333", 
                                                                                size=6),
                                                     panel.background = element_rect(
                                                       fill = 'white', colour = 'black'),
                                                     legend.title=element_blank(),
                                                     legend.position = c(0.85, 0.19),
                                                     legend.text=element_text(size=6),
                                                     panel.grid.minor.x = element_line(size = 0.4),
                                                     panel.grid.minor.y = element_line(size = 0.4),
                                                     panel.grid.major.x = element_line(size = 0.4),
                                                     panel.grid.major.y = element_line(size = 0.4),
                                                     panel.grid.major = element_line(colour = "grey"),
                                                     axis.title.x = element_blank(),
                                                     axis.title.y = element_blank())

#axis.title.y = element_text(face = "bold", size = rel(1.2), angle = 90),
#axis.title.x = element_text(face = "bold", size = rel(1.2), angle = 00),
#) + ylab(paste0("PDD at ",K + 1)) + xlab(paste0("PDD at ",K))


K = 7
DF.K = data.frame(X = pDD[[K]],Y = pDD[[K + 1]])


basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(colour = NA, aes(
  fill = cut(..count.., c(0, 30, 70,
                          100, Inf))),
  size = 3, bins = 100, alpha = 0.75,
) + coord_equal() +
  labs(fill = NULL) 


# 
# basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(aes(x = X, y = Y,fill = NULL,
#                              colour = cut(..count.., c(0, 10, 50, 100,
#                                                       500,  Inf))),
#                          geom = "point", size = 3, bins = 30, alpha = 0.75,
#                          data = DF.K) 


p2 = basic.hexmap + scale_colour_brewer(palette = "OrRd",
                                        labels = c("<50", "50-100",
                                                   "100-500", "500-1000",
                                                   ">1000")) + theme(
                                                     axis.text.x = element_text(face="bold", color="#993333", 
                                                                                size=6),
                                                     axis.text.y = element_text(face="bold", color="#993333", 
                                                                                size=6),
                                                     panel.background = element_rect(
                                                       fill = 'white', colour = 'black'),
                                                     legend.title=element_blank(),
                                                     legend.position = c(0.85, 0.19),
                                                     legend.text=element_text(size=6),
                                                     panel.grid.minor.x = element_line(size = 0.4),
                                                     panel.grid.minor.y = element_line(size = 0.4),
                                                     panel.grid.major.x = element_line(size = 0.4),
                                                     panel.grid.major.y = element_line(size = 0.4),
                                                     panel.grid.major = element_line(colour = "grey"),
                                                     axis.title.x = element_blank(),
                                                     axis.title.y = element_blank())

#axis.title.y = element_text(face = "bold", size = rel(1.2), angle = 90),
#axis.title.x = element_text(face = "bold", size = rel(1.2), angle = 00),
#) + ylab(paste0("PDD at ",K + 1)) + xlab(paste0("PDD at ",K))


pdf("sim.pdf")
ggarrange(p1, p2, ncol = 2,common.legend = TRUE, legend="right")
dev.off()



