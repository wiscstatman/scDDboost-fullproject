##code for Fig10 bursting param

load("empirical_data/TSC_D3E.RData")

tmp = abs(log(d3e$a1 + 1) - log(d3e$a2 + 1))

DF.K = data.frame(X = d3e$pB,Y = tmp / max(tmp))




basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(colour = NA, aes(
  fill = cut(..count.., c(0, 1, 5,
                          10, Inf))),
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


tmp = abs(log(d3e$b1 + 1) - log(d3e$b2 + 1))


DF.K = data.frame(X = d3e$pB,Y = tmp / max(tmp))


basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(colour = NA, aes(
  fill = cut(..count.., c(0, 1, 5,
                          10, Inf))),
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

tmp = abs(log(d3e$g1 + 1) - log(d3e$g2 + 1))

DF.K = data.frame(X = d3e$pB,Y = tmp / max(tmp))


basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(colour = NA, aes(
  fill = cut(..count.., c(0, 1, 5,
                          10, Inf))),
  size = 3, bins = 100, alpha = 0.75,
) + coord_equal() +
  labs(fill = NULL) 


# 
# basic.hexmap <- ggplot(DF.K,aes(X,Y)) + stat_bin_hex(aes(x = X, y = Y,fill = NULL,
#                              colour = cut(..count.., c(0, 10, 50, 100,
#                                                       500,  Inf))),
#                          geom = "point", size = 3, bins = 30, alpha = 0.75,
#                          data = DF.K) 


p3 = basic.hexmap + scale_colour_brewer(palette = "OrRd",
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


pdf("d3eplot.pdf")
ggarrange(p1, p2,p3, ncol = 3,common.legend = TRUE, legend = "right")
dev.off()