###roc curve plot

roc = function(est.P,DD,ED){
    TPR = rep(0, 1001)
    FPR = rep(0, 1001)
    for(i in 1:1001){
        tmp = which(est.P <= (i - 1) / 1000)
        TPR[i] = length(intersect(tmp, DD))/length(DD)
        FPR[i] = (length(tmp) - length(intersect(tmp, DD)))/nED
    }
    res = list()
    res$TPR = TPR
    res$FPR = FPR
    return(res)
}

#TPRvec = c(FPR_1,FPR_sc,FPR_m,FPR_DE)
#FPRvec = c(TPR_1,TPR_sc,TPR_m,TPR_DE)
#Namevec = c("scDDboost","scDD","MAST","DESeq2")

rocPlot = function(TPRvec, FPRvec, Namevec){
    if(length(TPRvec) != length(FPRvec)){
        print("number of TPR should equal number of FPR")
        return
    }
    
    df = data.frame(x = TPRvec,
    y = FPRvec,
    method = rep(Namevec,each = length(TPRvec[1])) )
    clr = c("red","green","blue","pink")
    ltp = factor(1:4)
    p = ggplot(df, aes(x,y,group = method)) + geom_line(aes(colour = factor(method),linetype = factor(method)), size = 2)
    pdf("roc_sim.pdf")
    p + geom_abline(intercept = 0, slope = 1, size = 2) + theme(
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
}

