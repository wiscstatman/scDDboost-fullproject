
#' generating new distance matrix between cells after EBSeq analysis
#' 
#' @param X cell clusters
#' @param part partition patterns, obtained through Posp function
#' @param Y index of map pattern
#' @return distance matrix of cells





new_D<-function(X,part,Y){
    .Call('new_D',X = X,part = part,Y = Y)
}
