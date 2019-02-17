
#' preprocess data, remove non expressed genes
#'
#' @param data normalized transcripts
#' @return preprocessed transcripts


prep = function(data){
    data = round(data)
    dm = apply(data, 1, max)
    nz = which(dm > 0)
    return(data[nz,])
}
