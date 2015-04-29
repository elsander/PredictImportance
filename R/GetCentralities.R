GetCentralities <- function(matfile){
    mat <- as.matrix(read.table(matfile, header = FALSE))
    mat <- (mat > 0) * 1 ##convert to adjacency matrix
    gmat <- graph.adjacency(mat)

    ##add detritus node for eigenvector centrality
    mat2 <- rbind(rep(1,ncol(mat)), mat)
    rsum <- rowSums(mat2)
    rbool <- (rsum == 0)*1
    mat2 <- cbind(rbool, mat2)
    ##convert network to column stochastic matrix
    mat2 <- apply(mat2, 2, function(x) return(x/sum(x)))
    
    centralities <- data.frame(LogCloseness = log(closeness(gmat)),
                               LogEigenvector = log(abs(Re(eigen(mat2)$vectors[-1,1]))),
                               TrophicLevel = getTL(mat))
    return(centralities)
}
