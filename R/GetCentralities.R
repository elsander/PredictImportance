GetCentralities <- function(matfile){
    mat <- as.matrix(read.table(matfile, header = FALSE))
    mat <- (mat > 0) * 1 ##convert to adjacency matrix
    gmat <- graph.adjacency(mat)

    centralities <- data.frame(LogCloseness = log(closeness(gmat)),
                               LogEigenvector = log(abs(Re(eigen(mat2)$vectors[-1,1]))),
                               TrophicLevel = getTL(mat))
    return(centralities)
}
