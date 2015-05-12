#' Calculate network metrics for a food web
#'
#' Read in a matrix and calculate log closeness centrality,
#' log eigenvector centrality, and trophic level for each species.
#' Note that this function is used internally and is not likely to be
#' useful on its own. Use igraph functionality and getTL() for centrality
#' calculations beyond replicating simulations from Wootton et al.
#'
#' @param matfile path to parameterized nework
#'
#' @return a data frame containing log(closeness centrality),
#' log(eigenvector centrality), and trophic level for each
#' species in the network.
#'
#' @export

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
