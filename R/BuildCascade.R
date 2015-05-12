#' Build a Cascade Model Network
#'
#' Generates a random network of a given size and approximate
#' connectance according to the cascade model
#'
#' @param S Number of species in network; defaults to 50
#' @param C Network connectance, defined as \eqn{2L/(S(S-1))},
#'   where L is the number of links. Defaults to .1. The resulting
#'   matrix will have a connectance of C +/- \eqn{2sqrt(.5C(1-C)S(S-1))}.
#'
#' @return A weakly connected, unsigned adjacency matrix of size (S,S).
#'
#' @examples
#' BuildCascade()
#'
#' @references Cohen, J.E. and C.M. Newman. 1985. A Stochastic Theory of
#'   Community Food Webs: I. Models and Aggregated Data. Proc. R. Soc.
#'   Lond. B. 224: 421-448; DOI: 10.1098/rspb.1985.0042.
#'
#' @export

########################################
## NOTE: CONNECTANCE IS DEFINED AS C = L / (S * (S-1) / 2)
## WHICH IS DIFFERENT FROM THE USUAL DEFINITION, BUT 
## MAKES MORE SENSE AS WE ARE MODELING PAIRS OF INTERACTIONS
########################################
BuildCascade <- function(S = 50, C = .1){
    Success <- FALSE
    DesiredL <- C * S * (S-1) / 2
    CutOff <- 2 * sqrt((1 - C) * C * S * (S-1) / 2)

    M <- matrix(0, S, S)
    while(Success == FALSE){
        M[upper.tri(M, diag = FALSE)] <- runif(S*(S-1)/2)
        M <- (M >= (1-C)) * 1

        ##check that:
        ## (1) connectance is about right
        ## (2) graph is weakly connected
        ## if these are met, return the network
        if (abs(sum(M) - DesiredL) < CutOff){
            g <- graph.adjacency(M)
            if (is.connected(g, "weak")){
                Success <- TRUE
            }
        } 

    }
    return(M)
}
