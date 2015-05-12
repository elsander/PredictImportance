#' Build a Niche Model Network
#'
#' Generates a random network of a given size and approximate
#' connectance according to the niche model
#'
#' @param S Number of species in network; defaults to 50
#' @param C Network connectance, defined as 2L/(S(S-1)),
#'   where L is the number of links. Defaults to .1. The resulting
#'   matrix will have a connectance of C +/- 2sqrt(.5C(1-C)S(S-1)).
#'
#' @return M A weakly connected, unsigned adjacency matrix of size (S,S).
#'
#' @examples
#' BuildNiche()
#'
#' @references Williams, R.J., and N.D. Martinez. 2000.
#'   Simple rules yield complex food webs. Nature 404(6774): 180-3.
#'
#' @export

########################################
## IMPORTANT: CONNECTANCE IS DEFINED AS C = L / (S * (S-1) / 2)
## WHICH IS DIFFERENT FROM THE USUAL DEFINITION, BUT 
## MAKES MORE SENSE AS WE ARE MODELING PAIRS OF INTERACTIONS
########################################
BuildNiche <- function(S = 50, C = .1){
    ## Build a connected network with the right connectance using the
    ## niche model
    ## Niche model Ref: Richard J. Williams & Neo D. Martinez (2000)
    ##"Simple rules yield complex food webs", Nature
    Success <- FALSE
    DesiredL <- C * S * (S-1) / 2
    CutOff <- 2 * sqrt((1 - C) * C * S * (S-1) / 2)
    while (Success == FALSE){
        M <- matrix(0, S, S)
        ## assign the niche value for all species
        ## Note: the niche values of the species are ordered
        ni <- sort(runif(S, 0, 1))
        ## determine the radius of the food spectrum for each species 
        ## the radius is drawn randomly from a ni * Beta(1, beta) distribution
        beta <- 1.0 / C - 1 
        ri <- rbeta(S, 1, beta) * ni
        ## Set species with the lowest niche value to be the basal species
        ri[1] <- 0
        ## The center of the food spectrum for each species are 
        ## drawn uniformly from an uniform distribution on [ri/2, ni]
        ci <- numeric(S)
        for(i in 1:S){
            ci[i] <- runif(1, ri[i] / 2, min(ni[i], 1 - ri[i] / 2))
            ## Correction as in Allesina et al Science 2008
            ## Determine the boundary for the food spectrum for each species:
            upper <- ci[i] + ri[i] / 2
            lower <- ci[i] - ri[i] / 2
            ## check which species are falling into that interval: [lower, upper]
            ## and set the corresponding value in the ajacency matrix M to be 1,
            ## indicating a link is established between the two species.
            for(j in 1:S){
                if(ni[j] > lower & ni[j] < upper)
                    M[j, i] <- 1
            }
        }
        diag(M) <- 0

        ##check that:
        ## (1) connectance is about right
        ## (2) graph is weakly connected
        ## (3) there are no species that both eat each other
        ## if these are met, return the network
        if (abs(sum(M) - DesiredL) < CutOff){
            g <- graph.adjacency(M)
            if (is.connected(g, "weak")){
                M2 <- M + t(M)
                ##if there are any M2 == 2, then there is
                ##at least one pair of spp that both eat
                ##each other, which is problematic for
                ##parameterization
                if(sum(M2 == 2) == 0){
                    Success <- TRUE
                }
            }
        } 
    }
    return(M)
}
