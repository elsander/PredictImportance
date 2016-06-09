#' Parameterize a Lotka-Volterra system given an adjacency matrix
#'
#' Parameterize an adjacency matrix with interaction strengths and
#' equilibrium abundances. Interaction strengths and abundances are drawn from a
#' lognormal distribution truncated 5 standard deviations above the
#' mean.
#'
#' @param Adj An adjacency matrix
#' @param Immigration A flag for whether we are parameterizing a system with
#' immigration or a closed system. Defaults to TRUE.
#'
#' @return a list with elements Mat (a matrix of interaction strengths) and
#'   Pop (a vector of equilibrium population sizes).
#'
#' @examples
#' mat <- BuildCascade()
#' LognormalParam(mat)
#'
#' @export

LognormalParam <- function(Adj, Immigration = TRUE){
    ##truncated lognormal distribution after 5 sd's from the mean
    ##the exp(.5) is the distribution mean
    cutoff <- sqrt(exp(1)*(exp(1)-1))*5 + exp(.5)

    ##now draw interaction strengths for prey
    S <- dim(Adj)[1]
    Adj <- Adj*rlnorm(S*S)
    while(sum((Adj > cutoff)*1) > 0){
        ##redraw values past the cutoff to make this a truncated lognormal
        Adj[Adj > cutoff] <- rlnorm(length(Adj[Adj > cutoff]))
    }

    ##prey are harmed by this interaction
    Adj <- Adj * -1
    ##now incorporate predator interactions
    ##predators don't benefit as much as prey are harmed
    Adj <- Adj + t(Adj)*rnorm(S*S, .2, .01)*-1

    ##get leading eigenvalue
    leadingev <- max(Re(eigen(Adj, only.values = TRUE)$values))

    if(Immigration == 0){
        ##If we aren't using immigration to stabilize, we need to stabilize
        ## with the diagonal entries
        
        ## if the system is unstable (positive leading eigenvalue),
        ## increase density-dependence on the diagonal to enforce
        ## stability.
        if(leadingev > 0){
            Adj <- Adj + diag(-1.1*leadingev, S, S)
        }
        Imm <- rep(0, S)
        ## draw equilibrium Ns
        Nstar <- rlnorm(S)
        while(sum((Nstar > cutoff)*1) > 0){
            ##redraw values past the cutoff to make this a truncated lognormal
            Nstar[Nstar > cutoff] <- rlnorm(length(Nstar[Nstar > cutoff]))
        }
    } else {
        ## Draw equilibrium Ns with immigration constraint
        ## cutoff for Nstars can be tightened based on max leading eigenvalue of Adj
        evalCutoff <- 1/leadingev
        if(evalCutoff < cutoff){
            cutoff <- evalCutoff
        }
        Nstar <- rlnorm(S)
        while(sum((Nstar > cutoff)*1) > 0){
            ##redraw values past the cutoff to make this a truncated lognormal
            Nstar[Nstar > cutoff] <- rlnorm(length(Nstar[Nstar > cutoff]))
        }

        ## Draw immigration terms
        Imm <- rep(0, S)
        for(i in 1:S){
            Imm[i] <- runif(1, leadingev*Nstar[i]^2, Nstar[i])
        }
    }
    return(list(Mat = Adj, Pop = Nstar, Imm = Imm))
}
