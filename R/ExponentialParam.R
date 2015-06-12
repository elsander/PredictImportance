ExponentialParam <- function(Adj){
    ##draw interaction strengths for prey
    ##have mean = exp(.5) to be similar to lognormal distribution
    ##it will be a little high because the lognormal is truncated
    S <- dim(Adj)[1]
    Adj <- Adj*rexp(S*S, rate = exp(.5))

    ##prey are harmed by this interaction
    Adj <- Adj * -1
    ##now incorporate predator interactions
    ##predators don't benefit as much as prey are harmed
    Adj <- Adj + t(Adj)*rnorm(S*S, .2, .01)*-1

    Nstar <- rlnorm(S)

    ##get leading eigenvalue
    leadingev <- max(Re(eigen(Adj, only.values = TRUE)$values))
    ## if the system is unstable (positive leading eigenvalue),
    ## increase density-dependence on the diagonal to enforce
    ## stability.
    if(leadingev > 0){
        Adj <- Adj + diag(-1.1*leadingev, S, S)
    }

    return(list(Mat = Adj, Pop = Nstar))
}
