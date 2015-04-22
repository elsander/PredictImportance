discreteLV <- function(rmat, alphas, n0s,  deltat, simtime){
    ## simulate discrete-time Lotka-Volterra time series given:
    ## rmat: a matrix with S (number of species) rows and simtime/deltat cols
    ##       containing the growth rates at each time step (potentially
    ##       perturbed at each time step)
    ## alphas: an S by S matrix of per capita interaction strengths
    ## n0s: a vector of length S containing species abundances
    ## deltat: the size of the time step
    ## simtime: the total number of time steps to simulate
    
    S <- dim(alphas)[1]
    ns <- matrix(0, S, simtime)
    ns[,1] <- n0s
    for(i in 1:(simtime-1)){
        ns[,i+1] <- ns[,i]*exp(deltat*(rmat[,i] + alphas %*% ns[,i]))
    }
    return(ns)
}
