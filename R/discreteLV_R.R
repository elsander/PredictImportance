#' Simulate a discrete-time Lotka-Volterra time series
#'
#' Simulate discrete-time Lotka-Volterra population dynamics
#' with perturbed growth rates
#'
#' @param rmat a matrix with S (number of species) rows and simtime
#' cols containing the growth rates at each time step (potentially
#' perturbed at each time step) 
#' @param alphas an S by S matrix of per capita interaction strengths
#' @param n0s a vector of length S containing species abundances
#' @param deltat the size of the time step
#' @param simtime the total number of time steps to simulate
#'
#' @return an S by simtime matrix of species abundances at each timestep
#'
#' @export

discreteLV <- function(rmat, alphas, n0s, deltat, simtime){
    S <- dim(alphas)[1]
    ## we only want to keep each integer time point
    ns <- matrix(0, S, simtime)
    ## we keep intermediate time step abundances here
    ntmp <- n0s
    for(i in 1:simtime){
        ns[,i] <- ntmp
        for(j in seq(i, i+1, by = deltat)){
            ntmp <- ntmp*exp(deltat*(rmat[,i] + alphas %*% ntmp))
        }
    }
    ns[,ncol(i)] <- ntmp
    return(ns)
}
