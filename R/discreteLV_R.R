#' Simulate a discrete-time Lotka-Volterra time series
#'
#' Simulate discrete-time Lotka-Volterra population dynamics
#' with perturbed growth rates
#'
#' @param rmat a list of length S, containing splined forcing functions for
#' each species
#' @param alphas an S by S matrix of per capita interaction strengths
#' @param n0s a vector of length S containing species abundances
#' @param deltat the size of the time step
#' @param simtime the total number of time steps to simulate
#'
#' @return an S by simtime matrix of species abundances at each timestep
#'
#' @export

discreteLV <- function(rspline, alphas, n0s, deltat, simtime){
    S <- dim(alphas)[1]
    ## we only want to keep each integer time point
    ns <- matrix(0, S, simtime)
    ## convert splining functions to actual values
    rmat <- matrix(0, S, (simtime - 1)/deltat +1)
    for(i in 1:S){
        rmat[i,] <- predict(rspline[[i]], seq(1, simtime, by = deltat))$y
    }
    ## we keep intermediate time step abundances here
    ntmp <- n0s
    k <- 1
    for(i in 1:(simtime-1)){
        ns[,i] <- ntmp
        for(j in seq(i, i+1-deltat, by = deltat)){
            if(k > ncol(rmat)) browser()
            ntmp <- ntmp*exp(deltat*(rmat[,k] + alphas %*% ntmp))
            k <- k + 1
        }
    }
    ## fill in final abundance
    ns[,ncol(ns)] <- ntmp
    return(ns)
}
