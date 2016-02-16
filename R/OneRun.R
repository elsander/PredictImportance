#' Simulate Lotka-Volterra population dynamics and get summary stats
#'
#' Read in a matrix and population file, simulate discrete-time
#' population dynamics with perturbed growth rates, and calculate
#' species importance, variability, and degree. Note that this function
#' is used internally, and many parameters (simtime, perturbation size)
#' are hardcoded. Use discreteLV_C() for your own simulations or
#' Step2_Discrete_LV() to replicate simulations from Wootton et al.
#'
#' @param matname path to parameterized network to be simulated
#' @param popname path to vector of equilibrium abundances for the network
#' @param simtime number of time steps to simulate
#' @param deltat step size for simulation (smaller value will allow the model to
#' more closely approximate a continuous-time Lotka-Volterra)
#'
#' @return a data frame containing Jaccard distance, perturbation importance,
#' mean abundance, standard deviation of abundance over the simulation, and
#' network degree for each species in the network.
#'
#' @export

OneRun <- function(matname, popname, simtime = 1000, deltat = .001){
    ## hardcoded jaccard simtime and burn-in
    jaccardSimtime <- 500
    burnIn <- 300
    
    ##matrix of alphas
    ## read in, if applicable
    if(is.character(matname)){
        mat <- as.matrix(read.table(matname, header = FALSE))
    } else {
        mat <- as.matrix(matname)
    }
    S <- dim(mat)[1]
    mat <- t(mat) ##to make it consistent with how we're simulating it

    ##equilibrium ns
    ##read in, if applicable
    if(is.character(popname)){
        nstars <- as.vector(as.matrix(read.table(popname, header = FALSE)))
    } else {
        nstars <- as.vector(as.matrix(popname))
    }
    r0s <- (-mat) %*% nstars

    rmat <- matrix(r0s, S, simtime)
    ##add external forcing
    rmat <- rmat + matrix(rnorm(S*simtime, mean = 0, sd = .01), S, simtime)

    ## This is the workhorse of the function.
    ## It actually carries out the simulation.
    ns <- discreteLV_C(rmat, mat, nstars, deltat, simtime)
    
    extinctionthreshold <- 10e-06
    nfinals <- ns[,dim(ns)[2]]
    nfinals[nfinals < extinctionthreshold] <- 0

    ##average abundance for simulation
    mus <- apply(ns[,(simtime/2):simtime], 1, mean)
    ##standard deviation for simulation
    sigmas <- apply(ns[,(simtime/2):simtime], 1, sd)

    ##remove each species to calculate effect on community
    rmMus <- matrix(0, S, S)
    jaccardvals <- rep(0, S)
    ##external forcing for jaccard distance calculation
    rmat2 <- matrix(r0s, S, jaccardSimtime) +
        matrix(rnorm(S*jaccardSimtime, mean = 0, sd = .01), S, jaccardSimtime)
    for(i in 1:S){
        nfinalstmp <- nfinals
        nfinalstmp[i] <- 0
        rmNs <- discreteLV_C(rmat2, mat, nfinalstmp, deltat, jaccardSimtime)
        rmNs[rmNs < extinctionthreshold] <- 0
        rmMus[,i] <- apply(rmNs[,burnIn:jaccardSimtime], 1, mean)

        ##calculate Jaccard distance
        tmpmu <- mus
        tmpmu[i] <- 0
        jaccardvals[i] <- jaccard(tmpmu, rmMus[,i])
    }

    ##calculate degree
    adj <- (mat > 0)*1
    degs <- rowSums(adj) + colSums(adj)

    matinv <- solve(diag(nstars, S, S)%*%mat)
    presspert <- colSums(abs(matinv))

    outdata <- data.frame(Jaccard = jaccardvals,
                          Perturbation = presspert,
                          Mean = mus,
                          Sigma = sigmas,
                          Degree = degs)

    return(outdata)
}
