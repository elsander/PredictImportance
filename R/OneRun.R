OneRun <- function(matname, popname, outfile = NA, seed = NA){
    if(!is.na(seed)) set.seed(seed)
    
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
    r0s <- (-mat)%*%nstars

    simtime <- 10000
    rmat <- matrix(r0s, S, simtime)
    ##add external forcing
    rmat <- rmat + matrix(rnorm(S*simtime, mean = 0, sd = .01), S, simtime)

    ## This is the workhorse of the function.
    ## It actually carries out the simulation.
    ns <- discreteLV_R(rmat, mat, nstars, deltat = .001, simtime)
    
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
    simtime <- simtime*10
    ##external forcing for jaccard distance calculation
    rmat2 <- matrix(r0s, S, simtime) + matrix(rnorm(S*simtime, mean = 0, sd = .01), S, simtime)
    for(i in 1:S){
        nfinalstmp <- nfinals
        nfinalstmp[i] <- 0
        rmNs <- discreteLV_R(rmat2, mat, nfinalstmp, deltat = .0001, simtime)
        rmMus[,i] <- apply(rmNs[,(simtime/2):simtime], 1, mean)

        ##calculate Jaccard distance
        tmpmu <- mus
        tmpmu[i] <- 0
        jaccardvals[i] <- jaccard(tmpmu, rmMus[,i])
    }

    ##calculate degree
    adj <- (mat > 0)*1
    degs <- rowSums(adj) + colSums(adj)

    matinv <- ginv(diag(nstars, S, S)%*%mat)
    presspert <- colSums(abs(matinv))

    outdata <- rbind(jaccardvals, mus, sigmas, degs,
                     presspert, Spenser)

    if(!is.na(outfile)){
        write.table(outdata, outfile, quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
        if(sum(is.na(outfile)) > 0){
            return(outfile)
        }
    } else {
        return(outdata)
    }
}
