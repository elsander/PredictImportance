multReg <- function(matfile, jcdfile, outpath){
    mat <- as.matrix(read.table(matfile, header = FALSE))
    mat <- (mat > 0) * 1 ##convert to adjacency matrix
    gmat <- graph.adjacency(mat)

    jcd <- as.matrix(read.table(jcdfile, header = FALSE))
    
    jacc <- as.vector(jcd[1,])
    ##print(jacc)
    jacc <- log(jacc/(1-jacc)) ##log odds transformation
    comm <- log(as.vector(jcd[5,])) ## log transformation
    mu <- as.vector(jcd[2,])
    sig <- as.vector(jcd[3,])
    CV <- log(sig/mu) ##log transformation
    CV[!is.finite(CV)] <- NA ##put NAs where mu is 0
    deg <- as.vector(jcd[4,])
    deg <- log(deg) ##log odds transformation
    close <- log(closeness(gmat)) ##log transformation

    ##Trophic Level
    ##it's the one variable that doesn't need to be transformed!
    TL <- getTL(mat)
    
    ##add detritus node for eigenvector centrality
    mat2 <- rbind(rep(1,ncol(mat)), mat)
    rsum <- rowSums(mat2)
    rbool <- (rsum == 0)*1
    mat2 <- cbind(rbool, mat2)
    ##convert network to column stochastic matrix
    mat2 <- apply(mat2, 2, function(x) return(x/sum(x)))
    ##log eigenvector centrality
    eig <- log(abs(Re(eigen(mat2)$vectors[-1,1])))

    lmjacc <- lm(jacc ~ CV + deg + close + eig + TL, na.action = na.omit)
    anovajacc <- anova(lmjacc)
    coeffs <- as.vector(lmjacc$coefficients[2:6])
    ##standardize coeffs
    coeffs[1] <- coeffs[1]*sd(CV, na.rm=TRUE)/sd(jacc, na.rm=TRUE)
    coeffs[2] <- coeffs[2]*sd(deg, na.rm=TRUE)/sd(jacc, na.rm=TRUE)
    coeffs[3] <- coeffs[3]*sd(close, na.rm=TRUE)/sd(jacc, na.rm=TRUE)
    coeffs[4] <- coeffs[4]*sd(eig, na.rm=TRUE)/sd(jacc, na.rm=TRUE)
    coeffs[5] <- coeffs[5]*sd(TL, na.rm=TRUE)/sd(jacc, na.rm=TRUE)
    ##print(anovaCV)
    pvals <- as.vector(anovajacc[[5]][1:5])

    lmcomm <- lm(comm ~ CV + deg + close + eig + TL, na.action = na.omit)
    anovacomm <- anova(lmcomm)
    commcoeffs <- as.vector(lmcomm$coefficients[2:6])
    ##standardize commcoeffs
    commcoeffs[1] <- commcoeffs[1]*sd(CV, na.rm=TRUE)/sd(comm, na.rm=TRUE)
    commcoeffs[2] <- commcoeffs[2]*sd(deg, na.rm=TRUE)/sd(comm, na.rm=TRUE)
    commcoeffs[3] <- commcoeffs[3]*sd(close, na.rm=TRUE)/sd(comm, na.rm=TRUE)
    commcoeffs[4] <- commcoeffs[4]*sd(eig, na.rm=TRUE)/sd(comm, na.rm=TRUE)
    commcoeffs[5] <- commcoeffs[5]*sd(TL, na.rm=TRUE)/sd(comm, na.rm=TRUE)
    ##print(anovaCV)
    commpvals <- as.vector(anovacomm[[5]][1:5])

    jaccoutmat <- rbind(coeffs, pvals)
    colnames(jaccoutmat) <- c('log(CV)', 'log(degree)', 'log(closeness)',
                              'log(eigenvector)', 'TrophicLevel')
    commoutmat <- rbind(commcoeffs, commpvals)
    colnames(commoutmat) <- c('log(CV)', 'log(degree)', 'log(closeness)',
                              'log(eigenvector)', 'TrophicLevel')

    jcdfile <- strsplit(jcdfile, '/')
    jcdfile <- jcdfile[[1]][length(jcdfile[[1]])]
    write.table(jaccoutmat, paste0(outpath, '/Regressions/', jcdfile, '.jaccreg'),
                quote = FALSE, row.names = FALSE)
    write.table(commoutmat, paste0(outpath, '/Regressions/', jcdfile, '.commreg'),
                quote = FALSE, row.names = FALSE)
    return(list(jacc = jaccoutmat, comm = commoutmat))
}

allReg <- function(inpath, outpath, by = c('web', 'run'), num){
    ##path is the path to the mat and jcd files
    ##we have multiple runs for multiple files. by chooses whether we regress
    ##all runs of one web (by='web'), or all webs of one run (by='run').
    ##num is the run number (if by='run') or web number (if by='web').
    
    fs <- list.files(inpath)
    if(by == 'web'){
        matinds <- grep(paste0('.*web-', num, '-run-.*-mat.txt'), fs)
    } else {
        matinds <- grep(paste0('.*-run-', num, '-mat.txt'), fs)
    }
    matfs <- fs[matinds]

    jacccoeffs <- matrix(0, length(matfs), 5)
    jaccpvals <- matrix(0, length(matfs), 5)
    commcoeffs <- matrix(0, length(matfs), 5)
    commpvals <- matrix(0, length(matfs), 5)
    
    for(i in 1:length(matfs)){
        matsplit <- strsplit(matfs[i], 'mat.txt')[[1]]
        jcdf <- paste0(matsplit, 'jcd.txt')
        results <- multReg(paste0(inpath,'/',matfs[i]),
                           paste0(inpath,'/',jcdf),
                           outpath)
        jacccoeffs[i,] <- results$jacc[1,]
        jaccpvals[i,] <- results$jacc[2,]
        commcoeffs[i,] <- results$comm[1,]
        commpvals[i,] <- results$comm[2,]
    }

    ##convert pvals to one-sided pvals
    ##we expect importance to be negatively correlated with CV, so negative -> low p
    ##we expect importance to be positively correlated with others, so positive -> low p
    jacccoeffsTF <- jacccoeffs < 0
    jaccpvals[,1][jacccoeffsTF[,1]] <- jaccpvals[,1][jacccoeffsTF[,1]]/2
    jaccpvals[,1][!jacccoeffsTF[,1]] <- 1 - jaccpvals[,1][!jacccoeffsTF[,1]]/2
    jaccpvals[,2:5][jacccoeffsTF[,2:5]] <- 1 - jaccpvals[,2:5][jacccoeffsTF[,2:5]]/2
    jaccpvals[,2:5][!jacccoeffsTF[,2:5]] <- jaccpvals[,2:5][!jacccoeffsTF[,2:5]]/2
    jaccsigTF <- jaccpvals < .05
    ##calculate combined pvals
    jaccfinalpvals <- apply(jaccpvals, 2, StoufferZTransform)
    jaccCoeffAvg <- apply(jacccoeffs, 2, mean)
    jaccCoeffSD <- apply(jacccoeffs, 2, sd)
    jaccpctsig <- apply(1*jaccsigTF, 2, mean)

    jaccmytable <- rbind(jaccCoeffAvg, jaccCoeffSD, jaccfinalpvals, jaccpctsig)
    jaccmytable <- as.data.frame(jaccmytable)
    names(jaccmytable) <- c('logCV', 'logDegree', 'logCloseness', 'logEigenvector', 'TrophicLevel')
    row.names(jaccmytable) <- c('Std.Beta', 'SD', 'Combined-p', 'Prop.Significant')

    ##now for the perturbation importance
    ##convert pvals to one-sided pvals
    ##we expect importance to be negatively correlated with CV, so negative -> low p
    ##we expect importance to be positively correlated with others, so positive -> low p
    commcoeffsTF <- commcoeffs < 0
    commpvals[,1][commcoeffsTF[,1]] <- commpvals[,1][commcoeffsTF[,1]]/2
    commpvals[,1][!commcoeffsTF[,1]] <- 1 - commpvals[,1][!commcoeffsTF[,1]]/2
    commpvals[,2:5][commcoeffsTF[,2:5]] <- 1 - commpvals[,2:5][commcoeffsTF[,2:5]]/2
    commpvals[,2:5][!commcoeffsTF[,2:5]] <- commpvals[,2:5][!commcoeffsTF[,2:5]]/2
    commsigTF <- commpvals < .05
    ##calculate combined pvals
    commfinalpvals <- apply(commpvals, 2, StoufferZTransform)
    commCoeffAvg <- apply(commcoeffs, 2, mean)
    commCoeffSD <- apply(commcoeffs, 2, sd)
    commpctsig <- apply(1*commsigTF, 2, mean)

    commmytable <- rbind(commCoeffAvg, commCoeffSD, commfinalpvals, commpctsig)
    commmytable <- as.data.frame(commmytable)
    names(commmytable) <- c('logCV', 'logDegree', 'logCloseness', 'logEigenvector', 'TrophicLevel')
    row.names(commmytable) <- c('Std.Beta', 'SD', 'Combined-p', 'Prop.Significant')

    prefix <- strsplit(jcdf, '-')[[1]][1]
    if(by == 'web'){
        write.table(jaccmytable, paste0(outpath, '/Combined/',
                                        prefix, '-allruns-web-',
                                        num, '.jacc'), quote = FALSE)
        write.table(commmytable, paste0(outpath, '/Combined/',
                                        prefix, '-allruns-web-',
                                        num, '.comm'), quote = FALSE)
    } else {
        write.table(jaccmytable, paste0(outpath, '/Combined/',
                                        prefix, '-allwebs-run-',
                                        num, '.jacc'), quote = FALSE)
        write.table(commmytable, paste0(outpath, '/Combined/',
                                        prefix, '-allwebs-run-',
                                        num, '.comm'), quote = FALSE)
    }
}
