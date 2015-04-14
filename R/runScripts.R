runScriptStep1LogNorm <- function(seed = NULL,
                                  web = 'MPN',
                                  foldname = web,
                                  path = 'Data',
                                  GapProb = .25){
    ##web type can be "Niche", "NichePlant", "MPN", "Cascade", or "Random"
    
    ##for reproducibility
    set.seed(seed)
    
    setwd(path)
    system(paste('mkdir', foldname))
    for(i in 1:30){
        if(web == 'Niche'){
            Adj <- BuildNiche(S = 50, C = .1)
        } else {
            if(web == 'Cascade'){
                Adj <- BuildCascade(S = 50, C = .1)
            } else {
                if(web == 'MPN'){
                    Adj <- BuildMPN(S = 50, C = .1, GapProb = GapProb)
                } else {
                    stop('web argument must be "Niche", "MPN", or "Cascade"')
                }
            }
        }

        print(i)
        for(j in 1:30){
            outfile1 <- paste0(foldname, '/', foldname,
                               '-web-', i, '-run-', j, '-mat.txt')
            outfile2 <- paste0(foldname, '/', foldname,
                               '-web-', i, '-run-', j, '-pop.txt')
            out <- LognormalParam(Adj)
            write.table(out$Mat, outfile1, row.names = FALSE, col.names = FALSE)
            write.table(out$Pop, outfile2, row.names = FALSE, col.names = FALSE)
        }
    }
}

runScriptStep1Empirical <- function(webfile, foldname = webfile, seed = NULL){
    warning('This script puts equal values on the diagonal to enforce stability')
    
    ##for reproducibility
    set.seed(seed)

    system(paste0('mkdir Data/empirical/', foldname))
    Adj <- data.matrix(read.table(webfile, header = FALSE))
    for(j in 1:30){
        outfile1 <- paste0('Data/empirical/', foldname, '/', foldname,
                           '-web-', 1, '-run-', j, '-mat.txt')
        outfile2 <- paste0('Data/empirical/', foldname, '/', foldname,
                           '-web-', 1, '-run-', j, '-pop.txt')
        out <- LognormalParam(Adj)
        write.table(out$Mat, outfile1, row.names = FALSE, col.names = FALSE)
        write.table(out$Pop, outfile2, row.names = FALSE, col.names = FALSE)
    }
}

runScriptStep2 <- function(path, seed = NULL){
    setwd(path)
    fs <- list.files()
    inds <- grep('*mat.txt', fs)
    matfs <- fs[inds]
    matsplit <- strsplit(matfs, 'mat.txt')
    matsplit <- lapply(matsplit, function(x) return(x[1]))
    matsplit <- unlist(matsplit)

    for(i in 1:length(matsplit)){
        if(i%%1 == 0) print(i)
        out <- OneRun(paste0(matsplit[i], 'mat.txt'),
                      paste0(matsplit[i], 'pop.txt'),
                      paste0(matsplit[i], 'jcd.txt'),
                      seed = seed + i)
        if(!is.null(out)) print(out)
    }
}

runScriptStep3empir <- function(inpath){
    setwd(inpath)

    ##get number of species
    tmp <- strsplit(inpath, '_')
    S <- as.numeric(strsplit(tmp[[1]][2], '/'))
    
    fs <- list.files()
    matinds <- grep(paste0('-mat.txt'), fs)
    matfs <- fs[matinds]

    ##preallocate data frame
    dflen <- length(matfs)*S ##S spp per mat
    discLV <- data.frame(Web = rep(0, dflen),
                         Run = rep(0, dflen),
                         LogOddsJaccard = rep(0, dflen),
                         LogPerturbation = rep(0, dflen),
                         LogCV = rep(0, dflen),
                         LogDegree = rep(0, dflen),
                         LogCloseness = rep(0, dflen),
                         LogEigenvector = rep(0, dflen),
                         TrophicLevel = rep(0, dflen),
                         ShapeReal = rep(0, dflen),
                         ShapeSurreal = rep(0, dflen))
    k <- 1 ##indexing for the data frame
    for(i in 1:length(matfs)){
        if(!i%%25) print(i)
        matsplit <- strsplit(matfs[i], 'mat.txt')[[1]]
        jcdf <- paste0(matsplit, 'jcd.txt')
        mat <- as.matrix(read.table(matfs[i], header = FALSE))
        mat <- (mat > 0) * 1 ##convert to adjacency matrix
        gmat <- graph.adjacency(mat)

        jcd <- as.matrix(read.table(jcdf, sep = ' ', header = FALSE))
        mu <- as.vector(jcd[2,])
        sig <- as.vector(jcd[3,])
        jacc <- as.vector(jcd[1,])
        discLV$ShapeReal[k:(k+(S-1))] <- jcd[6,]
        discLV$ShapeSurreal[k:(k+(S-1))] <- jcd[7,]
        discLV$LogOddsJaccard[k:(k+(S-1))] <- log(jacc/(1-jacc)) ##log odds
        discLV$LogPerturbation[k:(k+(S-1))] <- log(as.numeric(as.vector(jcd[5,]))) ##log 
        discLV$LogCV[k:(k+(S-1))] <- log(sig/mu) ##log transformation
        deg <- as.numeric(as.vector(jcd[4,]))
        discLV$LogDegree[k:(k+(S-1))] <- log(deg) ##log odds transformation
        discLV$LogCloseness[k:(k+(S-1))] <- log(closeness(gmat)) ##log transformation
        ##Trophic Level
        ##it's the one variable that doesn't need to be transformed!
        discLV$TrophicLevel[k:(k+(S-1))] <- getTL(mat)
        ##add detritus node for eigenvector centrality
        mat2 <- rbind(rep(1,ncol(mat)), mat)
        rsum <- rowSums(mat2)
        rbool <- (rsum == 0)*1
        mat2 <- cbind(rbool, mat2)
        ##convert network to column stochastic matrix
        mat2 <- apply(mat2, 2, function(x) return(x/sum(x)))
        ##log eigenvector centrality
        discLV$LogEigenvector[k:(k+(S-1))] <- log(abs(Re(eigen(mat2)$vectors[-1,1])))

        ##get web and run number
        tmp <- strsplit(matsplit, '-')
        web <- tmp[[1]][3]
        run <- tmp[[1]][5]
        discLV$Web[k:(k+(S-1))] <- web
        discLV$Run[k:(k+(S-1))] <- run

        ##remove extinct species, if any
        extinct <- which(mu == 0)
        if(length(extinct) > 0){
            extinct <- extinct + k
            discLV[extinct,3:9] <- NA
        }
        k <- k + S
    }
    dataname <- strsplit(inpath, '/')
    dataname <- dataname[[1]][2]
    write.csv(discLV, paste0('../../', dataname, '-allData.csv'),
              row.names = FALSE, quote = FALSE)
    return(discLV)
}

runScriptStep3 <- function(inpath){
    setwd(inpath)
    fs <- list.files()
    matinds <- grep(paste0('-mat.txt'), fs)
    matfs <- fs[matinds]

    ##preallocate data frame
    dflen <- length(matfs)*50 ##50 spp per mat
    discLV <- data.frame(Web = rep(0, dflen),
                         Run = rep(0, dflen),
                         LogOddsJaccard = rep(0, dflen),
                         LogPerturbation = rep(0, dflen),
                         LogCV = rep(0, dflen),
                         LogDegree = rep(0, dflen),
                         LogCloseness = rep(0, dflen),
                         LogEigenvector = rep(0, dflen),
                         TrophicLevel = rep(0, dflen),
                         ShapeReal = rep(0, dflen),
                         ShapeSurreal = rep(0, dflen))
    k <- 1 ##indexing for the data frame
    for(i in 1:length(matfs)){
        if(!i%%25) print(i)
        matsplit <- strsplit(matfs[i], 'mat.txt')[[1]]
        jcdf <- paste0(matsplit, 'jcd.txt')
        mat <- as.matrix(read.table(matfs[i], header = FALSE))
        mat <- (mat > 0) * 1 ##convert to adjacency matrix
        gmat <- graph.adjacency(mat)

        jcd <- as.matrix(read.table(jcdf, sep = ' ', header = FALSE))
        mu <- as.vector(jcd[2,])
        sig <- as.vector(jcd[3,])
        jacc <- as.vector(jcd[1,])
        discLV$ShapeReal[k:(k+49)] <- jcd[6,]
        discLV$ShapeSurreal[k:(k+49)] <- jcd[7,]
        discLV$LogOddsJaccard[k:(k+49)] <- log(jacc/(1-jacc)) ##log odds
        discLV$LogPerturbation[k:(k+49)] <- log(as.numeric(as.vector(jcd[5,]))) ##log 
        discLV$LogCV[k:(k+49)] <- log(sig/mu) ##log transformation
        deg <- as.numeric(as.vector(jcd[4,]))
        discLV$LogDegree[k:(k+49)] <- log(deg) ##log odds transformation
        discLV$LogCloseness[k:(k+49)] <- log(closeness(gmat)) ##log transformation
        ##Trophic Level
        ##it's the one variable that doesn't need to be transformed!
        discLV$TrophicLevel[k:(k+49)] <- getTL(mat)
        ##add detritus node for eigenvector centrality
        mat2 <- rbind(rep(1,ncol(mat)), mat)
        rsum <- rowSums(mat2)
        rbool <- (rsum == 0)*1
        mat2 <- cbind(rbool, mat2)
        ##convert network to column stochastic matrix
        mat2 <- apply(mat2, 2, function(x) return(x/sum(x)))
        ##log eigenvector centrality
        discLV$LogEigenvector[k:(k+49)] <- log(abs(Re(eigen(mat2)$vectors[-1,1])))

        ##get web and run number
        tmp <- strsplit(matsplit, '-')
        web <- tmp[[1]][3]
        run <- tmp[[1]][5]
        discLV$Web[k:(k+49)] <- web
        discLV$Run[k:(k+49)] <- run

        ##remove extinct species, if any
        extinct <- which(mu == 0)
        if(length(extinct) > 0){
            extinct <- extinct + k
            discLV[extinct,3:9] <- NA
        }
        k <- k + 50
    }
    return(discLV)
}

runScriptStep4 <- function(foldername, mkdirs = TRUE){
    if(mkdirs){
        system(paste0('mkdir Results/', foldername))
        system(paste0('mkdir Results/', foldername, '/Combined'))
        system(paste0('mkdir Results/', foldername, '/Regressions'))
    }

    datafold <- paste0('Data/', foldername)
    resultfold <- paste0('Results/', foldername)
    for(i in 1:30){
        out <- allReg(datafold, resultfold, by = 'web', i)
        out <- allReg(datafold, resultfold, by = 'run', i)
    }
}

runScriptStep4empir <- function(foldername, mkdirs = TRUE){
    if(mkdirs){
        system(paste0('mkdir Results/', foldername))
        system(paste0('mkdir Results/', foldername, '/Combined'))
        system(paste0('mkdir Results/', foldername, '/Regressions'))
    }

    datafold <- paste0('Data/empirical/', foldername)
    resultfold <- paste0('Results/', foldername)
    print(datafold)
    print(resultfold)
    for(i in 1:30){
        out <- allReg(datafold, resultfold, by = 'run', i)
    }
}

runScriptStep5 <- function(foldername){
    ##get file names
    patt <- '.*jacc.*'
    patt2 <- '.*comm.*'
    path <- paste0('Results/', foldername, '/Regressions/')
    fs <- list.files(path)
    inds <- grep(patt, fs, perl = TRUE)
    inds2 <- grep(patt2, fs, perl = TRUE)
    jaccfs <- fs[inds]
    commfs <- fs[inds2]
    
    ##jaccard first
    ##initialize data frame
    dflen <- length(jaccfs)*5
    jaccdf <- data.frame(web = rep(0, dflen), var = rep(0, dflen),
                         run = rep(0, dflen), Standardized.Beta = rep(0, dflen),
                         Pvalue = rep(0, dflen))
    ##var will repeat in a predictable way, so we can save time by
    ##doing it outside of the loop
    vars <- c('log(CV)', 'log(Degree)', 'log(Closeness)', 'log(Eigenvector)', 'TL')
    jaccdf$var <- rep(vars, length(jaccfs))
    jaccdf$var <- as.factor(jaccdf$var)
    
    ##create table
    k <- 1 ##index where we are in the table
    for(i in 1:length(jaccfs)){
        jcd <- as.matrix(read.table(paste0(path, '/', jaccfs[i]),
                                    header = TRUE))
        fsplit <- strsplit(jaccfs[i], '-')[[1]]
        jaccdf$web[k:(k+4)] <- fsplit[3] #web number
        jaccdf$run[k:(k+4)] <- fsplit[5] #run number
        jaccdf$Standardized.Beta[k:(k+4)] <- as.vector(jcd[1, 1:5]) ##coefficients
        jaccdf$Pvalue[k:(k+4)] <- as.vector(jcd[2, 1:5]) ##pvalues
        k <- k+5 ##update indexing
    }
    
    ##generate violin plot
    jaccviol <- ggplot(jaccdf, aes(web, Standardized.Beta)) +
        geom_violin(fill = '#4292c6', adjust=.5) + facet_grid(var ~ .)
    ##dev.new()
    ##plot(jaccviol)
    
    pdf(paste0(path, '../Combined/', foldername, '-Standardized.Beta-violins-jacc.pdf'))
    plot(jaccviol)
    dev.off()
    
    jaccviol <- ggplot(jaccdf, aes(web, Pvalue)) +
        geom_violin(fill = '#4292c6', adjust = .5) + facet_grid(var ~ .)
    ##dev.new()
    ##plot(jaccviol)
    
    pdf(paste0(path, '../Combined/', foldername, '-Pvalue-violins-jacc.pdf'))
    plot(jaccviol)
    dev.off()
    
    ###########################################################
    ##perturbation importance now
    dflen <- length(commfs)*5
    commdf <- data.frame(web = rep(0, dflen), var = rep(0, dflen),
                         run = rep(0, dflen), Standardized.Beta = rep(0, dflen),
                         Pvalue = rep(0, dflen))
    ##var will repeat in a predictable way, so we can save time by
    ##doing it outside of the loop
    vars <- c('log(CV)', 'log(Degree)', 'log(Closeness)', 'log(Eigenvector)', 'TL')
    commdf$var <- rep(vars, length(commfs))
    commdf$var <- as.factor(commdf$var)
    
    ##create table
    k <- 1 ##index where we are in the table
    for(i in 1:length(commfs)){
        jcd <- as.matrix(read.table(paste0(path, '/', commfs[i]),
                                    header = TRUE))
        fsplit <- strsplit(commfs[i], '-')[[1]]
        commdf$web[k:(k+4)] <- fsplit[3] #web number
        commdf$run[k:(k+4)] <- fsplit[5] #run number
        commdf$Standardized.Beta[k:(k+4)] <- as.vector(jcd[1, 1:5]) ##coefficients
        commdf$Pvalue[k:(k+4)] <- as.vector(jcd[2, 1:5]) ##pvalues
        k <- k+5 ##update indexing
    }

    ##generate violin plot
    commviol <- ggplot(commdf, aes(web, Standardized.Beta)) +
        geom_violin(fill = '#4292c6', adjust = .5) + facet_grid(var ~ .)
    ##dev.new()
    ##plot(commviol)
    
    pdf(paste0(path, '../Combined/', foldername, '-Standardized.Beta-violins-comm.pdf'))
    plot(commviol)
    dev.off()
    
    commviol <- ggplot(commdf, aes(web, Pvalue)) +
        geom_violin(fill = '#4292c6', adjust = .5) + facet_grid(var ~ .)
    ##dev.new()
    ##plot(commviol)
    
    pdf(paste0(path, '../Combined/', foldername, '-Pvalue-violins-comm.pdf'))
    plot(commviol)
    dev.off()
}
