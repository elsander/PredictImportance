Step1_Generate_Networks <- function(web = 'MPN',
                                    foldname = web,
                                    path = 'Data',
                                    seed = NULL,
                                    S = 50,
                                    C = .1,
                                    GapProb = .25){
    ##web type can be "Niche", "NichePlant", "MPN", "Cascade", or "Random"
    
    ##for reproducibility
    set.seed(seed)

    system(paste('mkdir', foldname))
    for(i in 1:30){
        if(web == 'Niche'){
            Adj <- BuildNiche(S = S, C = C)
        } else {
            if(web == 'Cascade'){
                Adj <- BuildCascade(S = S, C = C)
            } else {
                if(web == 'MPN'){
                    Adj <- BuildMPN(S = S, C = C, GapProb = GapProb)
                } else {
                    stop('web argument must be "Niche", "MPN", or "Cascade"')
                }
            }
        }

        print(i)
        for(j in 1:30){
            outfile1 <- paste0(path, foldname, '/', foldname,
                               '-web-', i, '-run-', j, '-mat.txt')
            outfile2 <- paste0(path, foldname, '/', foldname,
                               '-web-', i, '-run-', j, '-pop.txt')
            out <- LognormalParam(Adj)
            write.table(out$Mat, outfile1, row.names = FALSE, col.names = FALSE)
            write.table(out$Pop, outfile2, row.names = FALSE, col.names = FALSE)
        }
    }
}

Step1_Empirical_Parameterization <- function(webfile,
                                             foldname = webfile,
                                             seed = NULL){
    ## for reproducibility
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

Step2_Discrete_LV <- function(path, seed = NULL){

    ## for reproducibility
    set.seed(seed)
    
    fs <- list.files(path)
    inds <- grep('*mat.txt', fs)
    matfs <- fs[inds]
    ## get the root of the filename
    matsplit <- strsplit(matfs, 'mat.txt')
    matsplit <- unlist(lapply(matsplit, function(x) return(x[1])))

    ## file where we will put all of the results
    outfname <- paste0(path, matsplit, '-AllData.csv')

    for(i in 1:length(matsplit)){
        if(i%%1 == 0) print(i)
        matfile <- paste0(path, matsplit[i], 'mat.txt')
        out <- OneRun(matfile,
                      paste0(path, matsplit[i], 'pop.txt'))

        ## transform variables
        out$LogOddsJaccard <- log(out$Jaccard/(1-out$Jaccard))
        out$LogPerturbation <- log(out$Perturbation)
        out$LogCV <- log(out$Sigma/out$Mean)
        out$LogDegree <- log(out$Degree)
        
        ## get network centralities
        centralities <- GetCentralities(matfile)
        
        ## get web and run number
        tmp <- strsplit(matsplit, '-')
        web <- tmp[[1]][3]
        run <- tmp[[1]][5]
        vars <- c('LogOddsJaccard', 'LogPerturbation', 'LogCV', 'LogDegree')
        allData <- data.frame(Web = rep(web, ncol(out)),
                              Run = rep(run, ncol(out)),
                              vars,
                              centralities)
        ## append to file
        write.csv(allData, outfname, append = TRUE,
                  row.names = FALSE, quote = FALSE,
                  col.names = !file.exists(outfname))
    }
    return(outfname)
}

Step4_Hierarchical_Model <- function(datafname,
                                     PathToResults = 'Results',
                                     empirical = TRUE){
    ## run hierarchical model
    allData <- read.csv(datafname)
    allData$Run <- as.factor(allData$Run)
    allData$Web <- as.factor(allData$Web)

    ## extract web name
    webname <- strsplit(datafname, '/')[[1]]
    webname <- webname[length(webname)]
    webname <- strsplit(webname, '-AllData.csv')[[1]][1]

    ## standardize variables
    toStandardize <- c('LogCV', 'LogDegree', 'LogCloseness',
                       'LogEigenvector', 'TrophicLevel')
    empir[,toStandardize] <- apply(empir[,toStandardize],
                                   2, function(x) return((x-mean(x))/sd(x)))
    
    if(empirical == TRUE){
        model.jacc <- lmer(LogOddsJaccard ~
                           (-1+LogCV|Run) +
                           LogDegree +
                           LogCloseness +
                           LogEigenvector +
                           TrophicLevel, data = allData)
        
        model.pert <- lmer(LogPerturbation ~
                           (-1+LogCV|Run) +
                           LogDegree +
                           LogCloseness +
                           LogEigenvector +
                           TrophicLevel, data = allData)
    } else {
        model.jacc <- lmer(LogOddsJaccard ~
                           (-1+LogCV|Web:Run) +
                           (-1+LogDegree|Web) +
                           (-1+LogCloseness|Web) +
                           (-1+LogEigenvector|Web) +
                           (-1+TrophicLevel|Web), data = allData)

        model.pert <- lmer(LogPerturbation ~
                           (-1+LogCV|Web:Run) +
                           (-1+LogDegree|Web) +
                           (-1+LogCloseness|Web) +
                           (-1+LogEigenvector|Web) +
                           (-1+TrophicLevel|Web), data = allData)
    }

    write.table(coef(model.jacc)$Run, paste0(PathToResults, webname,
                                             '-AllCoefficients-Removal.txt'),
                row.names = FALSE, quote = FALSE)
    write.table(fixef(model.jacc), paste0(PathToResults, webname,
                                          '-FixedCoefficients-Removal.txt'),
                col.names = FALSE, quote = FALSE)
    save(model.jacc, file = paste0(PathToResults, webname,
                         '-HierarchicalModel-Removal.RData'))

    write.table(coef(model.pert)$Run, paste0(PathToResults, webname,
                                             '-AllCoefficients-Perturbation.txt'),
                row.names = FALSE, quote = FALSE)
    write.table(fixef(model.pert), paste0(PathToResults, webname,
                                          '-FixedCoefficients-Perturbation.txt'),
                col.names = FALSE, quote = FALSE)
    save(model.pert, file = paste0(PathToResults, webname,
                         '-HierarchicalModel-Perturbation.RData'))
}
