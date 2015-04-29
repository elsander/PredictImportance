Step3_Hierarchical_Model <- function(datafname,
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
