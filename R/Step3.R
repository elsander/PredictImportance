#' Step 3: Analyze data using a hierarchical model
#'
#' Takes data file created in Step 2 and analyzes the data using a hierarchical model
#'
#' @param datafname path to file generated in Step 2 (file name will end
#'   in '-AllData.csv')
#' @param PathToResults path to location where results should be stored.
#' @param empirical flags if the simulated data are based on an empirical
#' network structure. This will affect the hierarchical model structure
#'
#' @return This function does not return an object, but it writes hierarchical
#'   model output to text files, and saves the entire lmer object to a .RData file
#' 
#' @export

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

    ## add webname directory in Results if it doesn't exist already
    system(paste0('mkdir ', PathToResults, '/', webname))

    ## add slash to path if necessary
    if (!hasTrailingSlash(PathToResults)){
        PathToResults <- paste0(PathToResults, '/')
    }
    PathToResults <- paste0(PathToResults, webname, '/')
    
    ## standardize variables
    toStandardize <- c('LogCV', 'LogDegree', 'LogCloseness',
                       'LogEigenvector', 'TrophicLevel')
    allData[,toStandardize] <- apply(allData[,toStandardize],
                                     2, function(x) return((x-mean(x))/sd(x)))
    
    if(empirical == TRUE){
        model.jacc <- lme4::lmer(LogOddsJaccard ~
                                 (-1+LogCV|Run) +
                                 LogDegree +
                                 LogCloseness +
                                 LogEigenvector +
                                 TrophicLevel, data = allData)
        
        model.pert <- lme4::lmer(LogPerturbation ~
                                 (-1+LogCV|Run) +
                                 LogDegree +
                                 LogCloseness +
                                 LogEigenvector +
                                 TrophicLevel, data = allData)
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
    } else {
        model.jacc <- lme4::lmer(LogOddsJaccard ~
                                 (-1+LogCV|Web:Run) +
                                 (-1+LogDegree|Web) +
                                 (-1+LogCloseness|Web) +
                                 (-1+LogEigenvector|Web) +
                                 (-1+TrophicLevel|Web), data = allData)

        model.pert <- lme4::lmer(LogPerturbation ~
                                 (-1+LogCV|Web:Run) +
                                 (-1+LogDegree|Web) +
                                 (-1+LogCloseness|Web) +
                                 (-1+LogEigenvector|Web) +
                                 (-1+TrophicLevel|Web), data = allData)
        write.table(coef(model.jacc)$`Web:Run`, paste0(PathToResults, webname,
                                                       '-AllCoefficients-Removal.txt'),
                    quote = FALSE)
        write.table(coef(model.jacc)$Web, paste0(PathToResults, webname,
                                                 '-FixedCoefficients-Removal.txt'),
                    quote = FALSE)
        save(model.jacc, file = paste0(PathToResults, webname,
                             '-HierarchicalModel-Removal.RData'))
        write.table(coef(model.pert)$`Web:Run`, paste0(PathToResults, webname,
                                                       '-AllCoefficients-Perturbation.txt'),
                    quote = FALSE)
        write.table(coef(model.pert)$Web, paste0(PathToResults, webname,
                                                 '-FixedCoefficients-Perturbation.txt'),
                    quote = FALSE)
        save(model.pert, file = paste0(PathToResults, webname,
                             '-HierarchicalModel-Perturbation.RData'))
    }
}
