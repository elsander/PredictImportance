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
    toStandardize <- c('ResidualVar', 'LogDegree', 'LogCloseness',
                       'LogEigenvector', 'TrophicLevel')
    allData[,toStandardize] <- apply(allData[,toStandardize],
                                     2, function(x) return((x-mean(x))/sd(x)))
    
    if(empirical == TRUE){
        model.jacc <- lme4::lmer(LogOddsJaccard ~
                                     (-1+ResidualVar|Run) +
                                     LogDegree +
                                     LogCloseness +
                                     LogEigenvector +
                                     TrophicLevel, data = allData)
        model.jacc.varrm <- lm(LogOddsJaccard ~
                                    LogDegree +
                                    LogCloseness +
                                    LogEigenvector +
                                    TrophicLevel, data = allData)
        model.jacc.degrm <- lme4::lmer(LogOddsJaccard ~
                                            (-1+ResidualVar|Run) +
                                            LogCloseness +
                                            LogEigenvector +
                                            TrophicLevel, data = allData)
        model.jacc.closerm <- lme4::lmer(LogOddsJaccard ~
                                              (-1+ResidualVar|Run) +
                                              LogDegree +
                                              LogEigenvector +
                                              TrophicLevel, data = allData)
        model.jacc.eigrm <- lme4::lmer(LogOddsJaccard ~
                                            (-1+ResidualVar|Run) +
                                            LogDegree +
                                            LogCloseness +
                                            TrophicLevel, data = allData)
        model.jacc.tlrm <- lme4::lmer(LogOddsJaccard ~
                                           (-1+ResidualVar|Run) +
                                           LogDegree +
                                           LogCloseness +
                                           LogEigenvector, data = allData)
        
        write.table(coef(model.jacc)$Run, paste0(PathToResults, webname,
                                                 '-AllCoefficients-Removal.txt'),
                    row.names = FALSE, quote = FALSE)
        write.table(fixef(model.jacc), paste0(PathToResults, webname,
                                              '-FixedCoefficients-Removal.txt'),
                    col.names = FALSE, quote = FALSE)
    } else {
        model.jacc <- lme4::lmer(LogOddsJaccard ~
                                 (-1+ResidualVar|Web:Run) +
                                 (-1+LogDegree|Web) +
                                 (-1+LogCloseness|Web) +
                                 (-1+LogEigenvector|Web) +
                                 (-1+TrophicLevel|Web), data = allData)
        model.jacc.varrm <- lme4::lmer(LogOddsJaccard ~
                                 (-1+LogDegree|Web) +
                                 (-1+LogCloseness|Web) +
                                 (-1+LogEigenvector|Web) +
                                 (-1+TrophicLevel|Web), data = allData)
        model.jacc.degrm <- lme4::lmer(LogOddsJaccard ~
                                 (-1+ResidualVar|Web:Run) +
                                 (-1+LogCloseness|Web) +
                                 (-1+LogEigenvector|Web) +
                                 (-1+TrophicLevel|Web), data = allData)
        model.jacc.closerm <- lme4::lmer(LogOddsJaccard ~
                                 (-1+ResidualVar|Web:Run) +
                                 (-1+LogDegree|Web) +
                                 (-1+LogEigenvector|Web) +
                                 (-1+TrophicLevel|Web), data = allData)
        model.jacc.eigrm <- lme4::lmer(LogOddsJaccard ~
                                 (-1+ResidualVar|Web:Run) +
                                 (-1+LogDegree|Web) +
                                 (-1+LogCloseness|Web) +
                                 (-1+TrophicLevel|Web), data = allData)
        model.jacc.tlrm <- lme4::lmer(LogOddsJaccard ~
                                 (-1+ResidualVar|Web:Run) +
                                 (-1+LogDegree|Web) +
                                 (-1+LogCloseness|Web) +
                                 (-1+LogEigenvector|Web), data = allData)

        write.table(coef(model.jacc)$`Web:Run`, paste0(PathToResults, webname,
                                                       '-AllCoefficients-Removal.txt'),
                    quote = FALSE)
        write.table(coef(model.jacc)$Web, paste0(PathToResults, webname,
                                                 '-FixedCoefficients-Removal.txt'),
                    quote = FALSE)
    }
    save(model.jacc,
         model.jacc.varrm,
         model.jacc.degrm,
         model.jacc.closerm,
         model.jacc.eigrm,
         model.jacc.tlrm,
         file = paste0(PathToResults, webname, '-HierarchicalModel-Removal.RData'))
}
