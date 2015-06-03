#' Step 2: Simulate dynamics using discrete-time Lotka-Volterra
#'
#' Takes network parameterizations created in Step 1 and simulates community
#' dynamics using a discrete-time Lotka-Volterra. Transforms variables for
#' use in hierarchical model (Step 3).
#'
#' @param path path to folder where '-mat.txt' and '-pop.txt' files were
#'   generated in Step 1.
#' @param seed random seed, for reproducibility
#'
#' @return Returns the csv file name where results are stored. This file name
#'   can then be passed to Step 3. This results file is appended to file
#'   as the simulations are completed.
#' 
#' @export

Step2_Discrete_LV <- function(path, seed = NULL){
    fs <- list.files(path)
    inds <- grep('*mat.txt', fs)
    matfs <- fs[inds]
    ## get the root of the filename
    matsplit <- strsplit(matfs, 'mat.txt')
    matsplit <- unlist(lapply(matsplit, function(x) return(x[1])))

    ## add trailing slash if not already present
    if(!hasTrailingSlash(path)){
        path <- paste0(path, '/')
    }
    
    ## file where we will put all of the results
    filebase <- strsplit(matsplit[1], 'web.*$', perl = TRUE)[[1]][1]
    outfname <- paste0(path, filebase, 'AllData.csv')

    for(i in 1:length(matsplit)){
        print(i)
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
        tmp <- strsplit(matsplit[i], '-')
        web <- tmp[[1]][3]
        run <- tmp[[1]][5]
        vars <- c('LogOddsJaccard', 'LogPerturbation', 'LogCV', 'LogDegree')
        allData <- data.frame(Web = rep(web, nrow(out)),
                              Run = rep(run, nrow(out)),
                              out[,vars],
                              centralities)
        ## append to file
        write.table(allData, outfname, append = TRUE,
                    row.names = FALSE, quote = FALSE,
                    col.names = !file.exists(outfname),
                    sep = ',')
    }
    return(outfname)
}
