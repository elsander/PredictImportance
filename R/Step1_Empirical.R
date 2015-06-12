#' Step 1: Generate and/or parameterize Lotka-Volterra communities
#'
#' Generates parameterized networks and abundance data that can be
#' simulated in Step 2.
#'
#' @param web empirical network
#'   name ("caricaie", "otago", "serengeti", "sylt", "ythan", "flensburg",
#'   "reef", "stmarks", or "tatoosh")
#' @param foldname Folder name where matrices and abundance vectors should
#'   be written
#' @param path Path where data are kept. 'foldname' will be created as a folder
#'   here if it does not already exist
#' @param nruns Number of parameterizations to be generated.
#' @param Exponential TEMPORARY
#' @param seed Random seed for reproducibility
#'
#' @return This function does not return an object, but it writes parameterized
#'   matrices (with file names ending in '-mat.txt') and vectors of equilibrium
#'   abundances (with file names ending in '-pop.txt') to the folder
#'   path/foldname.
#' 
#' @export

Step1_Empirical_Parameterization <- function(web,
                                             foldname = web,
                                             path = 'Data',
                                             nruns = 30,
                                             Exponential = FALSE,
                                             seed = NULL){
    ## for reproducibility
    set.seed(seed)

    data(list = web)

    ## standardize path format
    if(!hasTrailingSlash(path)){
        path <- paste0(path, '/')
    }

    ###################################
    if(Exponential){
        foldname <- paste0(web, '-', 'exp')
    }
    ###################################
    
    system(paste0('mkdir ', path, foldname))
    filebase <- paste0(path, foldname, '/', web)
    for(j in 1:nruns){
        outfile1 <- paste0(filebase, '-web-', 1, '-run-', j, '-mat.txt')
        outfile2 <- paste0(filebase, '-web-', 1, '-run-', j, '-pop.txt')
        ## parameterize
        ###################################
        if(Exponential){
            out <- ExponentialParam(as.matrix(get(web)))
        } else {
            out <- LognormalParam(as.matrix(get(web)))
        }
        ###################################
        write.table(out$Mat, outfile1, row.names = FALSE, col.names = FALSE)
        write.table(out$Pop, outfile2, row.names = FALSE, col.names = FALSE)
    }
}
