#' Step 1: Generate and/or parameterize Lotka-Volterra communities
#'
#' Generates parameterized networks and abundance data that can be
#' simulated in Step 2.
#'
#' @param web Model name ('Cascade', 'Niche', or 'MPN') or empirical network
#'   name ("caricaie", "otago", "serengeti", "sylt", "ythan", "flensburg",
#'   "reef", "stmarks", or "tatoosh")
#' @param foldname Folder name where matrices and abundance vectors should
#'   be written
#' @param path Path where data are kept. 'foldname' will be created as a folder
#'   here if it does not already exist
#' @param S Number of species in the generated networks
#' @param C Network connectance, defined as 2L/(S(S-1)), where L
#'   is the number of links.
#' @param GapProb Probability of a gap in a niche. Only used if web == 'MPN'.
#' @param nwebs Number of random web structures to be generated.
#' @param nruns Number of parameterizations to be generated.
#' @param seed Random seed for reproducibility
#'
#' @return This function does not return an object, but it writes parameterized
#'   matrices (with file names ending in '-mat.txt') and vectors of equilibrium
#'   abundances (with file names ending in '-pop.txt') to the folder
#'   path/foldname.
#' 
#' @export

Step1_Generate_Networks <- function(web = 'Cascade',
                                    foldname = web,
                                    path = 'Data',
                                    S = 50,
                                    C = .1,
                                    GapProb = .25,
                                    nwebs = 30,
                                    nruns = 30,
                                    seed = NULL
                                    ){

    ## change folder name to include gap prob if appropriate
    if(web == 'MPN' && foldname == web){
        foldname <- paste0(web, GapProb*100)
    }
    
    ## for reproducibility
    set.seed(seed)

    ## standardize path format
    if(!hasTrailingSlash(path)){
        path <- paste0(path, '/')
    }
    
    system(paste0('mkdir ', path, foldname))
    for(i in 1:nwebs){
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
        for(j in 1:nruns){
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

Step1_Empirical_Parameterization <- function(web,
                                             foldname = web,
                                             path = 'Data',
                                             nruns = 30,
                                             seed = NULL){
    ## for reproducibility
    set.seed(seed)

    data(list = web)

    ## standardize path format
    if(!hasTrailingSlash(path)){
        path <- paste0(path, '/')
    }

    system(paste0('mkdir ', path, foldname))
    filebase <- paste0(path, foldname, '/', web)
    for(j in 1:nruns){
        outfile1 <- paste0(filebase, '-web-', 1, '-run-', j, '-mat.txt')
        outfile2 <- paste0(filebase, '-web-', 1, '-run-', j, '-pop.txt')
        ## parameterize
        out <- LognormalParam(as.matrix(get(web)))
        write.table(out$Mat, outfile1, row.names = FALSE, col.names = FALSE)
        write.table(out$Pop, outfile2, row.names = FALSE, col.names = FALSE)
    }
}
