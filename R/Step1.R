Step1_Generate_Networks <- function(web = 'MPN',
                                    foldname = web,
                                    path = 'Data',
                                    seed = NULL,
                                    S = 50,
                                    C = .1,
                                    GapProb = .25,
                                    nwebs = 30,
                                    nruns = 30){
    ## web type can be "Niche", "NichePlant", "MPN", "Cascade", or "Random"
    
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

Step1_Empirical_Parameterization <- function(webfile,
                                             nruns = 30,
                                             seed = NULL){
    ## for reproducibility
    set.seed(seed)

    ## get base of webfile name to use for outfile names
    filebase <- stripFileExtension(webfile)
    
    Adj <- data.matrix(read.table(webfile, header = FALSE))
    for(j in 1:nruns){
        outfile1 <- paste0(filebase, '-web-', 1, '-run-', j, '-mat.txt')
        outfile2 <- paste0(filebase, '-web-', 1, '-run-', j, '-pop.txt')
        ## parameterize
        out <- LognormalParam(Adj)
        write.table(out$Mat, outfile1, row.names = FALSE, col.names = FALSE)
        write.table(out$Pop, outfile2, row.names = FALSE, col.names = FALSE)
    }
}
