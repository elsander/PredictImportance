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
