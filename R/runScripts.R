#' @name runScripts
#' @aliases runScriptsModel
#' @aliases runScriptsEmpirical
#' 
#' @title Run all analyses for a food web model or network
#'
#' Run all data simulation and analysis steps for a given model or
#' empirical web(s). This may take a while!
#'
#' @param path path to where Data and Results folders are or should be created
#' 
#' @rdname runScripts
#' @param model food web model used to construct the network: options are
#'   'Cascade', 'Niche', 'MPN25', 'MPN35', and 'MPN45'.
#'
#' @examples
#' \dontrun{
#' runScriptsModel('Niche')
#' }
#'
#' @export

runScriptsModel <- function(model = 'Cascade', path = './'){
    old <- getwd()
    on.exit(setwd(old), add = TRUE)
    
    setwd(path)
    system('mkdir Data')
    system('mkdir Results')
    print(getwd())
    
    if(!(model %in% c('Cascade', 'Niche', 'MPN25', 'MPN35', 'MPN45')) ||
       length(model) > 1){
        stop("model must be 'Cascade', 'Niche', 'MPN25', 'MPN35', or 'MPN45'.")
    }

    system(paste0('mkdir Data/', model))

    ## get GapProb, if necessary
    GapProb <- .25
    if(model %in% c('MPN25', 'MPN35', 'MPN45')){
        GapProb <- as.numeric(strsplit(model, 'MPN')[[1]][2])/100
        model2 <- 'MPN'
    } else {
        model2 <- model
    }
        
    Step1_Generate_Networks(web = model2, GapProb = GapProb)
    fname <- Step2_Discrete_LV(path = paste0('Data/', model))
    Step3_Hierarchical_Model(fname, empirical = FALSE)
}

#' @rdname runScripts
#'
#' @param web empirical network to be parameterized and simulated: options are
#'   "caricaie", "otago", "serengeti", "sylt", "ythan", "flensburg",
#'   "reef", "stmarks", and "tatoosh". A vector of network names can be
#'   provided to simulate more than one network, or web = "All" can be used
#'   to simulate all networks in the package.
#'
#' @examples
#' \dontrun{
#' runScriptsEmpirical('tatoosh')
#' }
#'
#' @export

runScriptsEmpirical <- function(web = 'All', path = './'){
    setwd(path)
    system('mkdir Data')
    system('mkdir Results')
    nets <- c('caricaie', 'otago', 'serengeti', 'sylt', 'ythan',
              'flensburg', 'reef', 'stmarks', 'tatoosh')
    
    if(!all(web %in% nets) && web != 'All'){
        stop('web argument must be "All" or a vector containing one or more of the following: "caricaie", "otago", "serengeti", "sylt", "ythan", "flensburg", "reef", "stmarks", or "tatoosh".')
    }
    
    if(web != 'All'){
        nets <- web
    }

    ## run the three steps for each network in the list
    for(net in nets){
        system(paste0('mkdir Data/', net))
        Step1_Empirical_Parameterization(net)
        out <- Step2_Discrete_LV(paste0('Data/',net))
        Step3_Hierarchical_Model(out)
    }
}
