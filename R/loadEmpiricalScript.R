#' Parameterize empirical networks, run models, and plot results as in Wootton et al.
#'
#' Parameterizes empirical networks 30 times, simulates
#' population dynamics, runs hierarchical model, and generates violin plots
#' for the following food webs: caricaie, flensburg, otago, reef, serengeti,
#' stmarks, sylt, tatoosh, and ythan. This replicates main text results
#' in Wootton et al. (in prep). Note that this simulation and analysis will likely
#' days to complete!
#'
#' @export

loadEmpiricalScript <- function(){
    source(system.file(package = 'predictimportance',
                       'scripts/FullSimulationEmpirical.R'))
}
