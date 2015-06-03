loadModelScript <- function(){
    source(system.file(package = 'predictimportance',
                       'scripts/FullSimulationModels.R'))
}

loadEmpiricalScript <- function(){
    source(system.file(package = 'predictimportance',
                       'scripts/FullSimulationEmpirical.R'))
}
