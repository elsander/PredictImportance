SimulatedViolins <- function(pathToResults){
    ## rename pathToResults to be shorter
    p <- pathToResults
    if(!hasTrailingSlash(p)){
        p <- paste0(p, '/')
    }

    ## read in hierarchical model data
    CascAllRem <- read.table(paste0(p, 'Cascade/Cascade-AllCoefficients-Removal.txt'))
    CascAllPert <- read.table(paste0(p,
                              'Cascade/Cascade-AllCoefficients-Perturbation.txt'))
    NicheAllRem <- read.table(paste0(p, 'Niche/Niche-AllCoefficients-Removal.txt'))
    NicheAllPert <- read.table(paste0(p,
                               'Niche/Niche-AllCoefficients-Perturbation.txt'))
    MPN25AllRem <- read.table(paste0(p, 'MPN25/MPN25-AllCoefficients-Removal.txt'))
    MPN25AllPert <- read.table(paste0(p,
                               'MPN25/MPN25-AllCoefficients-Perturbation.txt'))
    MPN35AllRem <- read.table(paste0(p, 'MPN35/MPN35-AllCoefficients-Removal.txt'))
    MPN35AllPert <- read.table(paste0(p,
                               'MPN35/MPN35-AllCoefficients-Perturbation.txt'))
    MPN45AllRem <- read.table(paste0(p, 'MPN45/MPN45-AllCoefficients-Removal.txt'))
    MPN45AllPert <- read.table(paste0(p,
                               'MPN45/MPN45-AllCoefficients-Perturbation.txt'))

    
    CascFixedRem <- read.table(paste0(p,
                               'Cascade/Cascade-FixedCoefficients-Removal.txt'))
    CascFixedPert <- read.table(paste0(p,
                                'Cascade/Cascade-FixedCoefficients-Perturbation.txt'))
    NicheFixedRem <- read.table(paste0(p, 'Niche/Niche-FixedCoefficients-Removal.txt'))
    NicheFixedPert <- read.table(paste0(p,
                                 'Niche/Niche-FixedCoefficients-Perturbation.txt'))
    MPN25FixedRem <- read.table(paste0(p, 'MPN25/MPN25-FixedCoefficients-Removal.txt'))
    MPN25FixedPert <- read.table(paste0(p,
                                 'MPN25/MPN25-FixedCoefficients-Perturbation.txt'))
    MPN35FixedRem <- read.table(paste0(p, 'MPN35/MPN35-FixedCoefficients-Removal.txt'))
    MPN35FixedPert <- read.table(paste0(p,
                                 'MPN35/MPN35-FixedCoefficients-Perturbation.txt'))
    MPN45FixedRem <- read.table(paste0(p, 'MPN45/MPN45-FixedCoefficients-Removal.txt'))
    MPN45FixedPert <- read.table(paste0(p,
                                 'MPN45/MPN45-FixedCoefficients-Perturbation.txt'))

    ##make data frame for ggplot
    grouplist <- c('Cascade', 'Niche', 'MPN25', 'MPN35', 'MPN45')
    webvec <- c(rep(grouplist, each = 900),
                rep(grouplist, 4, each = 30))
    varvec <- c(rep('log(CV)', 900*length(grouplist)),
                rep(c('log(Degree)', 'Trophic Level',
                      'log(Closeness)', 'log(Eigenvector)'),
                    each = 30*length(grouplist)))
    rcoeff <- c(CascAllRem$LogCV,
                NicheAllRem$LogCV,
                MPN25AllRem$LogCV,
                MPN35AllRem$LogCV,
                MPN45AllRem$LogCV,
                CascFixedRem$LogDegree,
                NicheFixedRem$LogDegree,
                MPN25FixedRem$LogDegree,
                MPN35FixedRem$LogDegree,
                MPN45FixedRem$LogDegree,
                CascFixedRem$TrophicLevel,
                NicheFixedRem$TrophicLevel,
                MPN25FixedRem$TrophicLevel,
                MPN35FixedRem$TrophicLevel,
                MPN45FixedRem$TrophicLevel,
                CascFixedRem$LogCloseness,
                NicheFixedRem$LogCloseness,
                MPN25FixedRem$LogCloseness,
                MPN35FixedRem$LogCloseness,
                MPN45FixedRem$LogCloseness,
                CascFixedRem$LogEigenvector,
                NicheFixedRem$LogEigenvector,
                MPN25FixedRem$LogEigenvector,
                MPN35FixedRem$LogEigenvector,
                MPN45FixedRem$LogEigenvector)

    pcoeff <- c(CascAllPert$LogCV,
                NicheAllPert$LogCV,
                MPN25AllPert$LogCV,
                MPN35AllPert$LogCV,
                MPN45AllPert$LogCV,
                CascFixedPert$LogDegree,
                NicheFixedPert$LogDegree,
                MPN25FixedPert$LogDegree,
                MPN35FixedPert$LogDegree,
                MPN45FixedPert$LogDegree,
                CascFixedPert$TrophicLevel,
                NicheFixedPert$TrophicLevel,
                MPN25FixedPert$TrophicLevel,
                MPN35FixedPert$TrophicLevel,
                MPN45FixedPert$TrophicLevel,
                CascFixedPert$LogCloseness,
                NicheFixedPert$LogCloseness,
                MPN25FixedPert$LogCloseness,
                MPN35FixedPert$LogCloseness,
                MPN45FixedPert$LogCloseness,
                CascFixedPert$LogEigenvector,
                NicheFixedPert$LogEigenvector,
                MPN25FixedPert$LogEigenvector,
                MPN35FixedPert$LogEigenvector,
                MPN45FixedPert$LogEigenvector)

    violdf <- data.frame(web = webvec,
                         var = varvec,
                         RemCoeff = rcoeff,
                         PertCoeff = pcoeff)

    ## get ggplot to order categories correctly
    violdf$web <- ordered(violdf$web,
                          levels =  c('Niche', 'MPN25',
                              'MPN35', 'MPN45', 'Cascade'))

    ##generate violin plot
    remviol <- ggplot(violdf, aes(web, RemCoeff)) +
        geom_violin(fill = '#4292c6', adjust=.5) + facet_grid(. ~ var) +
        coord_flip() + scale_y_continuous(limits = c(-1.5,1.5)) + 
        theme_bw() + ylab('Regression Coefficient') + xlab('') +
        geom_hline(yintercept = 0, alpha = .5)

    
    pdf('Simulated-Violins-Removal.pdf', height = 2.25, width = 7)
    plot(remviol)
    dev.off()
    
    pertviol <- ggplot(violdf, aes(web, PertCoeff)) +
        geom_violin(fill = '#4292c6', adjust=.5) + facet_grid(. ~ var) +
        coord_flip() + scale_y_continuous(limits = c(-1.5,1.5)) + 
        theme_bw() + ylab('Regression Coefficient') + xlab('') +
        geom_hline(yintercept = 0, alpha = .5)

    
    pdf('Simulated-Violins-Perturbation.pdf', height = 2.25, width = 7)
    plot(pertviol)
    dev.off()
}

## Demo script
readline('This code will generate figures, data, and model results in the current directory. Press any key to continue:')

readline("Simulation and analysis for the Cascade model will run now. This may take a couple of days! Press any key to continue:")

runScriptsModel(model = 'Cascade')

print('Niche')
##readline("Simulation and analysis for the Niche model will run now. This may take a couple of days! Press any key to continue:")

runScriptsModel(model = 'Niche')

##readline("Simulation and analysis for the MPN (25% gap prob) model will run now. This may take a couple of days! Press any key to continue:")

print('MPN25')
runScriptsModel(model = 'MPN25')

##readline("Simulation and analysis for the MPN (35% gap prob) model will run now. This may take a couple of days! Press any key to continue:")

print('MPN35')
runScriptsModel(model = 'MPN35')

##readline("Simulation and analysis for the MPN (45% gap prob) model will run now. This may take a couple of days! Press any key to continue:")

print('MPN45')
runScriptsModel(model = 'MPN45')

print('Generating violin plots from all results')
SimulatedViolins('Results')
