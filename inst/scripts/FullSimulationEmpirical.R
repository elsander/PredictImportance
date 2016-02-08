EmpiricalViolins <- function(pathToResults){
    webnames <- c('caricaie',
                  'flensburg',
                  'otago',
                  'reef',
                  'serengeti',
                  'stmarks',
                  'sylt',
                  'tatoosh',
                  'ythan')

    ##initialize data frame vectors
    webfixed <- webRand <- character(0)
    RemFixed <- RemRand <- numeric(0)
    PertFixed <- PertRand <- numeric(0)
    RemSE <- numeric(0)
    PertSE <- numeric(0)

    for(i in 1:length(webnames)){
        web <- webnames[i]
        print(web)
        
        ##read in files
        runrem <- read.table(paste0(pathToResults, '/', web, '/', web,
                                    '-AllCoefficients-Removal.txt'),
                              header = TRUE)
        runpert <- read.table(paste0(pathToResults, '/', web, '/', web,
                                     '-AllCoefficients-Perturbation.txt'),
                              header = TRUE)
        load(paste0(pathToResults, '/', web, '/', web,
                    '-HierarchicalModel-Removal.RData'))
        load(paste0(pathToResults, '/', web, '/', web,
                    '-HierarchicalModel-Perturbation.RData'))
        
        serem <- sqrt(diag(vcov(model.jacc)))
        sepert <- sqrt(diag(vcov(model.pert)))

        ##Information about random effects (only CV)
        webRand <- c(webRand, rep(web, length(runrem$LogCV)))
        RemRand <- c(RemRand, runrem$LogCV)
        PertRand <- c(PertRand, runpert$LogCV)
        
        ##Information about fixed effects (everything but CV)
        webfixed <- c(webfixed, rep(web, 4))
        RemFixed <- c(RemFixed, as.vector(as.matrix(runrem[1,3:6])))
        PertFixed <- c(PertFixed, as.vector(as.matrix(runpert[1,3:6])))
        RemSE <- c(RemSE, serem[2:5])
        PertSE <- c(PertSE, sepert[2:5])
    }

    ##data frame for pointrange facets (fixed effects)
    sedf <- data.frame(web = webfixed,
                       var = rep(c('log(Degree)', 'log(Closeness)',
                           'log(Eigenvector)', 'Trophic Level'), length(RemSE)/4),
                       RemCoeff = RemFixed,
                       PertCoeff = PertFixed,
                       RemSE = RemSE,
                       PertSE = PertSE)

    ##data frame for violin facets (random effects)
    violdf <- data.frame(web = webRand,
                         var = rep('log(CV)', length(webRand)),
                         RemCoeff = RemRand,
                         PertCoeff = PertRand)

    ## plots violins for removal importance with error bars around estimates
    empdf <- data.frame(web = c(webfixed, webRand),
                        var = as.factor(c(as.character(sedf$var),
                            as.character(violdf$var))),
                        RemCoeff = c(RemFixed, RemRand),
                        PertCoeff = c(PertFixed, PertRand))

    remviol <- ggplot(empdf, mapping = aes(x = web, y = RemCoeff)) + coord_flip()
    remviol <- remviol + facet_grid(. ~ var) +
        scale_y_continuous(limits = c(-3,3)) + 
        theme_bw() + ylab('Standardized Regression Coefficient') + xlab('') +
            geom_hline(yintercept = 0, alpha = .5)
    remviol <- remviol + geom_violin(data = violdf, fill = '#4292c6', adjust=1.5)
    remviol <- remviol + geom_pointrange(data = sedf, fill = '#4292c6',
                                           aes(x = web, y = RemCoeff, 
                                               ymax = RemCoeff + RemSE,
                                               ymin = RemCoeff - RemSE),
                                           size = .3)
    pdf('Empirical-Violins-Removal-errorbars.pdf', height = 2, width = 8)
    plot(remviol)
    dev.off()

####################Perturbation plot#########################
    ## same plot, but for perturbation importance
    pertviol <- ggplot(empdf, mapping = aes(x = web, y = PertCoeff)) + coord_flip()
    pertviol <- pertviol + facet_grid(. ~ var) +
        scale_y_continuous(limits = c(-3,3)) + 
        theme_bw() + ylab('Standardized Regression Coefficient') + xlab('') +
        geom_hline(yintercept = 0, alpha = .5)
    pertviol <- pertviol + geom_violin(data = violdf, fill = '#4292c6', adjust=1.5)
    pertviol <- pertviol + geom_pointrange(data = sedf, fill = '#4292c6',
                                           aes(x = web, y = PertCoeff, 
                                               ymax = PertCoeff + PertSE,
                                               ymin = PertCoeff - PertSE),
                                           size = .3)
    pdf('Empirical-Violins2-Perturbation-errorbars.pdf', height = 2, width = 8)
    plot(pertviol)
    dev.off()
}

FullSimulationEmpirical <- function(){
    readline('This code will generate figures, data, and model results in the current directory. Simulation and analysis for all empirical networks will run. This may take multiplehours or days! Press any key to continue:')

    runScriptsEmpirical(web = 'All', path = './')

    print('Generating violin plots from empirical results')
    EmpiricalViolins('Results')
}
