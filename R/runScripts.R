runScriptsModel <- function(model = 'Cascade', path = './'){
    setwd(path)
    system('mkdir Data')
    system('mkdir Results')
    if(model == 'Cascade'){
        system('mkdir Data/Cascade')
        Step1_Generate_Networks(web = 'Cascade')
        fname <- Step2_Discrete_LV(path = 'Data/Cascade')
        Step3_Hierarchical_Model(fname, empirical = FALSE)
    } else {
        if(model == 'Niche'){
            system('mkdir Data/Niche')
            Step1_Generate_Networks(web = 'Niche')
            fname <- Step2_Discrete_LV(path = 'Data/Niche')
            Step3_Hierarchical_Model(fname, empirical = FALSE)
        } else {
            if(model == 'MPN25'){
                system('mkdir Data/MPN25')
                Step1_Generate_Networks(web = 'MPN25')
                fname <- Step2_Discrete_LV(path = 'Data/MPN25')
                Step3_Hierarchical_Model(fname, empirical = FALSE)
            } else {
                if(model == 'MPN35'){
                    system('mkdir Data/MPN35')
                    Step1_Generate_Networks(web = 'MPN35')
                    fname <- Step2_Discrete_LV(path = 'Data/MPN35')
                    Step3_Hierarchical_Model(fname, empirical = FALSE)
                } else {
                    if(model == 'MPN45'){
                        system('mkdir Data/MPN45')
                        Step1_Generate_Networks(web = 'MPN45')
                        fname <- Step2_Discrete_LV(path = 'Data/MPN45')
                        Step3_Hierarchical_Model(fname, empirical = FALSE)
                    } else {
                        stop('model must be "Cascade", "Niche", "MPN25", "MPN35", or "MPN45"')
                    }
                }
            }
        }
    }

}

runScriptsEmpirical <- function(web = 'All'){

}
