#' Regress log(Sigma) on log(Mean) and Calculate Residuals
#'
#' Run a linear regression to find the relationship between log mean and log
#' standard deviation for one simulation of a parameterization of a network structure.
#' The residuals can then be used as a measure of variability; that is,
#' deviation from expected variability under Taylor's Law.
#'
#' @param Mean A vector of mean population sizes for a single simulation
#' of a single parameterization of a single network structure
#' @param Sigma An associated vector of standard deviations of population
#' sizes for a single simulation of a single parameterization of a single
#' network structure.
#'
#' @return if the regression is significant, a vector of residual values. If the
#' log(Sigma) ~ log(Mean) linear relationship is not significant, returns NULL.
#'
#' @export

MeanSigmaResiduals <- function(Mean, Sigma){
    lsig <- log(Sigma)
    lmean <- log(Mean)
    ## if there is a population crash during the simulation,
    ## there will be 0s in Mean/Sigma (and infinities in the log).
    ## These 0s will skew results anyway, so throw these out of the
    ## model.
    if(all(is.finite(lsig)) && all(is.finite(lmean))){
        fit <- lm(log(Sigma) ~ log(Mean))
        ## pval <- broom::tidy(fit)$p.value[2]
        pval <- summary(fit)$coefficients[2,4]
        ## if the relationship isn't significant, throw out
        ## this set of results
        if(pval <= .05){
            return(as.vector(as.matrix(fit$residuals)))
        } else {
            return(NULL)
        }
    } else {
        return(NULL)
    }
}
