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
    fit <- lm(log(Sigma) ~ log(Mean))
    pval <- broom::tidy(fit)$p.value[2]
    if(pval > .05){
        return(NULL)
    } else {
        return(as.vector(as.matrix(fit$residuals)))
    }
}
