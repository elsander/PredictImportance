#' Calculates a combined p-value using the weighted Z-method
#'
#' Combines a vector of one-sided p-values from independent tests into
#' a single combined two-sided p-value
#'
#' @param ps A vector of independent one-sided p-values to be
#'   combined. Note that these p-values must be one-sided to account
#'   for the fact that different tests could suggest a difference from
#'   the null in both directions, which should result in less overall
#'   evidence against the null when combined.
#'
#' @return A two-sided combined p-value.
#'
#' @examples
#' ps <- c(.05, .007, .4, .2)
#' StoufferZTransform(ps)
#'
#' @references Whitlock, M.C. 2005. Combining probability from independent
#'   tests: the weighted Z-method is superior to Fisher's approach.
#'   J. evol. biol 18(5): 1368-1373.
#'
#' @export

StoufferZTransform <- function(ps){
  ##Note: these p-values must be converted to one-sided
  ##p-values before calling this function
  ##reference: "Combining probability from independent tests:
  ##the weighted Z-method is superior to Fisher's approach"
  ##Whitlock 2005
  k <- length(ps)
  Zs <- sum(qnorm(ps))/sqrt(k)
  ##convert to pvalue by finding the probability of a more extreme Z
  ##the values must be multiplied by 2 to account for the fact that
  ##it is a two-tailed test
  newp <- min(pnorm(Zs), 1 - pnorm(Zs))*2
  return(newp)
}
