#' Calculate Bray-Curtis similarity and Jaccard distance
#' between two ecological communities
#'
#' Calculate Bray-Curtis similarity and Jaccard distance between
#' two real-valued vectors of abundances. Bray-Curtis similarity is
#' calculated as sum(|x-y|/|x+y|), where x and y are real-valued vectors.
#' and Jaccard similarity is calculated as 2BC/(1+BC), where BC
#' is the Bray-Curtis similarity.
#'
#' @param x A vector of species abundances.
#' @param y A second vector of abundances, with the same length as x.
#'
#' @return The similarity or distance between the two abundance vectors.
#'
#' @export

braycurtis <- function(x, y){
    return(sum(abs(x-y))/sum(abs(x+y)))
}

jaccard <- function(x, y){
    return(2*braycurtis(x,y)/(1+braycurtis(x,y)))
}
