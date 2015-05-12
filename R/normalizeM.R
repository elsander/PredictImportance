#' Normalize matrix columns
#' 
#' Normalize the columns of the ajacency matrix M 
#' by its column sums. This is a utility function for
#' predictimportance which will likely not be useful
#' except internally.
#'
#' @param M a matrix to be normalized
#'
#' @return the matrix M with columns normalized by column sum
#'
#' @export

normalizeM <- function(M){
	colsum_M <- colSums(M);
	colsum_M[colsum_M==0] <- 1;
	return(t(t(M)/colsum_M));
}

