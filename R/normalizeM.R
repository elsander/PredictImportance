## Normalize the columns of the ajacency matrix M 
## by its column sums.
normalizeM <- function(M){
	colsum_M <- colSums(M);
	colsum_M[colsum_M==0] <- 1;
	return(t(t(M)/colsum_M));
}

