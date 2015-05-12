#' @name utilities
#' @aliases hasTrailingSlash
#' @aliases stripFileExtension
#' 
#' @title Utility functions for handling path and file names
#'
#' Check for a trailing slash '/' in a path string and strip
#' the final file extension (.txt, etc.) from a given file name
#'
#' @return hasTrailingSlash() returns 1 if mystring ends in '/', 0 otherwise.
#' stripFileExtension() returns the filepath character vector with the
#' final file extension removed.
#'
#' @rdname utilities
#' @param mystring a character vector
#'
#' @export

hasTrailingSlash <- function(mystring){
    slashTF <- grep('/$', mystring, perl = TRUE)

    ## test if a path string ends in '/'
    if(length(slashTF) > 0){
        return(1)
    } else {
        return(0)
    }
}

#' @rdname utilities
#' @param filepath a character vector containing a file name
#'
#' @export

stripFileExtension <- function(filepath){
    ## remove the final file extension from filepath
    split <- strsplit(filepath, '\\.', perl = TRUE)[[1]]
    split <- split[1:(length(split)-1)]
    file <- split[1]
    for(i in 2:length(split)){
        file <- paste0(file, '.', split[i])
    }
    return(file)
}
