#' Utility functions for handling path and file names
#'
#' hasTrailingSlash checks for a trailing slash '/' in a path string,
#' addTrailingSlash adds a slash if it isn't already present, and
#' stripFileExtension strips the final file extension (.txt, etc.)
#' from a given file name
#'
#' @param dirpath a string containing a path to a directory
#' @param filepath a string containing a path to a file
#'
#' @return hasTrailiingSlash returns a boolean flag. addTrailingSlash
#' returns the directory path, and stripFileExtension returns the file
#' with the extension removed.
#' 
#' @export

hasTrailingSlash <- function(dirpath){
    slashTF <- grep('/$', dirpath, perl = TRUE)

    ## test if a path string ends in '/'
    if(length(slashTF) > 0){
        return(1)
    } else {
        return(0)
    }
}

addTrailingSlash <- function(dirpath){
    if(hasTrailingSlash(dirpath)){
        return(dirpath)
    } else {
        return(paste0(dirpath, '/'))
    }
}

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

