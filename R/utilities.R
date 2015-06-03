## Utility functions for handling path and file names
##
## hasTrailingSlash: Check for a trailing slash '/' in a path string and
##
## stripFileExtension: strip the final file extension (.txt, etc.)
## from a given file name
##

hasTrailingSlash <- function(mystring){
    slashTF <- grep('/$', mystring, perl = TRUE)

    ## test if a path string ends in '/'
    if(length(slashTF) > 0){
        return(1)
    } else {
        return(0)
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
