#' Title A function that return list of file naes with specified patterns
#'
#' @param directory a string, absolute file path
#' @param pattern specified patterns e.g., "txt", "xlsx", and etc.
#'
#' @return a list, file path as the element while file name as the element name
#' @export
#'
#' @examples
get_files_dir <-function(directory=NULL,pattern=NULL){
  directory <- ifelse(is.null(directory),getwd(),directory)
  pattern <- ifelse(is.null(pattern),"txt",pattern)

  x.file <- dir(directory,pattern = pattern)
  x.dir <- paste0(directory,'/',x.file)
  names(x.dir) <- x.file
  return(x.dir)
}
