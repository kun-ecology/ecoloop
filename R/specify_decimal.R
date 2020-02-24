#' Specify decimal number for numeric data
#'
#' @param x a numeric vector
#' @param k numbers of digits
#'
#' @return a character
#' @export
#'
#' @examples
specify_decimal <- function(x, k=NULL){
  k <- ifelse(is.null(k),3,k)
  x.k <- trimws(try(format(round(x, k),silent = T,nsmall = k)))
}
