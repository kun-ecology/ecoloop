#' Title Customed log1p for data with negative value less than -1
#'
#' @param x a numeric vector
#' @param base exp(1) or 10, or 2 based log (x+a)
#'
#' @return a dataframe
#' @export
#'
#' @examples
my_log1p <- function(x,base=exp(1)){
  base <- ifelse(base==exp(1),base,10)
  min <- min(x,na.rm = T)
  if (min <=(-1)){
    b <- abs(round(min))
    return(log(x+b+1,base = base))
  } else {
    b <- 1
    log(x+b,base = base)
  }
}
