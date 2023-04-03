#' Specify decimal number for numeric data
#'
#' @param x a numeric vector
#' @param k numbers of digits
#' @param extra
#'
#' @return a character
#' @export
#'
#' @examples
specify_decimal <- function(x, k=NULL, extra = FALSE){


  if (extra){
     tmp <- 10^(-(k+1))*5
     x.k <- ifelse(abs(x)<tmp,
            round(x, digits = k+1) %>%  format(nsmall = k+1) %>%  trimws(),
            round(x, digits = k) %>%  format(nsmall = k) %>% trimws()
  )
  } else {
    k <- ifelse(is.null(k),3,k)
    x.k <- trimws(try(format(round(x, k),silent = T,nsmall = k)))
  }

  return(x.k)
}
