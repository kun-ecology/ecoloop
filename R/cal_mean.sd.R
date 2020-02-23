#' Calculated mean±sd
#'
#' @param df a numeric vector
#' @param digits how many digits
#' @param nsmall the minimum number of digits to the right of the decimal point in formatting real/complex numbers in non-scientific formats. Allowed values are 0 <= nsmall <= 20.
#'
#' @return a character
#' @export
#'
#' @examples usuallt it is used with dplyr::summarise and dplyr::mutate funciton
cal_mean.sd <- function(df,digits=NULL,nsmall=NULL){
  digits <- ifelse(is.null(digits),2,digits)
  nsmall <- ifelse(is.null(nsmall),2,nsmall)
  df.mean.sd <- paste0(mean(df,na.rm = T) %>% format(digits=digits,nsmall = nsmall),
                       "±",
                       sd(df,na.rm=T) %>% format(digits=digits,nsmall = nsmall))
}
