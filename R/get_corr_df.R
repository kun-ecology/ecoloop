#' A function return correlation and p value in dataframe formate
#'
#' @param df.x a dataframe, only numeric data will be used
#' @param df.y a dataframe, not necessary
#' @param type pearson (by defacult) or spearman correlation to be used
#'
#' @return
#' @export
#'
#' @examples
get_corr_df <- function(df.x,df.y=NULL,type=NULL){
  df.x <- dplyr::select_if(df.x,is.numeric)
  type <- ifelse(is.null(type),"pearson","spearman")

  specify_decimal <- function(x, k=NULL){
    if (is.null(k)) k=3
    else k=k
    trimws(format(round(x, k), nsmall=k))
  }

  if (is.null(df.y)){
    df.corr <- Hmisc::rcorr(df.x %>% as.matrix(),type=type)
    x.nm <- y.nm <- colnames(df.x)
  } else {
    df.y <- dplyr::select_if(df.y,is.numeric)
    df.corr <- Hmisc::rcorr(df.x %>% as.matrix(),df.y %>% as.matrix(),type=type)
    x.nm <-  colnames(df.x)
    y.nm <- colnames(df.y)
  }
  df.cor <- df.corr$r %>% as.data.frame() %>%
    mutate(x.nm=row.names(.)) %>%
    gather(key="y.nm",value = "Corr",1:(ncol(.)-1))
  df.cor.p <- df.corr$P %>% as.data.frame() %>%
    mutate(x.nm=row.names(.)) %>%
    gather(key="y.nm",value = "P",1:(ncol(.)-1))
  df.cor <- left_join(df.cor,df.cor.p,by=c("x.nm","y.nm"))

  # mutate(P=ifelse(is.na(P),0,P)) %>%
  # mutate(Corr=ifelse(P<=p.cut,Corr,NA))
  if (is.null(df.y)){
    x.df.dist <- combn(x.nm,2) %>%
      t() %>%
      as.data.frame()
    names(x.df.dist) <- c("x.nm","y.nm")
    df.cor.f <- x.df.dist %>%
      left_join(df.cor,by=c("x.nm","y.nm")) %>%
      mutate_if(is.numeric,specify_decimal)
  } else {
    df.cor.f <- data.frame(x.nm=rep(x.nm,each=length(y.nm)),y.nm=rep(y.nm,length(x.nm))) %>%
      left_join(df.cor,by=c("x.nm","y.nm")) %>%
      mutate_if(is.numeric,specify_decimal) %>%
      filter(x.nm!=y.nm)
  }

  df.cor.f$Corr <- as.numeric(df.cor.f$Corr)
  return(df.cor.f)
}
