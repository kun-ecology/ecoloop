#' Group taxa (mainly for OTU) into dominant and rare
#'
#' @param df species abundance data, see https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-019-0749-8#Sec28
#' for how the grouping things works
#' @param lower 0.0001
#' @param upper 0.01
#' @param simplify by default simpplify=F, which will return AAT, CAT, CRAT, MT, CRT, ART; otherwise only dominant and rare taxa will return
#'
#' @return
#' @export
#'
#' @examples
group_otu <- function(df,lower=NULL, upper=NULL, simplify=T){
  df.num <- select_if(df, is.numeric)
  if (any(df.num>1)){
    df.num.f <- df.num/rowSums(df.num)
    #df.num.f <- bind_cols(select_if(df, function(x)!is.numeric(x)),df.num.f)
    message(paste("Abundance is converted to relative abundance","\n"))
  } else {
    df.num.f <- df.num
  }

  lower <- ifelse(is.null(lower),0.0001,lower)
  upper <- ifelse(is.null(upper),0.01,upper)
  # Always abundant taxa (AAT), OTU with relative abundance always >=1% in all sample
  AAT <- df.num.f %>%
    select_if(function(x)all(x>=upper))
  # Conditionally abundant taxa (CAT) were defined as the OTUs with a relative abundance greater than 0.01% in all samples and ≥ 1% in some samples but never rare (< 0.01%)
  CAT <- df.num.f %>%
    select_if(function(x) all(x>lower) & any(x>=upper))
  # Conditionally rare and abundant taxa (CRAT) were defined as the OTUs with a relative abundance varying from rare (< 0.01%) to abundant (≥ 1%)
  CRAT <- df.num.f %>%
    select_if(function(x) any(x<lower) & any(x>=upper))
  # Moderate taxa (MT) were defined as the OTUs with relative abundance between 0.01% and 1% in all samples.
  MT <- df.num.f %>%
    select_if(function(x) all(x>lower) & all(x<upper))
  #  Conditionally rare taxa (CRT) were defined as the OTUs with a relative abundance < 0.01% in some samples but never ≥ 1% in all samples.
  CRT <- df.num.f %>%
    select_if(function(x) any(x>lower) & !any(x>=upper))
  # Always rare taxa (ART) were defined as the OTUs with a relative abundance always < 0.01% in all samples.
  ART <- df.num.f %>%
    select_if(function(x) all(x<lower))
  if(simplify==T){
    final <- list(AAT=AAT,CAT=CAT,CRAT=CRAT,MT=MT, CRT=CRT, ART=ART)
  } else {
    final1 <- list(AAT,CAT, CRAT) %>% do.call(bind_cols,.) #dominant
    final2 <- list(MT, CRT, ART) %>% do.call(bibind_cols,.) # rare
    final <- list(dominant=final1, rare=final2)
  }
  return(final)

}
