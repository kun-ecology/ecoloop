#' Species matrix transformation
#'
#' @param sp species abundance data, row as sites while column as species names.
#' Non-numeric data are allowed but will be droped in the final results
#' @param sp.trans methods for tranformation, allowed methods are: methods in vegan::decostand, log1p, log and
#' other self composed functions.
#'
#' @return a dataframe
#' @export
#'
#' @examples
spM_trans <- function(sp, sp.trans=NULL){
  sp.trans <- ifelse(is.null(sp.trans),"none",sp.trans)
  sp.num <- dplyr::select_if(sp,is.numeric)
  # method can be arguments in vegan::decostand
  decostand.method <-  c("total","max","frequency","normalize","range","rank",
                         "standardize","pa","chi.square","hellinger","log")
  if (sp.trans=="none"){
    sp.df <- sp
  } else if (sp.trans %in% decostand.method){
    sp.df <- sp.num %>%
      vegan::decostand(method = sp.trans)
  } else { # other transformation method i.e., log1p, sqrt etc.
    f <- match.fun(sp.trans)
   sp.df <- sp.num %>%
     mutate_all(f)
  }

  return(sp.df)
}
