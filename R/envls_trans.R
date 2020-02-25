#' Transformation for a list with environmental variables in dataframe
#'
#' @param env_ls a list with environmental variables in dataframe
#' @param env.trans methods for trasnformation
#'
#' @return a list
#' @export
#'
#' @examples
envls_trans <- function(env_ls,env.trans=NULL){
  env.trans <- ifelse(is.null(env.trans),"none",env.trans)
  # only numeric data will be used
  env_ls.num <- map(env_ls,function(df){
    df.num <- df %>%
      dplyr::select_if(is.numeric)
  })
  # transformation

  env_ls.trans <- map(env_ls.num,function(df){
    if (length(env.trans)==1 & env.trans=="none"){
      df.f <- df
    } else {
      f <- match.fun(env.trans)
      df.f <- df %>%
        dplyr::mutate_all(f)
    }
  })

  return(env_ls.trans)
}
