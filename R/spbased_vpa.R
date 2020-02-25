#' Conduct species based VPA
#'
#' @param sp a dataframe, row as sites while column as species names
#' @param env_ls a list, with environmental variables data as elements
#' @param sp.trans a string, method for species data transformation
#' @param env.trans a string, method for envrionmental data transformation
#' @param mode when species based VPA are conducted mode="rda"
#' @param force.mode force.mode=T, when species based VPAs are conducted
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
spbased_vpa <- function(sp, env_ls, sp.trans=NULL, env.trans=NULL,mode="rda",force.mode=T,...){
  sp.f.ls <- map(1:ncol(sp),function(x){
    sp %>%
      dplyr::select(x) %>%
      spM_trans(sp.trans=sp.trans)
  })
  names(sp.f.ls) <- colnames(sp)
  cat("=============================== \n")
  cat("Species based VPA: \n")
  cat("=============================== \n")
  n <- length(sp.f.ls)
  spbased.vpa.f <- map(1:n,function(x){
    cat(paste0("Processing ", names(sp.f.ls)[n], "(", x," otu of ",n, ") \n"))
    cat("=============================== \n")
    df <- sp.f.ls[[x]]
    my_vpa(sp=df, env_ls = env_ls, sp.trans = "none",env.trans = env.trans,
           mode = mode, force.mode=force.mode)

  })
  names(spbased.vpa.f) <- names(sp.f.ls)
  attr(spbased.vpa.f,"label") <- "species based VPA"
  return(spbased.vpa.f)

}
