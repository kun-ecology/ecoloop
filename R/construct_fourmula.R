#' Construct formula with multipe response and explanatory variabels
#'
#' @param Resp.Var response variables, should be character vectors
#' @param Exp.Vars explanatory variables, should be character vectors
#' @param type construct formula for each Resp.Vars and all each Exp.Vars (default) or
#' for each Resp.Vars and each Exp.Vars(set an arbitary value to type)
#' @param collapse build addititive (+) models or models with interactions (*)
#'
#' @return
#' @export
#'
#' @examples

construct_fourmula <- function(Resp.Var,Exp.Vars,type=NULL,collapse=NULL){
  collapse <- ifelse(is.null(collapse),"+","*")
  if (is.null(type)) {
    x.formula <- paste(Resp.Var,paste0(Exp.Vars,collapse =collapse),sep="~")
    x.formula <- sapply(x.formula,as.formula)
  } else {
    x.formula <- lapply(Resp.Var,function(x){
      temp.formula <- paste(x,Exp.Vars,sep = "~")
      names(temp.formula) <- paste(x,Exp.Vars,sep = "_")
      temp.formula <- lapply(temp.formula,as.formula)
    })

  }
  names(x.formula) <- Resp.Var
  return(x.formula)

}
