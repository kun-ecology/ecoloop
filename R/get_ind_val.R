#' get_ind.val: compute indicator species analsyis
#'
#' @param spe a data frame wirh row as observations (sites) while column as species
#' @param grp  a vector indicating the group of, it should be a factor
#' @param numitr number of iteration,deafult numitr=999
#' @param p.cut  significant levels, deafult p.cut=0.05
#' @param p.adj indicate whethe adjusted p value shoul be used, by default T
#' @param ... other param for function indval
#'
#' @return a dataframe with 5 columns i.e., species names (spe), group (group),
#' indicator values (indval), pvalue, freq
#' @export
#'
#' @examples
get_ind_val <- function(spe,grp,numitr=NULL,p.cut=NULL,p.adj=NULL,...){
  numitr <- ifelse(is.null(numitr),999,numitr)
  p.cut <- ifelse(is.null(p.cut),0.05,p.cut)
  p.adj <- ifelse(is.null(p.adj),T, p.adj)
  # only numeric data in spe will be considered as species data
  spe.df <- dplyr::select_if(spe,is.numeric)

  iva <- labdsv::indval(spe.df,grp,numitr=numitr,...)
  if (p.adj==T){
    pval.adj <- p.adjust(iva$pval)
  } else {
    pval.adj <- iva$pval
  }
  gr <- levels(grp)[iva$maxcls[pval.adj <= p.cut]]
  iv <- iva$indcls[pval.adj <= p.cut]
  pv <- iva$pval[pval.adj <= p.cut]
  fr <- apply(spe.df > 0, 2, sum)[pval.adj <= p.cut]
  fidg <- data.frame(
    group = gr,
    indval = iv,
    pvalue = pv,
    freq = fr
  )
  fidg <- fidg[order(fidg$group, -fidg$indval), ] %>%
    mutate(spe=row.names(.)) %>%
    select(5,1:4)

  return(list(indval.mod=iva,indval.df=fidg))


}
