#' A function to derive LCBD and SCBD from function beta.div
#'
#' @param df a dataframe containing species abundance/incidence with rows as sites and columns as species
#' @param ... further arguments for function beta.div
#'
#' @return a list containing two dataframes LCBD, SCBD
#' @export
#'
#' @examples
bdiv.f <- function(df, ...){
  df <- as.data.frame(df)
  row.names(df) <- df$site
  # delete site
  df <- df[, -1]
  # species must occur in at least 1 site
  df <- df[, apply(df, 2, function(x)sum(x)>0)]
  div <- beta.div(df, ...)

  SCBD <- div$SCBD %>%
    as.data.frame() %>%
    rownames_to_column("spe") %>%
    rename(SCBD = ".")

  LCBD <- div$LCBD %>%
    as.data.frame() %>%
    rownames_to_column("site") %>%
    rename(LCBD = ".") %>%
    bind_cols(p = div$p.LCBD, p.adj = div$p.adj)

  list(SCBD = SCBD, LCBD = LCBD)

}



#' Calculate taxonomic diversity
#'
#' @param spe.df a dataframe, row as site while column as species name, nonumeric data will be used in final results
#'
#' @return a dataframe, Species number, Shannon-wiener diversity, Simpson diversity and Evenness are calculated
#' @export
#'
#' @examples
cal_taxo.div <- function(spe.df){
  # species number
  spe <- spe.df[,sapply(spe.df,is.numeric)]
  S <- vegan::specnumber(spe)
  # H: Shannon-wiener diveristy,base=exp(1)
  H <- vegan::diversity(spe,index = "shannon")
  # D: Simpson diveristy
  D <- vegan::diversity(spe,index = "simpson")
  # E: eveness,Evenness is a measure of how homogeneous or even a community or ecosystem is in terms of the abundances of its species.
  # A community in which all species are equally common is considered even and has a high degree of evenness.
  E <- D/log(S)
  div <- data.frame(Shannon_div=H,Simpson_div=D,
                    Evenness=E, Species_num=S) %>%
    cbind(spe.df[,!(sapply(spe.df,is.numeric))],.)
  return(div)

}
