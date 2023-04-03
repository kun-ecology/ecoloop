#' variation partitioning for linear mixed model
#'
#' @param full_mod a model results returned by **nlme::lme**
#'
#' @return a data frame containing explained variation by each variable
#' @export
#'
#' @details see https://www.pnas.org/doi/abs/10.1073/pnas.1608980113 for details of how this algorithm works.
#' "The function r.squaredGLMM from the MuMIn-package (67) was used to calculate the marginal (fixed effects) and conditional
#'  (full model) R2 for all optimal and beyond-optimal models. We then followed a variation partitioning procedure to determine
#'  the relative proportion of variation for each response variable in the models with elevation, GDD, and plot productivity.
#'  To perform our variation partitioning, we constructed a series of models with (i) only one focal variable,
#'  (ii) all variables except that focal variable, or (iii) the full model with all explanatory variables.
#'  The proportion of variation explained by each fixed factor was represented by calculating the difference between
#'  the marginal R2 of the full model and of the model without the focal variable and dividing it by the marginal R2
#'   of the full model (68). For all factors in all models, the variation explained by the model without the focal variable
#'   (R2 of full model minus R2 of model with only the focal variable, divided by full model) and the shared variance
#'    (R2 of model of focal variable plus R2 of model without focal variable minus R2 of full model, divided by R2 of full model)
#'    were calculated as well (SI Appendix, Tables S7–S9). This variance partitioning could not be performed for the models
#'    of probability of flower production, because this probability was zero in the undisturbed subplots in the Scandes."
#'
#' @examples
lme_var.part <- function(full_mod){
  require(MuMIn)
  require(nlme)
  df <- full_mod$data
  r2 <- r.squaredGLMM(full_mod) %>% as.data.frame()
  # terms in the formula
  all.terms <- attr(full_mod$terms, "factors") %>% row.names()
  # at least two explanatory variables
  if(length(all.terms)==2) stop("At least 2 explanatory variables are required")

  # random var
  rnd.var <- names(full_mod$modelStruct$reStruct)
  # rep var
  resp.var <- all.terms[1]
  # exp var
  exp.var <- all.terms[2:length(all.terms)]

  # models for each variables
  mod.ls <- imap(exp.var, function(nm,n){
    # only focal varible
    focal.fmla <- paste0("lme(" ,resp.var, "~", nm, ", random= ~1|", rnd.var, ", data=df)")
    # only focal varible
    focal.ex.fmla <- paste0("lme(" ,resp.var, "~", paste(exp.var[-n], collapse = "+"),
                            ", random= ~1|",
                            rnd.var,", df)")
    list(include=focal.fmla, exclude=focal.ex.fmla) %>%
      map(rlang::parse_expr) %>%
      map(eval, envir = parent.env(current_env()))
  })
  names(mod.ls) <- exp.var
  mod.r2 <- map(mod.ls, ~ map(.x, r.squaredGLMM)) %>%
    map(~do.call(rbind.data.frame,.x)) %>%
    do.call(rbind.data.frame,.) %>%
    bind_cols(exp.var=word(row.names(.), 1, sep="\\."),
              focalOrNot=word(row.names(.), 2, sep="\\."),.) %>%
    bind_rows(data.frame(exp.var="all", focalOrNot="all", R2m=r2$"R2m", R2c=r2$"R2c")) %>% as_tibble()
  var.perc <- mod.r2 %>%
    pivot_wider(id_cols = exp.var, names_from = focalOrNot, values_from = c(R2m)) %>%
    mutate(all=sum(all, na.rm = T)) %>%
    filter(exp.var!="all") %>%
    mutate(perc=(all-exclude)/all)
  res <- list(r2=mod.r2, perc=var.perc)
  invisible(res)
}



#' extract explained variation for each nested level in lmer model
#'
#' @param mod a model results returned by **lme4::lmer**
#'
#' @return a dataframe contains explained variation for each level
#' @export
#'
#' @details see https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2435.13723 for details
#' "We applied a linear mixed model to partition the variance in root traits with a given trait as the
#' dependent factor and only random effects (‘elevation’ (i.e. conspecific plants between elevations) nested in ‘species’ nested in ‘growth form’].
#' The remaining variance was explained by trait differences between conspecific plants growing
#' in different replicate plots at the same elevation (Albert et al., 2010)."
#'
#' note that the random terms are ranked in this way: (1|Level3/Level2/Level1). Additionally, the order of these
#' terms do not affect the results.
#'
#'
#' @examples
#'
rand_var.part <- function(mod){
  # Get the variance explained per level:
  var.expl = as.data.frame(VarCorr(mod))
  # Calculate the total variance explained
  Tvar = sum(var.expl$vcov)
  # Calculate the % of the total variance explained by each level
  tmp <- var.expl$vcov
  var.expl <- map_dbl(tmp, ~ (.x/Tvar)*100)
  res <- data.frame(lvl= c(paste("Lvl", 1:(length(tmp)-1), sep="_"), "residual"),
                    terms = word(names(mod@cnms), sep = ":", 1),
                    var.explained = var.expl
  )
  return(res)
}

