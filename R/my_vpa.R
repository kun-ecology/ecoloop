#' Conduct VPA after CCA or RDA
#'
#' @param sp a dataframe, row as sites while column as species names
#' @param env_ls a list, with environmental variables data as elements
#' @param sp.trans a string, method for species data transformation
#' @param env.trans a string, method for envrionmental data transformation
#' @param mode a string, when force.mode=F (default), it indicates which model to use when 3<dca.axlength<4
#' when force.mode=T, it indicates which mode to conduct forcely
#' @param force.mode T or F (default)
#' @param ... further argumens for adespatial::forward.sel
#'
#' @return a list
#' @export
#'
#' @examples
my_vpa <- function(sp, env_ls, sp.trans=NULL, env.trans=NULL,mode=NULL,
                      force.mode=NULL,...){
  sp.trans <- ifelse(is.null(sp.trans), "none",sp.trans)
  env.trans <- ifelse(is.null(env.trans),"none", env.trans)
  force.mode <- ifelse(is.null(force.mode),F,T)
  # transformation for sp
  sp.f <- spM_trans(sp,sp.trans=sp.trans)

  # add names to env_ls
  if (is.null(names(env_ls))){
    names(env_ls) <- paste0("env",1:length(env_ls))
  }
  env_ls.f <-envls_trans(env_ls = env_ls, env.trans=env.trans)

  #############################################################
  # if mode is not set, then DCA will be conducted to decided whihc one will be used
  if (force.mode==F){
    # DCA analysis
    sp.dca <- sp.f %>%  vegan::decorana()
    # Axis length of the four axis
    sp.dca.axleng <- sp.dca$rproj %>%
      as.data.frame() %>%
      map_dbl(max) %>%
      max()
    # if sp.dca.axleng >= 4, then CCA model will be used
    # if 3<= sp.dca.axleng <4, then both CCA and RDA are fine to use
    # and RDA model will be used if model is not set
    # if 3<sp.dca.axleng, then RDA will be used.
    if (sp.dca.axleng>=4){
      mode <- "cca"
    } else if (sp.dca.axleng>=3){
      if (is.null(mode)){
        mode <- "rda"
      } else {
        mode <- "cca"
      }
    } else {
      mode <- "rda"
    }
    #print the mode will be used
    # print the axis length
    cat("=============================== \n")
    cat("Final model was selcted based on DCA: \n")
    cat(paste("Max of axis length in DCA is: ",sp.dca.axleng %>% specify_decimal(2),"\n"))
    cat(paste("Mode used: ", mode, "\n"))
    cat("=============================== \n")


  } else { # when force.mode is true,
    sp.dca <- NA
    sp.dca.axleng <- NA
    if (is.null(mode)){
      stop("When force.mode = T, either rda or cca should be assigned to mode")
    }
    cat("=============================== \n")
    cat(paste("Mode used: ", mode, "\n"))
  }
  #############################################################
  # conduct the mode
  mode.f <- match.fun(mode)
  sp_ls.mod <- try(map(env_ls.f,function(df){
    mod.f <- match.fun(mode)
    sp.mod <- try(mod.f(sp.f,df),silent = T)
  }),silent = T)

  if (ncol(sp.f)>2000){
    cat("=============================== \n")
    cat("Extract RsquareAdj (May take minutes if species number is high) \n")
    cat("=============================== \n")
  }
  # anova with the mode
  sp_ls.mod.aov <- try(map(sp_ls.mod,vegan::anova.cca),silent = T)
  sp_ls.mod.p <- try(map_dbl(sp_ls.mod.aov,function(mod.aov)mod.aov$`Pr(>F)`[1]),silent = T)
  sp_ls.mod.sig <- try(sp_ls.mod[sp_ls.mod.p<0.05],silent = T) # signifcant model

  # adjusted R2
  sp_ls.modR2 <- map_df(sp_ls.mod,vegan::RsquareAdj) %>%
    bind_cols(env.type=names(sp_ls.mod),.) %>%
    bind_cols(mod.pvalue=sp_ls.mod.p)
  ######################################################
  # model forward seleciton
  # only do forward selection with significant mod
  none.sel <- data.frame(variables="none",order=NA, R2=NA, R2Cum=NA, AdjR2Cum=NA, F= NA, pvalue=NA)

  if (length(sp_ls.mod.sig)==0){ # no models is significant
    cat("=============================== \n")
    cat("Warning: none of the mode is significant \n")
    cat("=============================== \n")
    mode.sel <- list(none=NA) # forward selection mode
    var.sel <- list(none=none.sel)
    vpa <- NA
  } else { # model is significant
    mode.sel <- names(sp_ls.mod.sig) %>%
      map(function(nm){
        cat("=============================== \n")
        cat(paste("Forward selection procedure for: ",nm, "\n"))
        md.sel <- try(adespatial::forward.sel(sp.f,env_ls.f[[nm]],...),silent = T)

      })
    names(mode.sel) <- names(sp_ls.mod.sig)
  }

  # only select valid mode.sel
  mode.sel <- discard(mode.sel, function(mod)class(mod)=="try-error") %>%
    discard(function(x)all(is.na(x)))

  # #####################################################
  ### extract selected variables
  if (length(mode.sel)==0) { # no variables was selected
    cat("=============================== \n")
    cat("Results of forward selections \n")
    cat("Warning: no variables are selected \n")
    cat("=============================== \n")
    mode.sel <- list(none=none.sel)
    var.sel <- list(none="none")
    vpa <- NA
  } else if (sum(!is.na(mode.sel))<2){ # variables are selected but only 1 data matrix
    cat("=============================== \n")
    cat("Warning: two to four explanatory variables were needed for varpart \n")
    cat("=============================== \n")
    var.sel <- map(mode.sel, function(mod){
      mod$variables %>% as.vector()
    })

    var.sel.df.ls <- names(var.sel) %>%
      map(function(nm){
        env_ls.f[[nm]] %>%
          dplyr::select(var.sel[[nm]])
      })
    names(var.sel.df.ls) <- names(var.sel)
    vpa <- NA

  } else {

    var.sel <- map(mode.sel, function(mod){
      mod$variables %>% as.vector()
    })

    var.sel.df.ls <- names(var.sel) %>%
      map(function(nm){
        env_ls.f[[nm]] %>%
          dplyr::select(var.sel[[nm]])
      })
    names(var.sel.df.ls) <- names(var.sel)
    vpa.expr <- paste0("varpart(","sp.f,",
                       paste0("env_ls.f[[",1:length(var.sel.df.ls),"]]",collapse = ","),
                       ")")
    vpa <- eval(rlang::parse_expr(vpa.expr))

  }
  #########################################################


  res.f <- list(rawdata=list(sp=sp,env_ls=env_ls),
                mode.info=list(sp.trans=sp.trans,env.trans=env.trans,sp.dca=list(dca=sp.dca)),
                mode.res=list(mode=sp_ls.mod,mode.aov=sp_ls.mod.aov,mode.aov.p=sp_ls.mod.p,
                              mode.sig=sp_ls.mod.sig,
                              mode.adjR2=sp_ls.modR2),
                var.sel=list(mode.sel=mode.sel, var.sel=var.sel),
                vpa=vpa)
  attr(res.f,"label") <- "community-based VPA"
  return(res.f)



}
