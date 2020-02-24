#' Community- and species- based variation partitioning
#'
#' @param sp species abundance data with sites as row while species as columns, raw data is needed if you want to run species-based VPA
#' @param env_ls a list containing at least two elements, while at most four elements;
#' meaningful name are highly recommended, otherwise env1,env2 and etc will be used as the name
#' @param sp.trans transformation for spe data
#' @param env.trans transformation methos for env data, can be a vector with length 1 or with same length with env-ls
#' @param mode which mode to use i.e., RDA or CCA, while do community- based VPA,
#'  the mode is chosed based on the max value of axis length in DCA analysis, i.e., axis length>4, CCA;
#'  3<= axis length <= 4, either RDA or CCA is fine, by default, here use RDA; axis length<3, RDA
#'
#' @param sp.based species-based VPA or not, by default community-based VPA will be conducted
#' @param ... arguments for forward.sel in package adespatial
#'
#' @return a list including raw data, model information, model results
#' @export
#'
#' @examples
my_vpa <- function(sp,env_ls,sp.trans=NULL,env.trans=NULL,mode=NULL,sp.based=NULL,...){
  require(tidyverse)
  require(vegan)
  require(adespatial)
  sp.trans <- ifelse(is.null(sp.trans),"none",sp.trans) # sp.trans should be one of the method in vegan::decostand
  env.trans <- ifelse(is.null(env.trans),"none",env.trans)
  sp.based <- ifelse(is.null(sp.based),F,sp.based)

  # only numeric data from sp data will be used
  sp.num <- sp %>% dplyr::select_if(is.numeric)

  specify_decimal <- function(x, k=NULL){
    k <- ifelse(is.null(k),3,k)
    x.k <- trimws(try(format(round(x, k),silent = T,nsmall = k)))
  }
  # transformation for sp
  cat("=============================== \n")
  cat(paste("Species data transformation: ", sp.trans, "\n"))
  cat("=============================== \n")
  if (sp.trans=="none"){
    sp.h <- sp.num
  } else {
    sp.h <- sp.num %>%
      vegan::decostand(method = sp.trans)
  }

  # initiate names for env list
  if (is.null(names(env_ls))){
    env_ls.nm <- paste0("env",1:length(env_ls))
  } else {
    env_ls.nm <- names(env_ls)
  }

  # only numeric data from env data will be used
  env_ls <- map(env_ls,function(df){
    df %>% dplyr::select_if(is.numeric)
  })
  # transformation for env
  cat("=============================== \n")
  cat("Environmental data transformation: \n")
  cat(paste("Transformation for",env_ls.nm,": ",env.trans, "\n"))
  cat("============================== \n")

  if(env.trans=="none"){
    env_ls.f <- env_ls
  } else if (all(env.trans!="none")){
    env_ls.f <- map2(env_ls,env.trans,function(df,f){
      df.f <- df %>%
        dplyr::select_if(is.numeric) %>%
        mutate(f)
      return(df.f)
    })
  }

  # a function for conducting model and do forward selection
  con.mod.sel <- function(sp.h,env_ls.f, mode){
    # conduct the mode
    sp.ls.mod <- try(map(env_ls.f,function(df){
      mod.f <- match.fun(mode)
      sp.mod <- try(mod.f(sp.h,df),silent = T)
    }),silent = T)

    # anova with the mode
    sp.ls.mod.aov <- try(map(sp.ls.mod,vegan::anova.cca),silent = T)
    sp.ls.mod.p <- try(map_dbl(sp.ls.mod.aov,function(mod.aov)mod.aov$`Pr(>F)`[1]),silent = T)
    sp.ls.mod.sig <- try(sp.ls.mod[sp.ls.mod.p<0.05],silent = T) # signifcant model

    # only do forward selection with significant mod
    # do forward selection
    if (length(sp.ls.mod.sig)==0){
      cat("=============================== \n")
      cat("Warning: none of the mode is significant \n")
      cat("=============================== \n")
      mod.sel <- NA
      sel.var.f <- NA
      vpa <- NA
    } else {
      sp.ls.mod.sel <- names(sp.ls.mod.sig) %>%
        map(function(nm){
          cat("=============================== \n")
          cat(paste("Forward selection procedure for: ",nm, "\n"))
          md.sel <- try(adespatial::forward.sel(sp.h,env_ls.f[[nm]],...),silent = T)
        })
      names(sp.ls.mod.sel) <- names(sp.ls.mod.sig)
      # extract selected variables
      sp.ls.mod.sel.var <- map(sp.ls.mod.sel,function(mod.sel){
        sel.var <- try(as.vector(mod.sel$variables),silent = T)
      })
      sp.ls.mod.sel.var.f <- sp.ls.mod.sel.var[!(map_lgl(sp.ls.mod.sel.var,function(x)class(x)=="try-error"))]

      # if no variables were selected
      if (length(sp.ls.mod.sel.var.f)==1 & is.null(sp.ls.mod.sel.var.f)){
        cat("=============================== \n")
        cat("Results of forward selections \n")
        cat("Warning: no variables are selected \n")
        cat("=============================== \n")
        mod.sel <- NA
        sel.var.f <- NA
        vpa <- NA
      } else {

        # final selected variables
        sel.var.f <- names(sp.ls.mod.sel.var.f) %>%
          map(function(nm){
            env_ls.f[[nm]] %>%
              dplyr::select(sp.ls.mod.sel.var.f[[nm]])
          })
        names(sel.var.f) <- names(sp.ls.mod.sel.var.f)
        if (length(sel.var.f)<2){
          cat("=============================== \n")
          cat("VPA model \n")
          cat("Warning: two to four explanatory variables were needed for varpart \n")
          cat("=============================== \n")
          mod.sel <- sp.ls.mod.sel
          sel.var.f <- sel.var.f
          vpa <- NA
        } else {

          vpa.expr <- paste0("varpart(","sp.h,",
                             paste0("env_ls.f[[",1:length(sel.var.f),"]]",collapse = ","),
                             ")")
          mod.sel <- sp.ls.mod.sel
          sel.var.f <- sel.var.f
          vpa <- eval(rlang::parse_expr(vpa.expr))}
      }
    }

    return(list(mode=sp.ls.mod,
                mod.aov=sp.ls.mod.aov,
                mod.sel=mod.sel,
                sel.var.f=sel.var.f,
                vpa=vpa))

  }

  # DCA will be conducted only with sp.based==F, i.e., for the whole community, not species based VPA
  if (sp.based==F){
    # DCA analysis
    sp.dca <- sp.h %>%  decorana()
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
      mod <- "cca"
    } else if (sp.dca.axleng>=3){
      if (is.null(mode)){
        mod <- "rda"
      } else {
        mod <- "cca"
      }
    } else {
      mod <- "rda"
    }
    #print the mode will be used
    # print the axis length
    cat("=============================== \n")
    cat("DCA analysis and model selection: \n")
    cat(paste("Max of axis length in DCA is: ",sp.dca.axleng %>% specify_decimal(2),"\n"))
    cat(paste("Mode used: ", mod, "\n"))
    cat("=============================== \n")
    mod.res <- con.mod.sel(sp.h = sp.h, env_ls.f = env_ls.f,mode=mod)

  } else {
    mod <-ifelse(is.null(mode),"rda",mode)
    sp.dca <- NA
    cat("=============================== \n")
    cat("Species based VPA: \n")
    cat(paste("DCA is not conducted: ","\n"))
    cat(paste("Mode used (use RDA by default): ", mod, "\n"))
    cat(paste("Transformation for species data (hellinger by default): ", sp.trans, "\n"))
    cat("=============================== \n")
    # construct sp list with each species as an element in a list
    sp.h.list <- map(1:ncol(sp.num),function(n){
      sp.n <- sp.num[,n,drop=F] %>%
        vegan::decostand(method = sp.trans)
    })
    names(sp.h.list) <- names(sp.num)
    n.length <- length(sp.h.list)
    mod.res <- map(1:n.length,function(n){
      nm <- names(sp.h.list)[n]
      cat("=============================== \n")
      cat(paste0("Species based VPA for ",nm, " (", n, " out of ", n.length,") \n"))
      sp.df <- sp.h.list[[n]]
      sp.df.mod <- con.mod.sel(sp.h = sp.df, env_ls.f = env_ls.f,mode=mod)
      return(sp.df.mod)
    })
    names(mod.res) <- names(sp.h.list)
  }

  return(list(data=list(sp=sp,env_ls=env_ls),
              mod.info=c(sp.trans=sp.trans,env.trans=env.trans,
                         mode.used=mod,sp.dca=list(sp.dca),
                         sp.based=sp.based),
              mod.res=mod.res))

  cat("==============END=============== \n")

}
