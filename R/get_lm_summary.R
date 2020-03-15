#' Get tidy results from linear, nls, lmer, lme,  sarlm mode
#'
#' @param lm.mod a model returned by function lm, nls, lmer, lme
#' @param CI whether return confidence interval for the estimates, only works for lm model
#' @param k numbers of digits
#'
#' @return
#' @export
#'
#' @examples
get_lm_summary <- function(lm.mod,CI=NULL,k=NULL){
  require(tidyverse)
  k <- ifelse(is.null(k),3,k)
  specify_decimal <- function(x){
    trimws(format(round(x, k) %>% try(silent = T), nsmall=k))
  }
  if (class(lm.mod)=="try-error"){
    lm.summary.res <- rep(NA,10)
    lm.summary.res <- t(lm.summary.res) %>% as.data.frame()
    colnames(lm.summary.res) <- c("Vars", "Estimate", "Std.Error", "t.value",
                                  "P", "R2", "AdjR2", "lm.P", "AIC", "BIC")
    return(lm.summary.res)
  } else if (class(lm.mod)=="lm") {

    lm.res <- summary(lm.mod)
    lm.summary.res <- lm.res$coefficients %>% as.data.frame()
    lm.summary.res$R2 <- lm.res$r.squared
    lm.summary.res$AdjR2<- lm.res$adj.r.squared
    lm.summary.res$lm.P <- pf(lm.res$fstatistic[1],
                              lm.res$fstatistic[2],
                              lm.res$fstatistic[3],
                              lower.tail=FALSE)
    lm.summary.res$AIC <- AIC(lm.mod)
    lm.summary.res$BIC <- BIC(lm.mod)
    lm.summary.res <- mutate(lm.summary.res,Vars=row.names(lm.summary.res))
    lm.summary.res <- mutate_if(lm.summary.res,is.numeric,specify_decimal)
    lm.summary.res <- mutate_at(lm.summary.res,1:8,function(x)suppressWarnings(as.numeric(x)))
    lm.summary.res <- lm.summary.res[,c(10,1:9)]
    names(lm.summary.res)[3:5] <- c("Std.Error","t.value","P")

  } else if (class(lm.mod)=="nls"){
    lm.res <- summary(lm.mod)
    lm.summary.res <- lm.res$coefficients %>% as.data.frame()
    lm.summary.res$R2 <- NA
    lm.summary.res$AdjR2<- NA
    lm.summary.res$lm.P <- NA
    lm.summary.res$AIC <- AIC(lm.mod)
    lm.summary.res$BIC <- BIC(lm.mod)
    lm.summary.res <- mutate(lm.summary.res,Vars=row.names(lm.summary.res))
    lm.summary.res <- mutate_if(lm.summary.res,is.numeric,specify_decimal)
    lm.summary.res <- mutate_at(lm.summary.res,1:8,function(x)suppressWarnings(as.numeric(x)))
    lm.summary.res <- lm.summary.res[,c(10,1:9)]
    names(lm.summary.res)[3:5] <- c("Std.Error","t.value","P")

  } else if (class(lm.mod)=="lme"){
    lm.res <- summary(lm.mod)
    lm.summary.res <- lm.res$tTable %>% as.data.frame()
    lm.summary.res <- lm.summary.res[,c(1,2,4,5)]
    lm.summary.res$R2 <- NA # do not know why can not get R2 with r.squaredGLMM function
    lm.summary.res$AdjR2<- NA
    lm.summary.res$lm.P <- NA
    lm.summary.res$AIC <- AIC(lm.mod)
    lm.summary.res$BIC <- BIC(lm.mod)
    lm.summary.res <- mutate(lm.summary.res,Vars=row.names(lm.summary.res))
    lm.summary.res <- mutate_if(lm.summary.res,is.numeric,specify_decimal)
    lm.summary.res <- mutate_at(lm.summary.res,1:8,function(x)suppressWarnings(as.numeric(x)))
    lm.summary.res <- lm.summary.res[,c(10,1:9)]
    names(lm.summary.res)[3:5] <- c("Std.Error","t.value","P")
  } else if (class(lm.mod)=="lmerMod"){
    lm.res <- summary(lm.mod)
    lm.summary.res <- lm.res$coefficients %>% as.data.frame()
    lm.summary.res$P <- NA
    lm.summary.res$P[2:nrow(lm.summary.res)] <- (car::Anova(lm.mod))[3] %>% unlist()
    lm.summary.res$R2 <- MuMIn::r.squaredGLMM(lm.mod)[2] #Marginal R_GLMM² represents the variance explained by the entire model
    lm.summary.res$AdjR2<- MuMIn::r.squaredGLMM(lm.mod)[1] #Marginal R_GLMM² represents the variance explained by the fixed effects
    lm.summary.res$lm.P <- NA # no such value for lmer package
    lm.summary.res$AIC <- AIC(lm.mod)
    lm.summary.res$BIC <- BIC(lm.mod)
    lm.summary.res <- mutate(lm.summary.res,Vars=row.names(lm.summary.res))
    lm.summary.res <- mutate_if(lm.summary.res,is.numeric,specify_decimal)
    lm.summary.res <- mutate_at(lm.summary.res,1:8,function(x)suppressWarnings(as.numeric(x)))
    lm.summary.res <- lm.summary.res[,c(10,1:9)]
    names(lm.summary.res)[3:5] <- c("Std.Error","t.value","P")
  } else if (class(lm.mod)=="sarlm") {

    lm.res <- summary(lm.mod, correlation=F, Nagelkerke=TRUE, Hausman=TRUE)
    lm.summary.res <- lm.res$Coef %>%
      as.data.frame() %>%
      bind_cols(Vars=row.names(.),.)
    lm.summary.res$R2 <- lm.res$NK
    lm.summary.res$AIC <- AIC(lm.mod)
    lm.summary.res$BIC <- BIC(lm.mod)
    lm.summary.res <- mutate_if(lm.summary.res,is.numeric,specify_decimal)
  }





  # only for linear model
  if (!is.null(CI)) {
    if (class(lm.mod)!="lme"){
      lm.ci <- confint(lm.mod) %>% as.data.frame() %>% mutate(Vars=row.names(.))
      lm.ci <- lm.ci[,c(3,1,2)]
      names(lm.ci)[2:3] <- c("lower_2.5","upper_97.5")
      lm.summary.res <- left_join(lm.summary.res,lm.ci,by="Vars")
      return(lm.summary.res)
    } else {
      lm.ci <- intervals(lm.mod)
      lm.ci <- as.data.frame(lm.ci$fixed)
      lm.ci <- dplyr::mutate(lm.ci,Vars=row.names(lm.ci))
      lm.summary.res <- left_join(lm.summary.res,lm.ci)
      return(lm.summary.res)
    }

  }

  return(lm.summary.res)



}

