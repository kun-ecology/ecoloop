#' Get multiple comparison letters from anova (HSD.test,package agricolae)
#'
#' @param resp.var a character regarding the numeric data
#' @param exp.var a character regarding the group factor
#' @param data a dataframe
#' @param p.adj methods for adjuted p value, see function aruguments methods in function p.adjust for more infomation
#'
#' @return
#' @export
#'
#' @examples
get_aov.letter <- function(resp.var,exp.var,data,p.adj=NULL){
  p.adj <- ifelse(is.null(p.adj),"none",p.adj)
  specify_decimal <- function(x, k=NULL){
    if (is.null(k)) k=3
    else k=k
    trimws(format(round(x, k), nsmall=k))
  }
  require(agricolae)
  df <- data[,c(resp.var,exp.var)]
  #df <- na.omit(df)
  df<- dplyr::mutate_at(df,exp.var,factor) %>% as.data.frame()
  aov.formula <- paste(resp.var,exp.var,sep="~") %>% as.formula
  names(aov.formula) <- resp.var
  #normality <- (shapiro.test(as.numeric(df$resp.var)))$p.value
  #normality <- (ks.test(as.numeric(df$resp.var)))$p.value
  #install.packages("nortest")
  #library(nortest)
  #normality <- nortest::ad.test(df[,1] %>% as.numeric)$p.value # why should use number as indice for subseting??
  normality <- (shapiro.test(as.numeric(df[,1])))$p.value

  if (normality>0.05){
    homogeneity <- (bartlett.test(aov.formula, df))$p.value
  } else {
    homogeneity <- (car::leveneTest(aov.formula, df))$'Pr(>F)'[1]
  }
  if (normality>=0.05 & homogeneity<=0.05){
    aov.res <- aov(aov.formula,df)
    hsd.res <- HSD.test(aov.res,exp.var)
    hsd.res <- hsd.res$groups %>% as.data.frame()
    hsd.res$Treatment <- row.names(hsd.res)
    aov.res.sum <- try(summary(aov.res),silent = T)
    hsd.res$anova.info <- paste0("ANOVA, F=",
                                 specify_decimal(aov.res.sum[[1]]$'F value'[1],2),
                                 ", P=",
                                 specify_decimal(aov.res.sum[[1]]$'Pr(>F)'[1]))
  } else {
    kruskal.res <- agricolae::kruskal(df[,resp.var],df[,exp.var],group=T,p.adj = p.adj) #df$resp.var not working?
    hsd.res <- kruskal.res$groups %>% as.data.frame()
    hsd.res$Treatment <- row.names(hsd.res)
    #aov.res.sum <- hsd.res$statistics %>% unlist %>% as.vector() # convert to numeric vector
    hsd.res$anova.info <- paste0("Kruskal-Wallis, Chi.sq=",
                                 specify_decimal(kruskal.res$statistics[1],2),
                                 ", P=",
                                 specify_decimal(kruskal.res$statistics[3],2))
  }
  mean.sd <- aggregate(df[,1],by=list(df[,2]),cal_mean.sd) %>% as.data.frame()
  colnames(mean.sd) <- c(exp.var,"mean.sd")
  colnames(hsd.res)[1:3] <- c("mean","aov_letter",exp.var)
  hsd.res <-  left_join(hsd.res,mean.sd,by=exp.var)
  hsd.res <- dplyr::select(hsd.res,4,3,1,5,2)
  return(hsd.res)
  #return(list(normality,homogeneity))
}
