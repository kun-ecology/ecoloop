#' Customed mantel
#'
#' @param y_mat a dist or dataframe
#' @param x_mat a dist or dataframe
#' @param y.dist methods for dissimilarity calculation of y_mat
#' @param x.dist methods for dissimilarity calculation of x_mat
#' @param foreach conduct mantel test for each variables in x.dist or not (by default)
#' @param method Correlation method, as accepted by cor: "pearson", "spearman" or "kendall"
#' @param ... further arguments for vegan::mantel
#'
#' @return a list contain 1) mode information. 2) pairwise dissimilarity of x_mat and y_mat, for visualization
#' 3) mode. 4) summary of the mode
#' @export
#'
#' @examples
my_mantel <- function(y_mat, x_mat,y.dist=NULL, x.dist=NULL,foreach=NULL,method=NULL,...){
  # intiate method
  method <- ifelse(is.null(method),"pearson",method)
  foreach <- ifelse(is.null(foreach),FALSE, TRUE)
  ######################
  # initiate y_mat
  if (class(y_mat)=="dist"){
    y_mat.f <- y_mat
    y.dist <- NA
  } else {
    if (is.null(y.dist)){
      stop("When y_mat is not a dissimilarity matrix, y.dist (method for dissimilarity) is needed \n")
    } else {
      y_mat.f <- vegan::vegdist(y_mat %>% dplyr::select_if(is.numeric), method = y.dist)
    }
  }
  #####################
  # initiate x_mat
  if (class(x_mat)=="dist"){
    x_mat.f <- x.mat
    x_mat.df <- ecoloop::dist2df(x_mat.f) # for saving 
    x.dist <- NA
    mantel.mod <- vegan::mantel(y_mat.f,x_mat.f,method = method,...)
    mod.res <- data.frame(Mantel_Cor=mantel.mod$statistic,
                          Mantel_p=mantel.mod$signif)
  } else {
    if (is.null(x.dist)){
      stop("When x_mat is not a dissimilarity matrix, x.dist (method for dissimilarity) is needed \n")
    } else {
      
      if (foreach==T){
        x_mat.num <- dplyr::select_if(x_mat,is.numeric)
        x_mat.f <- map(1:ncol(x_mat.num),function(n){
          x_mat.n <- x_mat.num[,n,drop=F] %>% 
            vegan::vegdist(method = x.dist)
        })
        names(x_mat.f) <- names(x_mat.num)
        x_mat.df0 <- map(x_mat.f,ecoloop::dist2df)
        x_mat.df1 <- do.call(rbind.data.frame,x_mat.df0)
        x_mat.df <-  dplyr::bind_cols(x_mat.nm=rep(names(x_mat.f),map_dbl(x_mat.df0,nrow)),
                                      x_mat.df1)
        
        ####################
        mantel.mod <- map(x_mat.f,function(dis){
          dis.mantel <- vegan::mantel(y_mat.f,dis,method = method,...)
        })
        ###################
        mod.res <- map(mantel.mod,function(mod){
          data.frame(Mantel_Cor=mod$statistic,
                     Mantel_p=mod$signif)
        }) %>% 
          do.call(rbind.data.frame,.) %>% 
          bind_cols(x_mat.nm=names(mantel.mod),.)
        
      } else {

        x_mat.f <- vegan::vegdist(x_mat %>% dplyr::select_if(is.numeric), method = x.dist)
        x_mat.df <- ecoloop::dist2df(x_mat.f)
        mantel.mod <- vegan::mantel(y_mat.f,x_mat.f,method = method,...)
        mod.res <- data.frame(Mantel_Cor=mantel.mod$statistic,
                              Mantel_p=mantel.mod$signif)
      }
      
    
    }
  }
  ######################
  # final output
  y_mat.df <- ecoloop::dist2df(y_mat.f)
  xy_mat.df <- dplyr::left_join(x_mat.df, y_mat.df,by=c("dist.x","dist.y")) %>% 
    dplyr::rename(Pair_x.nm="dist.x",Pair_y.nm="dist.y",x_mat.dist="value.x",y_mat.dist="value.y")

  res.f <- list(mode.info=list(y.dist=y.dist,x.dist=x.dist,method=method),
                xy_mat.df=xy_mat.df,
                mode=mantel.mod,
                mode.res=mod.res)
  
  
}
