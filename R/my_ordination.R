#' Customized function for PCoA and NMDS
#' @param sp.df Species abundance data, row as sites while column as species
#' @param sp.dist WHich distance should be used, by default, bray distance will be used
#' @param method Either PCoA or NMDS should be conducted, by default, PCoA will be performed with ape::pcoa.
#' Since it is recommended that PCoA are suitable for most of the applications.
#' @param axes Indicating the numbers of axes interested, by default axes=2.
#' @param ... More arguments can be passed to vegan::metaMDS or ape::pcoa
#'
#' @return A list including 1) the ordination model and 2) the summary for further visualizations
#' @export
#'
#' @examples
my_ordination <-
  function(sp.df,
           sp.dist = NULL,
           method = NULL,
           axes = NULL,
           ...) {
    sp.dist <- ifelse(is.null(sp.dist), "bray", sp.dist)
    method <- ifelse(is.null(method), "PCoA", method)
    axes <- ifelse(is.null(axes), 2, axes)
    axes <- 1:axes
    sp.df.nonnum <- dplyr::select_if(sp.df, function(x)
      ! is.numeric(x))
    sp.df.num <- dplyr::select_if(sp.df, is.numeric)
    msg <- paste0(method, " will be conducted \n")
    message(msg)
    if (method == "PCoA") {
      sp.df.bray <- vegan::vegdist(sp.df.num, method = sp.dist)
      sp.ord <- ape::pcoa(sp.df.bray, ...)
      site_coord <- sp.ord$vectors %>%
        as.data.frame() %>%
        dplyr::select(axes)
      names(site_coord) <- paste("PCoA", axes, sep = "")

      eig <- sp.ord$values %>%
        as.data.frame()
      # check if correction is applied
      corr.info <- sp.ord$correction
      if (corr.info[1] != "none" & corr.info[2] != 1) {
        eig.df <- eig %>%
          as.data.frame() %>%
          dplyr::select("Corr_eig", "Cum_corr_eig") %>%
          rename(eig = "Corr_eig", cum_eig = "Cum_corr_eig")
      } else {
        eig.df <- eig %>%
          as.data.frame() %>%
          dplyr::select("Eigenvalues", "Cumul_eig") %>%
          rename(eig = "Eigenvalues", cum_eig = "Cumul_eig")
      }
      eig.df <- eig.df[axes, ] %>%
        unlist %>%
        matrix(nrow = 1, byrow = T) %>%
        as.data.frame()
      names(eig.df) <- paste(rep(names(site_coord), 2),
                             rep(c("eig", "cum_eig"), each = length(axes)), sep =
                               "_")
      # merge site_corrd and eig.df
      sp.ord.df <- cbind(site_coord, eig.df[rep(1, nrow(site_coord)), ])
      #sp.ord.df <- eig.df
    } else if (method == "NMDS") {
      sp.ord <- vegan::metaMDS(sp.df.num, distance = "bray", ...)
      site_coord <- sp.ord$points %>%
        as.data.frame
      # add stress valur to site_coord
      site_coord$stress <- sp.ord$stress
      sp.ord.df <- site_coord
    } else {
      stop("Method should either PCoA or NMDS")
    }

    ##############
    sp.ord.df <- cbind(sp.df.nonnum, sp.ord.df) %>%
      as.data.frame() %>%
      bind_cols(nm = row.names(sp.df), .)
    ord.res <- list(mod = sp.ord, mod.summ = sp.ord.df)
    return(ord.res)
  }
