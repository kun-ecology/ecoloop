#' Fast computation of MPD and SES.MPD
#'
#' @param comm a data.frame containing species abundance, better set site name as rowname
#' @param phy_dist a matrix containing phylogenetic distance
#' @param null_model only implemented for random taxa.labels
#' @param nworkers numbers of cores used for computation
#' @param ses if SES.MPD should be calculated
#' @param abundance_weighted whether MPD/SES.MPD is abundance weighted
#' @param runs usually 999 random draws, can be set to to smaller number to reduce computation time
#'
#' @return
#' @export
#'
#' @examples
fst.ses.mpd <- function(comm,
                        phy_dist,
                        null_model = "taxaShuffle",
                        nworkers = NULL,
                        ses = F,
                        abundance_weighted = F,
                        runs = NULL) {
  require(tidyverse)
  require(furrr)
  options(future.rng.onMisuse="ignore")

  # a function for shuffling taxa labels of the dist
  taxaShuffle <- function(x) {
    if (!is.matrix(x))
      x <- as.matrix(x)
    rand.names <- sample(rownames(x))
    rownames(x) <- rand.names
    colnames(x) <- rand.names
    return(x)
  }
  # a function for mpd
  fn1 <- function(nm, n){
    ntaxa <- length(nm)

    # if ntaxa=1, mpd is not availiable
    if(ntaxa<=1){
      mpd.obs <- NA
    } else {
      # for obs.mpd
      comm.tmp <- comm[n, nm]
      obs.dist <- phy_dist[nm, nm]
      if (abundance_weighted) {
        obs.w <- t(as.matrix(comm.tmp)) %*% as.matrix(comm.tmp)
        mpd.obs <- weighted.mean(obs.dist, obs.w)
      } else {
        mpd.obs <- mean(obs.dist[lower.tri(obs.dist)])
      }
    }

    data.frame(obs=mpd.obs)
  }

  # a function for calculating ses.mpd
  fn2 <- function(nm, n){
    ntaxa <- length(nm)

    # if ntaxa=1, mpd is not availiable
    if(ntaxa<=1){
      mpd.obs <- NA
      rand.mpd.mean <- NA
      rand.mpd.sd <- NA
      mpd.obs.z <- NA
      mpd.obs.rank <- NA
      mpd.obs.p <- NA
    } else {
      # for obs.mpd
      comm.tmp <- comm[n, nm]
      obs.dist <- phy_dist[nm, nm]
      obs.w <- t(as.matrix(comm.tmp)) %*% as.matrix(comm.tmp)

      ###########
      # for rand mpd
      shuffled.dist <- map(1:runs, ~ taxaShuffle(phy_dist))
      names(shuffled.dist) <- paste0("r", 1:runs)
      shuffled.dist <- map(shuffled.dist, ~ .x[nm, nm])

      #######
      if (abundance_weighted) {
        mpd.obs <- weighted.mean(obs.dist, obs.w)
        shuffled.mpd <-map_dbl(shuffled.dist,  ~ weighted.mean(.x , obs.w))
      } else {
        mpd.obs <- mean(obs.dist[lower.tri(obs.dist)])
        shuffled.mpd <- map_dbl(shuffled.dist, ~ mean(.x[lower.tri(.x)]))
      }

      rand.mpd.mean <- mean(shuffled.mpd)
      rand.mpd.sd <- sd(shuffled.mpd)
      mpd.obs.z <- (mpd.obs - rand.mpd.mean)/rand.mpd.sd
      mpd.obs.rank <- rank(c(mpd.obs, shuffled.mpd))[1]
      mpd.obs.p <- mpd.obs.rank/(runs+1)
    }

    #list(rand.mpd.mean, rand.mpd.sd, mpd.obs.z, mpd.obs.rank, mpd.obs.p)
    data.frame(
      obs = mpd.obs,
      rand.mean = rand.mpd.mean,
      rand.sd = rand.mpd.sd,
      obs.rank = mpd.obs.rank,
      obs.z = mpd.obs.z,
      obs.z.p = mpd.obs.p
    )

  }

  ##################
  # initialize arguments
  runs <- ifelse(is.null(runs), 999, runs)
  nworkers <- ifelse(is.null(nworkers), future::availableCores()/2, nworkers)


  ############
  # species names
  #comm1 <- comm[map_lgl(comm, is.numeric)]
  spe.nm <- names(comm)

  plan("multisession", workers=nworkers)
  # species name for each site
  site.spe <- apply(comm,1, function(x)spe.nm[x>0])
  #########
  # loop for each site
  if (ses) {
    res.ls <- future_imap(site.spe, fn2, .progress = T)
  } else {
    res.ls <- future_imap(site.spe, fn1, .progress = T)
  }
  plan("sequential")
  gc()

  res.df <- do.call(bind_rows, res.ls)

  # add row.names if available
  if (!is.null(row.names(comm))) {
    res.df <- bind_cols(site = row.names(comm), res.df)
  }

  res.df <- as_tibble(res.df)
  message("\n")
  return(res.df)
}


#######################

#'  fast computation of MNTD/SES.MNTD
#'
#' @param comm a data.frame containing species abundance
#' @param phy_dist a matrix containing phylogenetic distance
#' @param null_model only implemented for random taxa.labels
#' @param nworkers numbers of cores used for computation
#' @param ses if SES.MNTD should be calculated
#' @param abundance_weighted whether MNTD/SES.MNTD is abundance weighted
#' @param runs usually 999 random draws, can be set to to smaller number to reduce computation time
#'
#' @return
#' @export
#'
#' @examples
fst.ses.mntd <- function(comm,
                         phy_dist,
                         null_model="taxaShuffle",
                         nworkers = NULL,
                         ses = NULL,
                         abundance_weighted=F,
                         runs=NULL){
  require(tidyverse)
  require(furrr)
  ses <- ifelse(is.null(ses), F, T)

  options(future.rng.onMisuse="ignore")
  # a function for shuffling taxa labels of the dist
  taxaShuffle <- function(x) {
    if (!is.matrix(x))
      x <- as.matrix(x)
    rand.names <- sample(rownames(x))
    rownames(x) <- rand.names
    colnames(x) <- rand.names
    return(x)
  }

  fn1 <- function(nm, n){
    ntaxa <- length(nm)
    # if ntaxa=1, mpd is not availiable
    if(ntaxa<=1){
      mntd.obs <- NA
    } else {
      # for obs.mntd
      comm.tmp <- comm[n, nm]
      obs.dist <- phy_dist[nm, nm]
      diag(obs.dist) <- NA
      obs.w <- as.matrix(comm.tmp)
      #######
      if (abundance_weighted) {
        mntds <- apply(obs.dist, 2, min, na.rm=T)
        mntd.obs <- weighted.mean(mntds, obs.w)
      } else {
        mntds <- apply(obs.dist, 2, min, na.rm=T)
        mntd.obs <- mean(mntds)
      }
    }

    data.frame(obs=mntd.obs)
  }

  # a function for calculating ses.mntd
  fn2 <- function(nm, n){
    ntaxa <- length(nm)

    if(ntaxa<=1){
      mntd.obs <- NA
      rand.mntd.mean <- NA
      rand.mntd.sd <- NA
      mntd.obs.z <- NA
      mntd.obs.rank <- NA
      mntd.obs.p <- NA
    } else {
      # for obs.mntd
      comm.tmp <- comm[n, nm]
      obs.dist <- phy_dist[nm, nm]
      diag(obs.dist) <- NA
      obs.w <- as.matrix(comm.tmp)

      ###########
      # for rand mntd
      shuffled.dist <- map(1:runs, ~ taxaShuffle(phy_dist))
      names(shuffled.dist) <- paste0("r", 1:runs)
      shuffled.dist <- map(shuffled.dist, ~ .x[nm, nm])
      shuffled.dist <- map(shuffled.dist, function(x){diag(x) <- NA; x})
      # randomed mntds
      shuffled.mntds <- map(shuffled.dist, ~ apply(.x, 2, min, na.rm=T))

      #######
      if (abundance_weighted) {
        mntds <- apply(obs.dist, 2, min, na.rm=T)
        mntd.obs <- weighted.mean(mntds, obs.w)
        shuffled.mntd <-map_dbl(shuffled.mntds,  ~ weighted.mean(.x , obs.w))
      } else {
        mntds <- apply(obs.dist, 2, min, na.rm=T)
        mntd.obs <- mean(mntds)
        shuffled.mntd <-map_dbl(shuffled.mntds, mean)
      }

      rand.mntd.mean <- mean(shuffled.mntd)
      rand.mntd.sd <- sd(shuffled.mntd)
      mntd.obs.z <- (mntd.obs - rand.mntd.mean)/rand.mntd.sd
      mntd.obs.rank <- rank(c(mntd.obs, shuffled.mntd))[1]
      mntd.obs.p <- mntd.obs.rank/(runs+1)
    }

    #list(rand.mntd.mean, rand.mntd.sd, mntd.obs.z, mntd.obs.rank, mntd.obs.p)
    data.frame(
      obs = mntd.obs,
      rand.mean = rand.mntd.mean,
      rand.sd = rand.mntd.sd,
      obs.rank = mntd.obs.rank,
      obs.z = mntd.obs.z,
      obs.z.p = mntd.obs.p
    )
  }

  ##################
  # initialize arguments
  runs <- ifelse(is.null(runs), 999, runs)
  nworkers <- ifelse(is.null(nworkers), future::availableCores()/2, nworkers)


  ############
  # species names
  #comm1 <- comm[map_lgl(comm, is.numeric)]

  plan("multisession", workers=nworkers)
  # species name for each site
  #site.spe <- future_map(1:nrow(comm1), ~ spe.nm[comm1[.x, ]>0])

  spe.nm <- names(comm)
  site.spe <- apply(comm, 1, function(x)spe.nm[x>0])


  # loop for each site
  if (ses){
    res.ls <- future_imap(site.spe, fn2, .progress=T)
  } else {
    res.ls <- future_imap(site.spe, fn1, .progress=T)
  }

  plan("sequential")
  gc()

  res.df <- do.call(bind_rows, res.ls)
  # add row.names if availiable
  if (!is.null(row.names(comm))) {
    res.df <- bind_cols(site = row.names(comm), res.df)
  }

  res.df <- as_tibble(res.df)
  message("\n")
  return(res.df)

}



#'  fast computation of inter-community mean nearest taxon distance (bNTI)
#'
#' @param comm a data.frame containing species abundance
#' @param phy_dist a matrix containing phylogenetic distance
#' @param abundance_weighted whether bNTI is abundance weighted
#' @param exclude_conspecifics Should conspecific taxa in different communities be exclude from MNTD calculations? (default = FALSE)
#' @param ses if ses.bNTI or NTI be calculated
#' @param runs 999 random draws, can be set to to smaller number to reduce computation time
#' @param nworkers numbers of cores used for computation
#'
#' @return
#' @export
#'
#' @examples
fst.comdistnt <- function(comm,
                          phy_dist,
                          exclude_conspecifics=FALSE,
                          nworkers = NULL,
                          ses = F,
                          abundance_weighted=FALSE,
                          runs = NULL
){
  # make sure species names in comm and dist match to each other
  # make sure dist is of dist class
  # make sure site is row.names
  # dat <- match.comm.dist(comm, dis)
  require(furrr)
  require(tidyverse)
  options(future.rng.onMisuse="ignore")

  # a function for shuffling taxa labels of the dist
  taxaShuffle <- function(x) {
    if (!is.matrix(x))
      x <- as.matrix(x)
    rand.names <- sample(rownames(x))
    rownames(x) <- rand.names
    colnames(x) <- rand.names
    return(x)
  }

  # initialize arguments
  runs <- ifelse(is.null(runs), 999, runs)
  nworkers <- ifelse(is.null(nworkers), future::availableCores()/2, nworkers)

  # numbers of site
  n.site <- nrow(comm)
  spe.nm <- names(comm)

  # species occurred in each site
  # spe4each <- map(1:n.site, ~ comm[.x, ]>0)
  # spe4each <- map(spe4each, ~ spe.nm[.x])
  # names(spe4each) <- row.names(comm)

  spe4each <- apply(comm, 1, function(x)spe.nm[x>0])

  # pairwise species names
  site.pair <- as.data.frame(combn(names(spe4each), 2))

  # a function that calculate mean nearest taxon distance for two set of species names and
  # a given phylogenetic distance
  fn <- function(nm.ls, phy_dist, abundance_weighted, exclude_conspecifics){
    nm1 <- nm.ls[1]
    nm2 <- nm.ls[2]
    spe.nm1 <- spe4each[[nm1]]
    spe.nm2 <- spe4each[[nm2]]

    # at least 1 species is needed for each site
    if(length(spe.nm1)>=1 & length(spe.nm2)>=1){
      tmp.dis <- phy_dist[spe.nm1, spe.nm2, drop=FALSE]

      # should conspecific taxa be excluded in different communities
      if (exclude_conspecifics) {
        tmp.dis[sample.dis == 0] <- NA
      }

      sample1NT <- apply(tmp.dis, 1, min, na.rm = TRUE)
      sample1NT[sample1NT == Inf] <- NA
      sample2NT <- apply(tmp.dis, 2, min, na.rm = TRUE)
      sample2NT[sample2NT == Inf] <- NA

      # pb based or abun based
      if (abundance_weighted) {
        sample1.weights <- as.numeric(comm[nm1, spe.nm1])
        sample2.weights <- as.numeric(comm[nm2, spe.nm2])
        # delete NA in either community
        if (any(is.na(sample1NT))) {
          sample1NT <- sample1NT[!is.na(sample1NT)]
          sample1.weights <- sample1.weights[!is.na(sample1NT)]
          sample1.weights <- sample1.weights / sum(sample1.weights)
        }


        if (any(is.na(sample2NT))) {
          sample2NT <- sample2NT[!is.na(sample2NT)]
          sample2.weights <- sample2.weights[!is.na(sample2NT)]
          sample2.weights <- sample2.weights / sum(sample2.weights)
        }

        # calculate the mean
        sampleNT <- c(sample1NT, sample2NT)
        sample.weights <- c(sample1.weights, sample2.weights)
        mntd <-weighted.mean(sampleNT, sample.weights, na.rm = TRUE)

      } else {mntd <- mean(c(sample1NT,sample2NT), na.rm=TRUE)  }

    } else {mntd <- NA}
  }

  ######################
  plan("multisession", workers=nworkers)
  if (ses) {
    message("calculating observed mntd\n")
    mntd.obs <-future_map_dbl(site.pair, ~ fn(.x, phy_dist = phy_dist,
                                              abundance_weighted =abundance_weighted,
                                              exclude_conspecifics = exclude_conspecifics
    ), .progress = T)



    message("\ncalculating randomized mntd\n")

    rand.commntd <- future_map(site.pair, function(nm.ls) {
      shuffled.dist <- map(1:runs, ~ taxaShuffle(phy_dist))
      names(shuffled.dist) <- paste0("r", 1:runs)
      # shuffled mntd
      rand.ls <- map_dbl(shuffled.dist, ~ fn(nm.ls, .x,
                                             abundance_weighted =abundance_weighted,
                                             exclude_conspecifics = exclude_conspecifics
      ))
    },.progress = T)

    # calculating means
    message("\ncalculating mean and sd\n")
    rand.commntd.mean <- map_dbl(rand.commntd, mean, na.rm=T)
    rand.commntd.sd <- map_dbl(rand.commntd, sd, na.rm=T)

    message("\nmerging data")
    # merge data
    rand.commntd.df <- data.frame(obs = mntd.obs,
                                  rand.mean = rand.commntd.mean,
                                  rand.sd = rand.commntd.sd) %>%
      mutate(obs.z = (obs - rand.mean) / rand.sd)

    res <- bind_cols(data.frame(site1 = as.character(site.pair[1, ]),
                                site2 = as.character(site.pair[2, ])),
                     rand.commntd.df) %>%
      as_tibble()
  } else {
    commntd <-future_map_dbl(site.pair, ~ fn(.x, phy_dist = phy_dist,
                                             abundance_weighted =abundance_weighted,
                                             exclude_conspecifics = exclude_conspecifics
    ), .progress = T)
    res <- tibble (
      site1 = as.character(site.pair[1, ,]),
      site2 = as.character(site.pair[2, ]),
      obs = commntd)

  }

  plan("sequential")

  invisible(res)


}

#' A function to calculate MPD or MNTD between native and aliens species
#'
#' @param comm a data.frame containing species abundance
#' @param phy_dist a matrix containing phylogenetic distance
#' @param native_sp.ls a vector containing native species
#' @param invasive_sp.ls a vector containing native species
#' @param ses if ses.bNTI or bNTI be calculated
#' @param method MPD (default) or MNTD
#' @param abundance_weighted should be indcies be abundance-weighted (default FALSE)
#' @param nworkers numbers of cores used for computation
#' @param runs 999 random draws, can be set to to smaller number to reduce computation time
#'
#' @return
#' @export
#'
#' @examples
md2nat <- function(comm,
                   phy_dist,
                   native_sp.ls,
                   invasive_sp.ls,
                   ses = F,
                   method="mpd",
                   abundance_weighted=F,
                   nworkers=NULL,
                   runs = NULL
){
  require(furrr)

  # initialize arguments
  nworkers <- ifelse(is.null(nworkers), future::availableCores()/2, nworkers)
  runs <- ifelse(is.null(runs), 999, runs)

  # a function for shuffling taxa labels of the dist
  taxaShuffle <- function(x) {
    if (!is.matrix(x)) x <- as.matrix(x)
    rand.names <- sample(rownames(x))
    rownames(x) <- rand.names
    colnames(x) <- rand.names
    return(x)
  }

  # new phy_dist which only contains native_sp.ls and invasive_sp.ls
  all.sp <- c(native_sp.ls, invasive_sp.ls)
  phy_dist <- phy_dist[all.sp, all.sp]

  # numbers of site
  n.site <- nrow(comm)
  spe.nm <- names(comm)

  # nm = species names
  # n = rows
  # fn1 = function for mntd.obs
  fn1 <- function(nat_abun, nat_spe, inv_abun, inv_spe){
    if(length(nat_spe)==0 | length(inv_spe)==0){
      mntd.obs <- NA
    } else {
      # for obs.mntd
      obs.dist <- phy_dist[nat_spe, inv_spe, drop=F]
      dist.min <- apply(obs.dist, 2, function(x)x==min(x, na.rm = T))
      mntds <- obs.dist[dist.min]
      ######
      if (abundance_weighted) {
        obs.w <- t(matrix(nat_abun, nrow = 1)) %*% matrix(inv_abun, nrow=1)
        obs.w <- obs.w[dist.min]
        mntd.obs <- weighted.mean(mntds, obs.w)
      } else {
        mntd.obs <- mean(mntds)
      }
    }
    data.frame(obs=mntd.obs)
  }

  # a function for calculating ses.mntd
  fn2 <- function(nat_abun, nat_spe, inv_abun, inv_spe){

    if(length(nat_spe)==0 | length(inv_spe)==0){
      mntd.obs <- NA
      rand.mntd.mean <- NA
      rand.mntd.sd <- NA
      mntd.obs.z <- NA
      mntd.obs.rank <- NA
      mntd.obs.p <- NA
      shuffled.dist <- NA
    } else {
      # for obs.mntd
      obs.dist <- phy_dist[nat_spe, inv_spe, drop=F]
      dist.min <- apply(obs.dist, 2, function(x)x==min(x))
      mntds <- obs.dist[dist.min]
      obs.w.tmp <- t(matrix(nat_abun, nrow = 1)) %*% matrix(inv_abun, nrow=1)
      obs.w <- obs.w.tmp[dist.min]

      ###########
      # for rand mntd
      shuffled.dist <- map(1:runs, ~ taxaShuffle(phy_dist))
      names(shuffled.dist) <- paste0("r", 1:runs)
      shuffled.dist <- map(shuffled.dist, ~ .x[nat_spe, inv_spe, drop=F])
      #randomed mntds
      shuffled.dist.min <- map(shuffled.dist, ~ apply(.x, 2, function(x)x==min(x)))
      shuffled.mntds <- map2(shuffled.dist, shuffled.dist.min, ~ .x[.y])
      shuffled.w <- map2(replicate(runs, obs.w.tmp, simplify = F),
                         shuffled.dist.min, ~ .x[.y])

      #######
      if (abundance_weighted) {

        mntd.obs <- weighted.mean(mntds, obs.w)
        shuffled.mntd <-map2_dbl(shuffled.mntds, shuffled.w,  ~ weighted.mean(.x , .y, na.rm=T))
      } else {
        mntd.obs <- mean(mntds)
        shuffled.mntd <-map_dbl(shuffled.mntds, mean, na.rm=T)
      }

      rand.mntd.mean <- mean(shuffled.mntd)
      rand.mntd.sd <- sd(shuffled.mntd)
      mntd.obs.z <- (mntd.obs - rand.mntd.mean)/rand.mntd.sd
      mntd.obs.rank <- rank(c(mntd.obs, shuffled.mntd))[1]
      mntd.obs.p <- mntd.obs.rank/(runs+1)
    }

    data.frame(obs=mntd.obs, rand.mean=rand.mntd.mean,
               rand.sd=rand.mntd.sd, obs.rank=mntd.obs.rank,
               obs.z=mntd.obs.z, obs.z.p=mntd.obs.p )
    #shuffled.dist
  }

  # function for MPD.obs
  fn3 <- function(nat_abun, nat_spe, inv_abun, inv_spe){

    if(length(nat_spe)==0 | length(inv_spe)==0){
      mpd.obs <- NA
    } else {
      # for obs.mpd
      obs.dist <- phy_dist[nat_spe, inv_spe, drop=F]
      if (abundance_weighted) {
        obs.w <- t(matrix(nat_abun, nrow = 1)) %*% matrix(inv_abun, nrow=1)
        mpd.obs <- weighted.mean(obs.dist, obs.w)
      } else {
        mpd.obs <- mean(obs.dist)
      }
    }

    data.frame(obs=mpd.obs)
  }

  # a function for calculating ses.mpd
  fn4 <- function(nat_abun, nat_spe, inv_abun, inv_spe){

    if(length(nat_spe)==0 | length(inv_spe)==0){
      mpd.obs <- NA
      rand.mpd.mean <- NA
      rand.mpd.sd <- NA
      mpd.obs.z <- NA
      mpd.obs.rank <- NA
      mpd.obs.p <- NA
    } else {
      # for obs.mpd
      obs.dist <- phy_dist[nat_spe, inv_spe, drop=F]


      ###########
      # for rand mpd
      shuffled.dist <- map(1:runs, ~ taxaShuffle(phy_dist))
      names(shuffled.dist) <- paste0("r", 1:runs)
      shuffled.dist <- map(shuffled.dist, ~ .x[nat_spe, inv_spe, drop=F])

      #######
      if (abundance_weighted) {
        obs.w <- t(matrix(nat_abun, nrow = 1)) %*% matrix(inv_abun, nrow=1)
        mpd.obs <- weighted.mean(obs.dist, obs.w)
        shuffled.mpd <-map_dbl(shuffled.dist,  ~ weighted.mean(.x , obs.w))
      } else {
        mpd.obs <- mean(obs.dist)
        shuffled.mpd <- map_dbl(shuffled.dist, mean)
      }

      rand.mpd.mean <- mean(shuffled.mpd)
      rand.mpd.sd <- sd(shuffled.mpd)
      mpd.obs.z <- (mpd.obs - rand.mpd.mean)/rand.mpd.sd
      mpd.obs.rank <- rank(c(mpd.obs, shuffled.mpd))[1]
      mpd.obs.p <- mpd.obs.rank/(runs+1)
    }

    data.frame(obs=mpd.obs, rand.mean=rand.mpd.mean,
               rand.sd=rand.mpd.sd, obs.rank=mpd.obs.rank,
               obs.z=mpd.obs.z, obs.z.p=mpd.obs.p )
  }

  plan("multisession", workers=nworkers)
  # species occurred in each site
  message("\nextracting invasive/native species for each community\n")

  # native species in each site
  nat_comm <- comm[,native_sp.ls, drop=F]
  nat_abun.ls <- apply(nat_comm, 1, function(x)x[x>0])
  nat_spe.ls <- apply(nat_comm, 1, function(x)native_sp.ls[x>0])

  # invasive species in each site
  inv_comm <- comm[,invasive_sp.ls, drop=F]
  inv_abun.ls <- apply(inv_comm, 1, function(x)x[x>0])
  inv_spe.ls <- apply(inv_comm, 1, function(x)invasive_sp.ls[x>0])

  df.tmp <- tibble(
    nat_abun.ls = apply(nat_comm, 1, function(x)x[x>0]),
    nat_spe.ls = apply(nat_comm, 1, function(x)native_sp.ls[x>0]),
    inv_abun.ls = apply(inv_comm, 1, function(x)x[x>0]),
    inv_spe.ls = apply(inv_comm, 1, function(x)invasive_sp.ls[x>0])
  )

  message("\ncalculating md2nat for each community\n")
  if (method=="mpd"){
    if (ses){
      res.tmp <- future_pmap(df.tmp, ~ fn4(..1, ..2, ..3, ..4), .progress = T)
    } else {
      res.tmp <- future_pmap(df.tmp, ~ fn3(..1, ..2, ..3, ..4), .progress = T)
    }
  } else {
    if (ses){
      res.tmp <- future_pmap(df.tmp, ~ fn2(..1, ..2, ..3, ..4), .progress = T)
    } else {
      res.tmp <- future_pmap(df.tmp, ~ fn1(..1, ..2, ..3, ..4), .progress = T)
    }
  }

  plan("sequential")

  gc()


  res.df <- res.tmp %>%
    do.call(rbind.data.frame,.) %>%
    bind_cols(site=row.names(comm),
              native.ntaxa=map_int(nat_spe.ls, length),
              invasive.ntaxa=map_int(inv_spe.ls, length),.) %>%
    as_tibble()

  return(res.df)

}
####################


#' fst.ses.pd a function for calculating PD
#'
#' @param comm a matrix contains site * species information, note that it must be a matrix
#' (will be automatically converted to sparse matrix) or sparse matrix.
#' @param phy a phylogenetic (or generated from functional distance) tree
#' @param null_model only implemented for taxaShuffle
#' @param nworkers numbers of cores for calculating, by default it only 1 core will be used.
#' @param ses a Boolean value indicates whether standardized form of PD will be calculated
#' @param runs an integer indicates numbers of randomization, by default 999, can be set to smaller values for tesing.
#'
#' @return a tibble contains observed PD and standardized PD (when ses = T)
#' @export
#'
#' @examples
fst.ses.pd <- function(comm,
                       phy,
                       null_model = "taxaShuffle",
                       nworkers = NULL,
                       ses = NULL,
                       runs = NULL) {
  # check if comm is a sparse Matrix
  if ( is(comm, "matrix") & !is(comm, "sparseMatrix")) {
    comm <- as(as.matrix(comm), "sparseMatrix")
    message("community matrix has been converted to sparseMatrix")
  } else if (!is(comm, "matrix") & !is(comm, "sparseMatrix") ) {
    stop("community matrix needs to be a (sparse) matrix!")
  }

  # check the input
  if (length(setdiff(colnames(comm), phy$tip.label)) > 0) {
    stop("There are species labels in community matrix missing in the tree!")
  }

  if (length(setdiff(phy$tip.label, colnames(comm))) > 0) {
    phy <- keep.tip(phy, intersect(phy$tip.label, colnames(comm)))
  }
  comm <- comm[, intersect(phy$tip.label, colnames(comm))]

  #########################
  # a function from phyloregion
  get_pd <-  function (comm, phy) {
    require(Matrix)
    require(phangorn)
    require(ape)
    el <- numeric(max(phy$edge))
    el[phy$edge[, 2]] <- phy$edge.length
    comm <- comm[, phy$tip.label]
    anc <- Ancestors(phy, seq_along(phy$tip.label))
    anc <- mapply(c, seq_along(phy$tip.label), anc, SIMPLIFY = FALSE)
    M <- sparseMatrix(as.integer(rep(seq_along(anc),lengths(anc))),
                      as.integer(unlist(anc)), x = 1L)
    commphylo <- comm %*% M
    commphylo@x[commphylo@x > 1e-08] <- 1
    list(Matrix = commphylo, edge.length = el)
    pd <- (commphylo %*% el)[, 1]
  }


  # a function for shuffle phylogenetic tree (only implemented for tipshuffle)
  tipshuffle <- function(phy) {
    phy$tip.label <- phy$tip.label[sample(length(phy$tip.label))]
    return(phy)
  }
  #########################
  # initiate arguments
  nworkers <- ifelse(is.null(nworkers), 1, nworkers)
  runs <- ifelse(is.null(runs), 999, runs)

  ######################################
  # calculate observed PD
  if (ses) {
    pd.obs <- get_pd(comm, phy)

    plan("multisession", workers = nworkers)
    pd.rand <- tibble(rep = paste0("rep", 1:runs)) %>%
      mutate(phy.rand = future_map(rep, ~ tipshuffle(phy), .progress = T)) %>%
      mutate(pd = future_map(phy.rand, ~ get_pd(comm, .x), .progress = T) ) %>%
      select(-phy.rand) %>%
      pivot_wider(names_from = rep, values_from = pd) %>%
      unnest(cols = 1:ncol(.))

    pd.rand.mean <- rowMeans(pd.rand, na.rm = T)
    pd.rand.sd <- apply(pd.rand, 1, function(x) sd(x, na.rm = T))
    pd.obs.rank <- apply(bind_cols(pd.obs = pd.obs, pd.rand), 1, function(x)rank(x)[1])

    plan("sequential")
    gc()
    ###########
    pd.obs.z <- (pd.obs - pd.rand.mean) / pd.rand.sd
    pd.obs.p <- pd.obs.rank / (runs + 1)

    res <-
      tibble(
        site = row.names(comm),
        obs = pd.obs,
        rand.mean = pd.rand.mean,
        rand.sd = pd.rand.sd,
        obs.rank = pd.obs.rank,
        obs.z = pd.obs.z,
        obs.z.p = pd.obs.p
      )
  } else {
    pd <- get_pd(comm, phy)
    res <- tibble(site = row.names(comm), obs = pd)
  }

  return(res)
}
