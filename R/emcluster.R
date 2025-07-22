..densityEMCluster <- function(x, nclass,init=NULL) {
  if (isS4(init)) { # from dMixmod to EMCluster
    if (nclass==length(pi <- init@proportions)) {
      sigmalist  <- init@variance
      LTSigma <- t(sapply(sigmalist,function(v) v[upper.tri(v,diag=TRUE)]))
      emobj <- list(pi=init@proportions,
                    Mu=init@mean,
                    LTSigma=LTSigma)
    } else init <- NULL
  } else if (inherits(init,"list")) { # from dMclust to EMCluster
    if (nclass==length(pi <- init$pro)) {
      sigma  <- init$variance$sigma
      LTSigma <- t(apply(sigma,3L,function(v) v[upper.tri(v,diag=TRUE)]))
      emobj <- list(pi=init$pro,
                    Mu=t(init$mean),
                    LTSigma=LTSigma)
    } else init <- NULL
  }
  if (is.null(init)) {
    simple.init <- .get_wrap("simple.init", pack="EMCluster")
    emobj <- simple.init(x=x, nclass=nclass)
    shortemcluster <- .get_wrap("shortemcluster", pack="EMCluster")
    emobj <- shortemcluster(x=x, emobj=emobj)
  } 
  emcluster <- .get_wrap("emcluster", pack="EMCluster")
  emobj <- emcluster(x=x, emobj=emobj)
  emobj
}

.get_best_emcluster_by_IC <- function(cluObject, x, which=Infusion.getOption("criterion")) {
  # if (inherits(cluObject,"try-error")) return(cluObject) ## passes original error info rather than hiding it under another error
  if (length(cluObject)==1L) return(cluObject[[1L]])
  ICs <- numeric(length(cluObject))
  if (which=="BIC") {
    em.bic <- get("em.bic", asNamespace("EMCluster"))
    for (it in seq_along(cluObject)) {
      ICs[it] <- em.bic(x=x,emobj=cluObject[[it]])
    }
    return(cluObject[[which.min(ICs)]])
  } else {
    em.aic <- get("em.aic", asNamespace("EMCluster"))
    for (it in seq_along(cluObject)) {
      ICs[it] <- em.aic(x=x,emobj=cluObject[[it]])
    }
    return(cluObject[[which.min(ICs)]])
  }
}

.vec2cov <- function(v,p){
  covcorr <- diag(p)
  covcorr[upper.tri(covcorr,diag = TRUE)] <- v
  covcorr[lower.tri(covcorr)] <- t(covcorr)[lower.tri(covcorr)]
  covcorr
}

.post_process_emresu <- function(resu, blob) {
  if (inherits(resu,"try-error")) return(resu)
  VARs <- apply(resu$LTSigma, 1, .vec2cov, p=resu$p, simplify = FALSE)
  chk <- try(lapply(VARs,chol), silent=TRUE)
  if (inherits(chk,"try-error")) return(chk)
  VARs <- simplify2array(VARs)
  dimnames(VARs) <- list(blob$varNames, blob$varNames, NULL)
  MEANs <- t(resu$Mu)
  rownames(MEANs) <- blob$varNames
  resu <- list(G=resu$nclass, 
               modelName="VVV", # maybe d,n,data, call, loglik...
               parameters=list(
                 pro=resu$pi,
                 mean=MEANs,
                 variance=list(modelName="VVV",
                               d=length(blob$varNames),
                               G=resu$nclass,
                               sigma=VARs)
               ))
}


# Alternative to (mixmodCluster+selection by AIC), 
# returns a dMclust object as the structure is simple and methods exist to handle it
# Does this create a dependence on the mclust package? Apparently not:
# the tests run with mclust removed from the installation.
.densityEMCluster <- function(data,stat.obs,nbCluster,
                              # modelNames=Infusion.getOption("mclustModel"), 
                              init=NULL,
                              verbose=FALSE, plot=FALSE) {
  if ( ! length(nbCluster)) return(NULL)
  # stat.obs useful for boundaries and to handle degenerate distributions:
  blob <- .process_boundaries(data, stat.obs, boundaries=attr(stat.obs,"boundaries")) # original role probably irrelevant
  #     but blob is used 
  if ( ! nrow(blob$data)) { # ie if length(blob$goodrows)==0L in .process_boundaries()
    ## code has long been invalid, but not tested. See further comments on parallel code for .densityMixmod()
    # This still ooks as an ugly patch since Rmixmod is potentially used in this mClust based function.
    resu <- new(structure("dMixmod", package="Infusion")) # line individually OK; no 'data' nor 'nbCluster' slots
    resu@strategy <- NULL
    resu@varNames <- blob$varNames
    resu@statNames <- blob$statNames
    resu@simuls_activeBoundaries <- blob$simuls_activeBoundaries
    resu@Sobs_activeBoundaries <- blob$Sobs_activeBoundaries
    resu@freq <- blob$freq # 1/(2*(nrow(data+1)))
  } else {
    # data lines matching the 'boundaries' values have now been removed. There may still be NA's in the data (notably, lines of NA's). The doc explains
    # "the inference method ... must thus ignore all empirical distributions including NA/NaN/inf.": 
    data <- .check_data_post_boundaries(blob$data)
    ## 
    if ((nc <- ncol(data))>(nr <- nrow(data))) {
      stop(paste0("Clustering expected to fail on data with more columns (",nc,
                  " retained parameters and/or statistics)) than rows (",nr," retained points)."))
    }
    if (length(nbCluster)==1L) {
      resu <- ..densityEMCluster(x=data, nclass=nbCluster,init=init)
      resu <- .post_process_emresu(resu, blob=blob)
      if (inherits(resu,"try-error")) { 
        resu <- ..densityEMCluster(x=data, nclass=nbCluster,init=init)
        resu <- .post_process_emresu(resu, blob=blob)
      }
      if (inherits(resu,"try-error")) { # Then we try any smaller size.
        if (.Infusion.data$options$mixturing_errorfn()) {
          cat(
            cli::col_green(
              paste0("emcluster(nclass=", nbCluster, ") failed. Trying decrease of numbers of clusters...\n")))
        }
        for (nb in rev(seq_len(min(nbCluster)-1L))) {
          resu <- ..densityEMCluster(x=data, nclass=nb)
          resu <- .post_process_emresu(resu, blob=blob)
          if ( ! inherits(resu,"try-error")) break
        }
      } 
    } else {
      emrets <- vector("list",length(nbCluster))
      for (it in nbCluster) {
        resu <- ..densityEMCluster(x=data, nclass=it) # of class "emret"
        if ( ! inherits(resu,"try-error")) emrets[[it]] <- resu 
      }
      resu <- .get_best_emcluster_by_IC(x=data, emrets) 
      resu <- .post_process_emresu(resu, blob=blob) # converts to dMclust format. class is c("dMclust", "list")
    }
    if (inherits(resu,"try-error")) stop("All emcluster() attempts failed.")
    # LTSigma is 'lower triangle' but read row-first so it's upper triangle when read col-first...
    resu <- structure(resu,
                      statNames=blob$statNames,
                      simuls_activeBoundaries=blob$simuls_activeBoundaries,
                      Sobs_activeBoundaries=blob$Sobs_activeBoundaries,
                      freq=blob$freq)
    class(resu) <- c("dMclust", class(resu))
  }
  resu
}


