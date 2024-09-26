.check_data_post_boundaries <- function(data) {
  data <- na.omit(data)
  if (nNAlines <- length(attr(data,"na.action"))) {
    message(paste0(nNAlines," lines of reference table still contained NA's after handling 'boundaries', and were removed. see ?handling_NAs for details."))
  }
  checkDegenerate <- cov(data)
  #if ( ! anyNA(checkDegenerate)) { # cf  use="na.or.complete" => NA if no complete cases
  eigvals <- eigen(checkDegenerate, only.values = TRUE)$values
  if (any(abs(eigvals<1e-14))) warning(
    paste("The covariance matrix of the (parameters,statistics) table seems singular,\n",
          "implying vanishing variances, or linear relationships between the variables. Problems will likely happen.\n",
          "Redundant variables, if any, should be eliminated."
    ),immediate. = TRUE)
  #}
  data
}

.log2pi <- 1.837877066409345 # 

## fast version of dmvnorm that uses precomputed version of chol(sigma) and minimal checking
.fast_dmvnorm <- function (x, mean, solve_t_chol_sigma, log = FALSE) {
  if ( ! is.matrix(mean)) { # originally written for single mean vector and possibly multiple x vectors
    if ( ! is.matrix(x)) { # single x vector
      tmp <- solve_t_chol_sigma %*% (x-mean)
      rss <- sum(tmp*tmp)
    } else { # multiple x vectors
      tmp <- solve_t_chol_sigma %*% mean
      dim(tmp) <- NULL
      tmp <- tcrossprod(solve_t_chol_sigma, x) - tmp # solve_t_chol_sigma %*% (t(x)-mean) # backsolve(chol_sigma, t(x) - mean, transpose = TRUE)
      rss <- matrixStats::colSums2(tmp*tmp) # .colSums(tmp^2,m=nrow(tmp), n=ncol(tmp))
    }
  } else { # multiple mean vectors
    if ( ! is.matrix(x)) { # single x vector
      tmp <- solve_t_chol_sigma %*% (x-mean) # substract x from each col of 'mean'
    } else { # multiple x vectors or 1-row matrix
      if (nrow(x)==1L) x <- x[1,]
      tmp <- solve_t_chol_sigma %*% (x-mean) # backsolve(chol_sigma, t(x) - mean, transpose = TRUE)
    }
    rss <- .colSums(tmp*tmp,m=nrow(tmp), n=ncol(tmp))
  }
  logdet <- attr(solve_t_chol_sigma,"logdet")
  if (is.null(logdet)) logdet <- sum(log(diag(solve_t_chol_sigma))) # in (2017) method, precomputation is not interesting 
  logretval <- logdet - 0.5 * ncol(solve_t_chol_sigma) * .log2pi - 0.5 * rss
  # names(logretval) <- rownames(x)
  if (log) { logretval } else exp(logretval)
}

.marginalize_Rmixmod <- function(jointdens, 
                                 colNames, # names of all dims of the 'jointdens'; bc it is not in the object 
                                 For # retained params
                                 ) {
  nbCluster <- jointdens@nbCluster
  margdens <- jointdens
  MEAN <- margdens@parameters@mean
  colnames(MEAN) <- colNames
  margdens@parameters@mean <- MEAN[,For,drop=FALSE]
  for (clu_it in seq_len(nbCluster)) {
    COV <- margdens@parameters@variance[[clu_it]]
    colnames(COV) <- rownames(COV) <- colNames
    margdens@parameters@variance[[clu_it]] <- COV[For,For,drop=FALSE]
  }
  margdens@varNames <- For
  margdens@simuls_activeBoundaries <- NULL
  #margdens_fP@criterionValue <- margdens_fP@likelihood <- "removed for safety"  
  return(margdens)
}

.reduce_dMixmod <- function(dMixmod_obj, thr=1e-14) {
  parameters <- dMixmod_obj@parameters
  good <- (parameters@proportions>thr)
  parameters@mean <- parameters@mean[good,,drop=FALSE]
  parameters@variance <- parameters@variance[good]
  parameters@proportions <- parameters@proportions[good]
  attr(dMixmod_obj,"nbCluster") <- sum(good)
  dMixmod_obj@parameters <- parameters
  dMixmod_obj
}

# https://stats.stackexchange.com/questions/348941 if in doubt about the meaning
.conditional_Rmixmod <- function(jointdens, #fittedPars, 
                                 given, 
                                 # expansion is inflation factor for cov matrix:
                                 # expansion=1 to get the conditional distribution.
                                 expansion=Infusion.getOption("expansion")) {  
  nbCluster <- jointdens@nbCluster
  conddens <- jointdens
  MEAN <- conddens@parameters@mean
  givenNames <- names(given)
  colnames(MEAN) <- colNames <- jointdens@varNames
  For <- setdiff(colNames,givenNames) 
  conddens@parameters@mean <- MEAN[,For,drop=FALSE] # resizes, but will be modified
  condprop <- conddens@parameters@proportions
  for (clu_it in seq_len(nbCluster)) {
    # In the ricker example, a number of simulations are extinctions of the population. This creates a proba mass
    # of (projected) summstats, and the sig22 block is the zero matrix => regularize has no effect.
    # The (log) condprop 'should' then be (log)(0). So we first compute it.
    # Actually it is not certain that condprop will be zero, so check is not perfect (___F I X M E___) although it is effective. 
    COV <- conddens@parameters@variance[[clu_it]]
    colnames(COV) <- rownames(COV) <- colNames
    mean2 <- MEAN[clu_it,givenNames]
    sig22 <-  COV[givenNames,givenNames,drop=FALSE]
    condprop[clu_it] <- log(condprop[clu_it])+dmvnorm(t(given), # dmvnorm() tests is.vector(x) which returns FALSE if x has attributes other than names.
                                                      mean2, sigma= sig22, log=TRUE)
    if (condprop[clu_it]> -46) { # approx log(1e-20)
      rhs <- try(solve(sig22,given-mean2), silent=TRUE)
      if (inherits(rhs,"try-error")) {
        sig22  <- regularize(sig22)
        rhs <- solve(sig22,given-mean2)
      }
      sig12 <-  COV[For,givenNames,drop=FALSE]
      conddens@parameters@mean[clu_it,] <- MEAN[clu_it,For] + sig12 %*% rhs
      sig11 <- COV[For,For,drop=FALSE]
      conddens@parameters@variance[[clu_it]] <- expansion* (sig11 - sig12 %*% solve(sig22,t(sig12))) 
    } else { # probabily irrelevant, but let us keep a clean structure:
      conddens@parameters@mean[clu_it,] <- MEAN[clu_it,For]
      conddens@parameters@variance[[clu_it]] <- COV[For,For,drop=FALSE] 
    }
  }
  maxlog <- max(condprop)
  condprop <- condprop-maxlog
  condprop <- exp(condprop)
  condprop <- condprop/sum(condprop)
  conddens@parameters@proportions <- condprop
  conddens@varNames <- For
  conddens@simuls_activeBoundaries <- NULL
  #margdens_fP@criterionValue <- margdens_fP@likelihood <- "removed for safety"  
  
  # some cluster have vanishingly low conditional proba (... and possibly variance, which makes a good reason for removing them)
  conddens <- .reduce_dMixmod(conddens) 
  
  return(conddens)
}

.assign_Rmixmod <- function(jointdens, given) {  
  nbCluster <- jointdens@nbCluster
  conddens <- jointdens
  MEAN <- conddens@parameters@mean
  givenNames <- names(given)
  colnames(MEAN) <- colNames <- jointdens@varNames
  For <- setdiff(colNames,givenNames) 
  condprop <- prop <- conddens@parameters@proportions
  for (clu_it in seq_len(nbCluster)) {
    # In the ricker example, a number of simulations are extinctions of the population. This creates a proba mass
    # of (projected) summstats, and the sig22 block is the zero matrix => regularize has no effect.
    # The (log) condprop 'should' then be (log)(0). So we first compute it.
    # Actually it is not certain that condprop will be zero, so chekc is not perfect (___F I X M E___) although it is effective. 
    COV <- conddens@parameters@variance[[clu_it]]
    colnames(COV) <- rownames(COV) <- colNames
    mean2 <- MEAN[clu_it,givenNames]
    sig22 <-  COV[givenNames,givenNames,drop=FALSE]
    condprop[clu_it] <- log(condprop[clu_it])+dmvnorm(t(given), # dmvnorm() tests is.vector(x) which returns FALSE if x has attributes other than names.
                                                      mean2, sigma= sig22, log=TRUE)
  }
  maxlog <- max(condprop)
  condprop <- condprop-maxlog
  condprop <- exp(condprop)
  condprop <- condprop/sum(condprop)
  clu_id <- which.max(condprop)
  list(clu_id=clu_id,prop=prop[clu_id],condprop=condprop[clu_id])
}

.is_clustering_suspect <- function(slik_j) {
  if (inherits(slik_j$jointdens,"dMixmod")) {
    asgn <- .assign_Rmixmod(slik_j$jointdens, c(slik_j$MSL$MSLE, get_from(slik_j,"proj_data")))
    asgn$prop<0.01 && asgn$condprop>0.85 # 0.8989995 on a bad replicate...
  } else return(FALSE) # not checked rather than FALSE....
}


.conditional_mclust <- function(jointdens, given,
                                expansion=Infusion.getOption("expansion")) { # expansion=1 to get the conditional distribution
  nbCluster <- jointdens$G
  conddens <- jointdens
  MEAN <- conddens$parameters$mean 
  varNames <- rownames(MEAN)
  givenNames <- names(given)
  For <- setdiff(varNames, givenNames)
  conddens$data <- conddens$data[,For]
  conddens$parameters$variance$d <- length(For)
  conddens$parameters$mean <- MEAN[For,,drop=FALSE] # resizes, but will be modified
  sigma11 <- conddens$parameters$variance$sigma[For,For,,drop=FALSE] # resizes, but will be modified
  condprop <- conddens$parameters$pro
  for (clu_it in seq_len(nbCluster)) {
    sigma_it <- conddens$parameters$variance$sigma[,,clu_it] # from single array for all clusters # drops the clu_it dimension
    sig22 <-  sigma_it[givenNames,givenNames, drop=FALSE]
    sig12 <-  sigma_it[For,givenNames, drop=FALSE]
    mean2 <- MEAN[givenNames,clu_it]
    rhs <- try(solve(sig22,given-mean2), silent=TRUE)
    if (inherits(rhs,"try-error")) {
      sig22  <- regularize(sig22)
      rhs <- solve(sig22,given-mean2)
    }
    conddens$parameters$mean[,clu_it] <- MEAN[For,clu_it] + sig12 %*% rhs
    sigma11[,,clu_it] <- expansion* (sigma_it[For,For, drop=FALSE] - sig12 %*% solve(sig22,t(sig12))) 
    condprop[clu_it] <- condprop[clu_it]*dmvnorm(t(given), # dmvnorm() tests is.vector(x) which returns FALSE if x has attributes other than names.
                                                 mean2, sigma= sig22, log=FALSE)
  }
  maxlog <- max(condprop)
  condprop <- condprop-maxlog
  condprop <- exp(condprop)
  condprop <- condprop/sum(condprop) # normalisation (even in non-safe-log-exp version) was missing for a long time!
  conddens$parameters$pro <- condprop
  if (length(For)==1L) {
    if (length(sigma11)==1L) {
      conddens$parameters$variance <- list(modelName="X", d=1, G=1L, sigma=sigma11)
    } else {
      # "V" may not be implied by the original fit, but we don't count dfs on this object
      conddens$parameters$variance <- list(modelName="V", d=1, G=length(sigma11), 
                                           sigma=sigma11, scale=sigma11) # was sigmasq previously to v2.1.129. also in sigma2decomp call. Something to check...
    }
  } else {
    conddens$parameters$variance <- .do_call_wrap("sigma2decomp", 
                                                  list(sigma=sigma11,tol=-1), # tol to avoid a simplification in the structure of $orientation that may hit a but in mclust...
                                                  pack="mclust") 
  }
  if (conddens$modelName=="VVV") {
    cholsigma <- array(0,dim=dim(sigma11))
    for (it in seq_len(dim(cholsigma)[3L])) cholsigma[,,it] <- chol(sigma11[,,it]) # ___F I X M E____ unsafe chol()
    conddens$parameters$variance$cholsigma <- cholsigma
  }
  return(conddens)
}


.marginalize_mclust <- function(jointdens, colNames, For, over) {
  nbCluster <- jointdens$G
  margdens <- jointdens
  margdens$data <- margdens$data[,For, drop=FALSE]
  margdens$parameters$variance$d <- length(For)
  margdens$parameters$mean <- margdens$parameters$mean[For,,drop=FALSE]
  COV <- margdens$parameters$variance$sigma # single array for all clusters 
  COV <- COV[For,For, ,drop=FALSE]
  if (length(For)==1L) {
    if (length(COV)==1L) {
      variance <- list(modelName="X", d=1, G=1L, sigma=COV)
    } else {
      # "V" may not be implied by the original fit, but we don't count dfs on this object
      variance <- list(modelName="V", d=1, G=length(COV), sigma=COV, scale=COV) 
    }
  } else {
    variance <- .do_call_wrap("sigma2decomp", 
                                                  list(sigma=COV,tol=-1), # tol to avoid a simplification in the structure of $orientation that may hit a but in mclust...
                                                  pack="mclust") 
    # => this loses the names on the sigma's...
    dimnames( variance$sigma) <- list(For, For, NULL)
    if (margdens$modelName=="VVV") {
      cholsigma <- array(0,dim=dim(COV))
      for (it in seq_len(dim(cholsigma)[3L])) cholsigma[,,it] <- chol(COV[,,it])
      variance$cholsigma <- cholsigma
    }
  }
  margdens$parameters$variance <- variance
  # Not necess useful here, but following thesamelogic as in .conditional_mclust():
  margdens$modelName <- variance$modelName
  return(margdens)
}

.get_best_mixmod_by_IC <- function(cluObject, which=Infusion.getOption("criterion"))  {
  if (inherits(cluObject,"try-error")) return(cluObject) ## passes original error info rather than hiding it under another error
  results <- cluObject@results # assume that an nbCluster range was tried
  if (length(results)==1L) return(results[[1L]])
  anyNaN <- chkpartition <- logical(length(results))
  # patch for odd bug of Rmixmod: NaN in parameters but no error reported, and lik & BIC are real.
  for (it in seq_along(results)) {
    anyNaN[it] <- anyNA(results[[it]]@proba)
    chkpartition[it] <- min(table(results[[it]]@partition))<2L
  }
  results <- results[ ! (anyNaN | chkpartition)]
  BICs <- logLs <- numeric(length(results))
  for (it in seq_along(results)) {
    BICs[it] <- results[[it]]@criterionValue
    logLs[it] <- results[[it]]@likelihood
  }
  if (which=="BIC") return(results[[which.min(BICs)]])
  # ELSE
  dfs <- (2*logLs+BICs)/(log(cluObject@nbSample))
  AICs <- -2*logLs+2*dfs
  return(results[[which.min(AICs)]])
}

# .get_best_clu_by_AIC <- function(cluObject) { ## cluObject: mixmodCluster object
#   if (inherits(cluObject,"try-error")) return(cluObject) ## passes original error info rather than hiding it under another error
#   BICs <- unlist(lapply(cluObject@results,slot,name="criterionValue"))
#   logLs <- unlist(lapply(cluObject@results,slot,name="likelihood"))
#   dfs <- (2*logLs+BICs)/(log(cluObject@nbSample))
#   AICs <- -2*logLs+2*dfs
#   return(cluObject@results[[which.min(AICs)]])
# }

.solve_t_cholfn <- function(mat, smoothing_mat=NULL, condnum=1e12) { # smoothing_mat arg added to allow optimization via iterative smmothing (which does not work)
  if (!is.null(smoothing_mat)) mat <- mat+smoothing_mat
  resu <- try(solve(t(chol(mat))), silent=TRUE)
  if (inherits(resu,"try-error")) {
    # met this when sampling many points on a facet of the parameter space
    .Infusion.data$options$chol_error("Presumably nearly-singular input matrix in .solve_t_cholfn().")
    # Assuming the error is due to high condition number:
    mat <- regularize(mat) 
    if (.Infusion.data$options$mixturing_errorfn()) { # default set to be TRUE only for me [but not when running the checks]
      resu <- try(solve(t(chol(mat))))
      if (inherits(resu,"try-error")) {
        utils::dump.frames(dumpto = "latest.solve_t_cholfn.dump", 
                           include.GlobalEnv = TRUE, to.file=TRUE ) 
        stop("solve(t(chol(mat))) failed in .solve_t_cholfn. See latest.solve_t_cholfn.dump.")
      }
    } else resu <- solve(t(chol(mat)))
  }
  return(structure(resu,logdet=sum(log(diag(resu)))))
}

seq_nbCluster <- function(nr) {do.call(Infusion.getOption("seq_nbCluster"), list(nr=nr))}

get_nbCluster_range <- function(projdata, nr= nrow(projdata), nc=ncol(projdata), 
                                nbCluster=seq_nbCluster(nr), verbose=TRUE) { # data must be the mixmod data
  maxnb <- do.call(.Infusion.data$options$maxnbCluster, list(nr=nr,nc=nc))
  maxnb <- max(maxnb,1L) # handles 0L
  if (any(nbCluster> maxnb) && verbose) {
    message(paste0("Gaussian mixture modelling constrained to at most ", 
                   maxnb," clusters by the number of columns."))
    if (identical(attr(nbCluster,"max"),TRUE)) { # then nbCluster was seq_nbCluster ... hence:
      maxnb
    } else return(nbCluster[nbCluster<=maxnb])
  } else if (identical(attr(nbCluster,"max"),TRUE)) {
    max(nbCluster)
  } else nbCluster 
}

.data_autopsy <- function(data, stat.obs) {
  maxnb <- do.call(.Infusion.data$options$maxnbCluster, list(nr=nrow(data),nc=ncol(data)))
  if (maxnb<1L) warning("Too few simulations: Gaussian mixture modelling may fail", immediate. = TRUE) 
  if (inherits(data,"data.frame")) {
    checkranges <- sapply(data,range)
  } else checkranges <- t(matrixStats::colRanges(data))
  rangewidths <- checkranges[2L,]-checkranges[1L,]
  degenerates <- rangewidths<1e-6
  if (any(degenerates)) {
    which_dg <- which(degenerates)
    if (length(degen_pars <- setdiff(names(which_dg), names(stat.obs)))) {
      warning(paste0("Parameter(s) '",paste(degen_pars,collapse="','"),"' have very narrow ranges:\n",
                     "clustering is likely to fail or be unusable. Rethink the parametrization." ),
              immediate. = TRUE)
    }
    if (length(degen_stats <- intersect(names(which_dg), names(stat.obs)))) {
      warning(paste0("Statistic(s) '",paste(degen_stats,collapse="','"),"' have very narrow ranges:\n",
                     "clustering is likely to fail or be unusable. Rethink the parametrization." ),
              immediate. = TRUE)
      ### Highly dubious effect:
      # which_out <- (stat.obs[degen_stats]>checkranges[2L,degen_stats]) | (stat.obs[degen_stats]<checkranges[1L,degen_stats])
      # if (any(which_out)) {
      #   for (st in degen_stats[which_out]) {
      #     dst <- max(abs(stat.obs[st]-data[,st]))
      #     data[,st] <- data[,st] + rnorm(n=nrow(data),sd=dst/5) ## makes it reasonably unlikely that obs will be within noised data
      #   }
      # } else  {
      #   ## stat.obs in degenerate simuls but mixmodCluster needs non-degenerate distribution 
      #   for (st in degen_stats[which_out]) data[,st] <- data[,st] + rnorm(n=nrow(data),sd=1e-06)
      # }
    }
    for (st in names(stat.obs)) { 
      tabstat <- table(data[,st])
      if (any(whichmass <- (tabstat>Infusion.getOption("repeat_stat_thr")))) { 
        print(paste("Repeated values for statistic",st,":"))
        print(tabstat[,whichmass])
      }
    }
  }
}

.safe_densityMixmod <- function(locarglist, stat.obs) {
  jointdens <- try(do.call(".densityMixmod",c(locarglist,list(stat.obs=stat.obs))),silent=TRUE) # using seed in locarglist ie that from geoOption
  if (.any_mixmodResult_issue(jointdens))  {
    if (.Infusion.data$options$mixturing_errorfn()) { # default set to be TRUE only for me
      message(cli::col_green("Error in .densityMixmod(): entering browser session as controlled by option 'mixturing_errorfn'."))
      browser()  # controlled by mixturing_errorfn = .is_devel_session
      message(cli::col_green("Dumping frames from infer_SLik_joint() for inspection..."))
      utils::dump.frames(dumpto = "latest.infer_SLik_joint.dump", to.file=TRUE ) 
      if (FALSE) { ## older code; maybe still useful for devel though debug(.densityMixmod)...
        locarglist$strategy <- .Infusion.data$options$get_mixModstrategy(nc=ncol(locarglist$data))
        jointdens <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
        jointdens <- .get_best_mixmod_by_IC(jointdens)
      }
    }
    .data_autopsy(data=locarglist$data, stat.obs=stat.obs)
    # All attempts in .densityMixmod() => no more ideas to fix problem
    if (anyNA(jointdens@parameters@mean)) stop("Rmixmod returned NaN's in parameter estimates.") # (But in RMSE bootstrap, new bootstrap samples may be analyzed if a few failed)
    # Other errors:
    stop(jointdens) # i.e stop(<error object>) # (But in RMSE bootstrap, new bootstrap samples may be analyzed if a few failed) 
  }  
  jointdens
}

.Rmixmodelize <- function(data, locnbCluster, inferredVars, statNames, initParam, 
                          stat.obs, latentVars, fittedPars, verbose, marginalize,
                          nbCluster) {
  if (verbose$most) cat(paste0("Joint density modeling: ",nrow(data)," points"))
  cat_xpctd_nbClu <- verbose$most && length(locnbCluster)==1L
  models <- .do_call_wrap("mixmodGaussianModel",list(listModels=Infusion.getOption("mixmodGaussianModel")))
  if (FALSE && length(locnbCluster)==1L && locnbCluster==1L) {
    # The idea is to allow a modeling by two clusters with identical 'orientation' (identical D_k's)
    # when two unconstrained clusters is still not 'safely' possible  
    # (this constraints is the one that economizes more dfs than other ones.)
    #
    # This runs, at least in a preliminary sense [screen messages are not fully consistent], 
    #  but with no obvious benefits so far.
    nc <- length(c(inferredVars,statNames))
    if (3L *nrow(data)> (nc+1L)*(nc+2L)) { # threshold ~(3/2) * dfs of the constrainedmodel with two clusters 
      locarglist <- list(data=data[,c(inferredVars,statNames)],nbCluster=2L, 
                         seed=Infusion.getOption("mixmodSeed") , 
                         models=.do_call_wrap("mixmodGaussianModel",
                                              list(listModels="Gaussian_pk_Lk_D_Ak_D")))
    }
  } else {
    locarglist <- list(data=data[,c(inferredVars,statNames)],nbCluster=locnbCluster, 
                       seed=Infusion.getOption("mixmodSeed") , models=models, initParam=initParam)
    if (cat_xpctd_nbClu) cat(paste0(", and given nbCluster=",locnbCluster))
  }
  # locarglist$strategy is controlled within .densityMixmod()
  jointdens <- .safe_densityMixmod(locarglist, stat.obs)
  # plotCluster(jointdens,data=locarglist$data,variable1="theta_p",variable2="theta") # var1: stat (prediction of projection); var2: actual param
  if (cat_xpctd_nbClu) {
    if (jointdens@nbCluster!=locnbCluster) cat(paste0("; only ",jointdens@nbCluster," clusters inferred"))
  } else if (verbose$most) cat(paste0("; ",jointdens@nbCluster," clusters"))
  if (verbose$pedantic) if (jointdens@nbCluster==max(nbCluster$jointdens)) message("Inferred # of clusters= max of allowed range.")
  if ( length(latentVars)) {
    completedens <- jointdens 
    jointdens <- .marginalize_Rmixmod(completedens, colNames= c(inferredVars,statNames),
                                      For=c(fittedPars,statNames)) # marginalize over  latentvars
  } else completedens <- NULL
  if (marginalize) {
    pardens <- .marginalize_Rmixmod(jointdens, colNames= c(fittedPars,statNames),
                                    For=fittedPars) # marginalize over statNames
    # defective but should be sufficient for prediction. Some garde-fou:
    pardens@proba <- matrix(NA_real_)
    pardens@parameters@nbFreeParam <- NA_integer_
  } else {
    locarglist <- list(data=data[,fittedPars,drop=FALSE], nbCluster=nbCluster$pardens, 
                       seed=123 , models=models)
    pardens <- try(do.call(".densityMixmod",c(locarglist,list(stat.obs=stat.obs))),silent = TRUE)
    if (inherits(pardens,"try-error")) {
      # ignore any non-default infer_SLik_joint(.,nbCluster) argument
      locarglist$nbCluster <- get_nbCluster_range(projdata=data)
      pardens <- do.call("mixmodCluster",locarglist)
      pardens <- .get_best_mixmod_by_IC(pardens) 
    }                               
    if (verbose$most) cat(paste("; parameter modeling: ",pardens@nbCluster," clusters"))
  }
  #plotCluster(pardens,data=data[,fittedPars]) to plot @results[[1L]] which is a 'MixModResults', not a 'mixmodCluster' object.
  # Rmixmod::plot(<mixmodCluster object>) with Rmixmod:: to plot from any envir, not only the global one
  clu_params <- list(logproportions=log(jointdens@parameters["proportions",]),
                     jointdens_means=t(jointdens@parameters["mean",]),
                     pardens_means=t(pardens@parameters["mean",]),
                     solve_t_chol_sigma_lists=list(
                       pardens=lapply(pardens@parameters["variance"], .solve_t_cholfn),
                       jointdens=lapply(jointdens@parameters["variance"], .solve_t_cholfn)
                     ) 
  )
  if (length(latentVars)) {
    clu_params$completedens_means <- t(completedens@parameters["mean",])
    clu_params$solve_t_chol_sigma_lists$completedens <- lapply(completedens@parameters["variance"], .solve_t_cholfn)
  }
  if (verbose$most) cat(paste0(".\n"))
  list(jointdens=jointdens, pardens=pardens, completedens=completedens, clu_params=clu_params)
}

.infer_SLik_joint <- function(data, ## reference table ~ abc
                            stat.obs,
                            logLname=Infusion.getOption("logLname"), ## not immed useful
                            Simulate=attr(data,"Simulate"), ## may be NULL
                            nbCluster=seq_nbCluster(nr=nrow(data)),
                            using=Infusion.getOption("mixturing"),
                            verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                         pedantic=FALSE,
                                         final=FALSE),
                            marginalize=TRUE,
                            constr_crits=NULL,
                            projectors=NULL,
                            is_trainset=NULL,
                            initParam=NULL,
                            latentVars # info kept in $colTypes element of return value. 
) {
  if ( ! is.data.frame(data)) {
    stop(paste("'object' is not a data.frame.\n Did you mean to call infer_logLs() rather than infer_Slik_joint() ?"))
  }
  if (is.null(attr(data,"LOWER"))) { # 
    warning('Some required attributes, such as "LOWER", seem to be missing from the "data" [see "Value" in help("add_reftable")].\n Further execution could fail.', 
            immediate. = TRUE)
  }
  
  proj_attr <- attr(data,"projectors")
  if ( ! is.null(proj_attr)) {
    attr(data,"projectors") <- NULL # remove 'pointer' 
    if ( ! is.null(projectors)) {
      warning(
        paste("call with both projected 'data' and 'projectors' argument is ambiguous.\n",
              "'projectors' will be ignored."), immediate. = TRUE)
    }
    projectors_here <- proj_attr
    # and the data must be the projected ones...
    reftable_raw <- attr(data,"raw_data") 
    attr(data,"raw_data") <- NULL # remove this binding as we have the 'reftable_raw' binding
    proj_data <- data
  } else if (is.null(projectors)) {
    # no projection in workflow
    projectors_here <- NULL
    reftable_raw <- NULL
    proj_data <- data
  } else {
    # projector argument, and the reftable has not been previously projected
    projectors_here <- projectors
    reftable_raw <- data
    if (is.null(is_trainset)) stop("'is_trainset' argument must be set.")
    proj_data <- .project_reftable_raw(data, projectors = projectors_here, use_oob=TRUE, 
                         is_trainset = is_trainset)
  }
  
  if (is.null(attr(stat.obs,"raw_data") )) {# if stat.obs has not been previously projected
    if (! is.null(projectors_here)) stat.obs <- .project_reftable_raw(stat.obs, projectors = projectors_here, 
                                                        use_oob=FALSE, is_trainset = FALSE)
  } else if ( ! is.null(projectors)) warning(
    paste("call with both projected 'stat.obs' and 'projectors' argument is ambiguous.\n",
          "'projectors' will be ignored."), immediate. = TRUE)
  
  
  if ( ! is.null( cn <- colnames(stat.obs))) {
    if (nrow(stat.obs)>1L) stop("'stat.obs' has several rows, but it should instead be a vector.")
    message("Note: 'stat.obs' should be a numeric vector, not a matrix or data.frame. Converting...")
    raw_data <- attr(stat.obs,"raw_data")
    if (is.data.frame(stat.obs)) {
      stat.obs <- unlist(stat.obs)
      raw_data <- unlist(raw_data)
    } else {
      stat.obs <- drop(stat.obs)
      raw_data <- drop(raw_data)
    }
    attr(stat.obs,"raw_data") <- raw_data
  }
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$pedantic)) verbose$pedantic <- FALSE
  if (is.null(verbose$final)) verbose$final <- FALSE
  if (is.null(verbose$MAF)) verbose$MAF <- verbose$most # MAF=2L => most verbose
  if (length(unique(colnames(data)))<ncol(data)) {
    stop("Some columns of the 'data' have identical names. Please correct.")
  }
  statNames <- names(stat.obs)
  if (is.null(statNames)) stop("'stat.obs' should be a named numeric vector: provide names.")
  if (nrow(data)<20000L) check_raw_stats(na.omit(data),statNames = statNames,remove = FALSE, verbose=FALSE) # silent if not problem detected, but still verbose otherwise
  if (is.null(allPars <- attr(data,"allPars"))) {  # the ONE point where allPars attr is used.
    allPars <- setdiff(colnames(data),c(statNames,latentVars)) ## first iteration
  } ## else do not add the cumul_iter variable to the parameters !
  isVar_Pars <- apply(data[,allPars,drop=FALSE], 2L, function(v) length(unique(range(v)))>1L)
  fittedPars <- names(which(isVar_Pars))
  inferredVars <- c(fittedPars, latentVars) # actual pars and latent variables
  fixedPars <- names(which( ! isVar_Pars))
  # Syntax to remove non-essential attributes otherwise inherited from the data:
  fixedPars <- data.frame(as.matrix(data[1L,fixedPars,drop=FALSE])) ## has both names and values
  # This makes it easier to create a light 'subobject'
  
  rownames(fixedPars) <- NULL## avoids warning on rownames when later used in cbind()
  
  if (using=="mafR") {
    message('using="mafR" is interpreted as ="MAFmix".')
    using <- "MAFmix"
  }
  if ( ! length(grep("u.mafR|c.mafR",using))) {
    nbCluster <- eval(nbCluster) # allowing evaluation of quoted default arg passed from refine.default()
    if ( identical(nbCluster,"max")) nbCluster <- max(get_nbCluster_range(data))
    # 
    if (! is.list(nbCluster)) nbCluster <- list(jointdens=nbCluster, pardens=nbCluster)
    locnbCluster <- nbCluster$jointdens
    if (length(locnbCluster)==1L) { # typical of RMSE bootstrap sample. Reduce by 1 or more
      locnbCluster <- max(get_nbCluster_range(projdata=data, nbCluster=seq(locnbCluster)))
    } else locnbCluster <- get_nbCluster_range(projdata=data, nbCluster=seq(max(locnbCluster))) # generate a more standard range
  } else locnbCluster <- NULL
  if (using=="mclust") {
    if ( ! "package:mclust" %in% search()) stop("'mclust' should be loaded first.")
    locarglist <- list(data=data[,c(inferredVars,statNames)],nbCluster=locnbCluster, 
                       seed=Infusion.getOption("mixmodSeed"))
    if (verbose$most) cat(paste0("Joint density modeling: ",nrow(data)," points"))
    jointdens <- .densityMclust(data=data[,c(inferredVars,statNames)],
                                stat.obs,nbCluster=locnbCluster)
    if (verbose$most) cat(paste0("; ",jointdens$G," clusters"))
    if (verbose$pedantic) if (jointdens$G==max(nbCluster$jointdens)) message("Inferred # of clusters= max of allowed range.")
    if ( length(latentVars)) {
      completedens <- jointdens 
      jointdens <- .marginalize_mclust(completedens, colNames= c(inferredVars,statNames),
                                       For=c(fittedPars,statNames), over=c(latentVars,statNames)) # marginalize over  latentvars
    } else completedens <- NULL
    if (marginalize) {
      pardens <- .marginalize_mclust(jointdens, colNames= c(inferredVars,statNames),
                              For=fittedPars, over=statNames)
      # defective but should be sufficient for prediction. Some garde-fou:
      #pardens@proba <- matrix(NA_real_)
    } else {
      pardens <- .densityMclust(data=data[,inferredVars],
                                  stat.obs,nbCluster=nbCluster$pardens)
      if (verbose$most) cat(paste0("; parameter modeling: ",pardens$G," clusters"))
    }
    resu <- list(jointdens=jointdens, pardens=pardens, completedens=completedens,
                 clu_params=list(logproportions=log(jointdens$parameters$pro)))
    pdl <- vector("list", pardens$G)
    sigma <- resu$pardens$parameters$variance$sigma
    for (it in seq_len(pardens$G)) pdl[[it]] <- .solve_t_cholfn(sigma[,,it])
    jdl <- vector("list", jointdens$G)
    sigma <- resu$jointdens$parameters$variance$sigma
    for (it in seq_len(jointdens$G)) jdl[[it]] <- .solve_t_cholfn(sigma[,,it])
    solve_t_chol_sigma_lists <- list(pardens=pdl, jointdens=jdl) 
    if (length(latentVars)) {
      cdl <- vector("list", completedens$G)
      sigma <- resu$completedens$parameters$variance$sigma
      for (it in seq_len(completedens$G)) cdl[[it]] <- .solve_t_cholfn(sigma[,,it])
      solve_t_chol_sigma_lists$completedens <- cdl
    }
    resu$clu_params$solve_t_chol_sigma_lists <- solve_t_chol_sigma_lists
    if (verbose$most) cat(paste0(".\n"))
  } else if (using=="xLLiM") {
    if (verbose$most) cat(paste0("Conditional density modeling: ",nrow(data)," points"))
    summstats <- data[,statNames,drop=FALSE]
    RGPpars <- data[,inferredVars,drop=FALSE]
    
    # conddens <- .gllim(responses=t(summstats), predictors=t(RGPpars), nbCluster=nbCluster$jointdens)
    
    gllimobj <- .wrap_gllim(RGPpars=t(RGPpars), summstats=t(summstats), nbCluster=nbCluster$jointdens)
    
    if (verbose$most) cat(paste0("; ",ncol(gllimobj$c)," clusters.\n"))
    #
    resu <- list(gllimobj=gllimobj)    
  } else if (length(grep("u.mafR|c.mafR|MAFmix",using))) {
    time1 <- Sys.time()
    resu <- .calc_all_MAFs(data, statNames, latentVars, fittedPars, inferredVars, 
                           verbose=verbose, using=using) 
    resu$all_MAF_time <- round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1) 
    if (length(grep("MAFmix",using))) {
      resu$MGMjointdens <- .Rmixmodelize(
        data=data, locnbCluster=locnbCluster, inferredVars=inferredVars, 
        statNames=statNames, initParam=initParam, stat.obs=stat.obs, 
        latentVars=latentVars, fittedPars=fittedPars, verbose=verbose, 
        marginalize=marginalize, nbCluster=nbCluster)$jointdens
    }
    resu$load_MAFs_info <- new.env(parent = emptyenv())
  } else {   # default: Rmixmod
    resu <-  .Rmixmodelize(
      data=data, locnbCluster=locnbCluster, inferredVars=inferredVars, 
      statNames=statNames, initParam=initParam, stat.obs=stat.obs, 
      latentVars=latentVars, fittedPars=fittedPars, verbose=verbose, 
      marginalize=marginalize, nbCluster=nbCluster)
  }
  attr(resu,"EDFstat") <- "[see this string in infer_SLik_joint()]" ## 
  resu$logLs <- structure(data,stat.obs=stat.obs,Simulate=Simulate) ## as in infer_surface.logLs
  
  if ( ! is.null(reftable_raw)) {
    attr(reftable_raw,"cumul_n") <- c(0L, nrow(reftable_raw))
    resu$reftable_raw <- reftable_raw
  }
  resu$projectors <- projectors_here
  
  resu$`Infusion.version` <- packageVersion("Infusion")
  resu$colTypes <- list(allPars=allPars, ## keeps the order of the columns
                        inferredVars=inferredVars,
                        fittedPars=fittedPars,
                        fixedPars=fixedPars,
                        latentVars=latentVars,
                        logLname=logLname,
                        statNames=statNames)
  resu$lower <- sapply(data[,fittedPars,drop=FALSE],min) ## so now $lower excludes latent vars while $LOWER includes them
  resu$upper <- sapply(data[,fittedPars,drop=FALSE],max) # ... used by MSL -> optim
  if (any(is.na(c(resu$lower,resu$upper)))) stop("NA(s) in c(lower,upper))") 
  resu$LOWER <- attr(data,"LOWER") # used in many places
  if (is.null(resu$LOWER)) resu$LOWER <- sapply(data[,inferredVars,drop=FALSE],min)
  resu$UPPER <- attr(data,"UPPER") # ... used by confintAll
  if (is.null(resu$UPPER)) resu$UPPER <- sapply(data[,inferredVars,drop=FALSE],max)
  attr(resu,"Infusion.version") <- packageVersion("Infusion")
  class(resu) <- c("SLik_j",class(resu))
  resu$using <- using
  resu$constr_crits <- constr_crits
  return(resu)
}

# API version has no initParam argument and sets the latentVars arg for the data attribute:
infer_SLik_joint <- function(data, ## reference table ~ abc
                              stat.obs,
                              logLname=Infusion.getOption("logLname"), ## not immed useful
                              Simulate=attr(data,"Simulate"), ## may be NULL
                              nbCluster=seq_nbCluster(nr=nrow(data)),
                              using=Infusion.getOption("mixturing"),
                              verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                           pedantic=FALSE,
                                           final=FALSE),
                              marginalize=TRUE,
                              constr_crits=NULL,
                              projectors=NULL,
                              is_trainset=NULL
) {
  mc <- match.call(expand.dots=TRUE) 
  mc[[1L]] <- get(".infer_SLik_joint", asNamespace("Infusion"), inherits=FALSE) 
  mc["latentVars"] <- list(attr(data,"latentVars")) # [] on mc seems to work as on a list
  resu <- eval(mc,envir = parent.frame())
  if (length(resu$colTypes$statNames) < 
      length(resu$colTypes$fittedPars)) warning("Fewer statistics than fitted parameters:\n  some parameters may not be identifiable.")
  resu
}

## was infer_SLik_joint -> predict -> predict.SLik_j -> predict.MixmodResults
# except that now                                    -> .get_dens_from_GMM.Rmixmod -> predict.dMixmod
# But predict.MixmodResults is this still used: predict(gofdens...)
predict.MixmodResults <- function(object, newdata,log=TRUE, solve_t_chol_sigma_list,
                                  logproportions=log(object@parameters@proportions), ...) {
  Sobs_activeBoundaries <- atb <- freqs <- NULL ## FR->FR tempo
  nbCluster <- object@nbCluster
  if (is.null(nrow(newdata)) ) {
    newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  } else if (inherits(newdata,"data.frame")) {
    newdata <- matrix(unlist(newdata,recursive = FALSE, use.names = FALSE), nrow=nrow(newdata),
                      dimnames=list(NULL,colnames(newdata)) ) # newdata <- as.matrix(newdata)
  }
  if (log) { 
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(newdata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- logproportions[k] + 
          .fast_dmvnorm(newdata, object@parameters["mean", k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=log)
      }
      mixture <- matrixStats::rowLogSumExps(density)
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture + atb*log(freqs)
    } else mixture <- atb*log(freqs)
  } else {
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(newdata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- object@parameters["proportions", k] * 
          .fast_dmvnorm(newdata, object@parameters["mean", k], solve_t_chol_sigma_list[[k]],log=log)
      }
      mixture <- .rowSums(density,m=nrow(density), n=ncol(density)) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture) # not fully protected from NaN's - for newdata *far* from the represented hypervolume.
}

.get_dens_from_GMM.Rmixmod <- function(X, # parameters only 
                                  object, 
                                  tstat.obs, # 1-row matrix as otherwise more cases should be considered for cbind'ing
                                  log, 
                                  which, # "safe" version ignores, by correcting, spuriously high logL in area of low parameter density.
                                  solve_t_chol_sigma_lists=object$clu_params$solve_t_chol_sigma_lists, # added to allow optimization via iterative smmothing (which does not work),
                                  thr_info=.get_thr_info(object),
                                  ...)  {
  if (is.null(dim(X))) {
    # dim(X) <- c(1L, length(X)) drops dimnames (sigh)
    X <- t(X)
    newjointX <- cbind(X,tstat.obs) # cbind two 1-row matrices
  } else {
    if (length(intersect(colnames(X),colnames(tstat.obs)))) stop("'X' should contain only parameters, not summary statistics")
    newjointX <- cbind(X,tstat.obs[rep(1,nrow(X)),,drop=FALSE]) # cbind two nrow(X)-row matrices
  }
  if (which!="parvaldens") {
    jointvaldens <- predict.dMixmod(object$jointdens,
                            newdata=newjointX,
                            solve_t_chol_sigma_list=solve_t_chol_sigma_lists$jointdens,
                            logproportions=object$clu_params$logproportions,
                            clu_means=object$clu_params$jointdens_means,
                            log=log,...)
  }
  if (which=="jointvaldens") return(jointvaldens) # not used...
  parvaldens <- predict.dMixmod(object$pardens,
                        newdata=newjointX, # statistics will be ignored
                        solve_t_chol_sigma_list=solve_t_chol_sigma_lists$pardens,
                        logproportions=object$clu_params$logproportions,
                        clu_means=object$clu_params$pardens_means,
                        log=log,...) 
  if (which=="parvaldens") return(parvaldens)
  # ELSE: "lik", or "safe" for safe version of "lik" using thr_dpar
  if (log) {
    condvaldens <- jointvaldens - parvaldens
    if (which=="safe") {
      ddens <- parvaldens-thr_info$thr_dpar
      negddens <- (ddens<0)
      if (any(negddens)) {
        dlogl <- condvaldens-thr_info$reft_maxlogl # perhaps by setting a sightly higher thr
        # high would allow a bit more extrapol and thus exploration? bc currently it may induce self-reinforcing 
        # high pardens in suboptimal regions.... but the effect seems opposite on exploration...
        posdlogl <- dlogl>0
        #lowextrapol <- negddens & ! posdlogl
        #condvaldens[lowextrapol] <- condvaldens[lowextrapol]+ ... # not clear what to do
        highextrapol <- negddens & posdlogl # ...................... potentially conservative intervals.............
        # make it continuous wrt to dlogl, but strongly compensating in most cases
        condvaldens[highextrapol] <- condvaldens[highextrapol]+ sqrt(dlogl[highextrapol])*ddens[highextrapol]
        # extrapol <-  posdlogl & ! negddens
        # condvaldens[extrapol] <- condvaldens[extrapol]- sqrt(dlogl[highextrapol])
        attr(condvaldens,"lowdens") <- negddens
      }
    }
  } else {
    condvaldens <- jointvaldens/parvaldens
    if (which=="safe") condvaldens <- condvaldens*pmin(1,parvaldens/exp(object$thr_info$thr_dpar))
  }
  return(condvaldens) # vector if X was a matrix
}

.get_densv <- function(X, object,  # SLik_j object OR subset (useful for parallelisation)
                      tstat.obs, log, which, solve_t_chol_sigma_lists=object$clu_params$solve_t_chol_sigma_lists) {
  if (inherits(object$jointdens,"dMixmod")) {
    somedens_value <- .get_dens_from_GMM.Rmixmod(X=X, object=object, tstat.obs=tstat.obs, log=log, which=which, 
                                            solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)
    ## RMSEs computation uses names of CI bounds, to be retained here, to name its results
    # apply(newdata,1L, .get_dens_from_GMM.Rmixmod, object=object, tstat.obs=tstat, log=log, which=which)
  } else if (inherits(object$jointdens,"MAF")) { # using="mafR" => .calc_all_MAFs() was called
    somedens_value <- .predict_MAF(object=object, X=X, tstat.obs=tstat.obs, which=which)
  # } else if (inherits(object$conddens,"MAF")) { # using="devel" => 
  #   # only conddens <- .wrap_MAF_cond() + Rmixmod were called, 
  #   # so that previous test inherits(object$jointdens,"MAF") was FALSE
  #   somedens_value <- .predict_MAF_devel(object=object, X=X, tstat.obs=tstat.obs, which=which)
  } else if (inherits(object$jointdens,"dMclust")) {
    somedens_value <- .predict_SLik_j_mclust(object=object, X=X, tstat.obs=tstat.obs, log=log, which=which, 
                                        solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)
  } else if (inherits(object$gllimobj,"gllim")) {
    somedens_value <- .get_dens_from_GMM.gllim(X=X, 
                                          object=object,
                                          tstat.obs=tstat.obs, log=log, which=which)
  } 
  names(somedens_value) <- rownames(X)
  somedens_value # logl or safe logl or RGPdens...
}

.predict_SLik_j_mclust <- function(
    object, X, tstat.obs, log, which, 
    solve_t_chol_sigma_lists=object$clu_params$solve_t_chol_sigma_lists,
    thr_info=.get_thr_info(object))  {
  if (is.null(dim(X))) {
    # dim(X) <- c(1L, length(X)) drops dimnames (sigh)
    X <- t(X)
    newjointX <- cbind(X,tstat.obs) # cbind two 1-row matrices
  } else {
    if (length(intersect(colnames(X),colnames(tstat.obs)))) stop("'X' should contain only parameters, not summary statistics")
    newjointX <- cbind(X,tstat.obs[rep(1,nrow(X)),,drop=FALSE]) # cbind two nrow(X)-row matrices
  }
  
  if (which!="parvaldens") jointvaldens <- predict(
    object$jointdens, newdata=newjointX, 
    solve_t_chol_sigma_list=solve_t_chol_sigma_lists$jointdens,
    logproportions=object$clu_params$logproportions,
    log=log) # ! order of elements in newdata must be that of fittedPars as in object$jointdens$data
  if (which=="jointvaldens") return(jointvaldens)
  parvaldens <- predict(object$pardens,newdata=X, 
                        solve_t_chol_sigma_list=solve_t_chol_sigma_lists$pardens,
                        logproportions=object$clu_params$logproportions,
                        log=log) 
  if (which=="parvaldens") return(parvaldens)
  if (log) {
    condvaldens <- jointvaldens - parvaldens
    if (which=="safe") {
      ddens <- parvaldens-thr_info$thr_dpar
      negddens <- (ddens<0)
      if (any(negddens)) {
        dlogl <- condvaldens-thr_info$reft_maxlogl
        posdlogl <- dlogl>0
        highextrapol <- negddens & posdlogl
        # make it continuous wrt to dlogl, but strongly compensating in most cases
        condvaldens[highextrapol] <- condvaldens[highextrapol]+ sqrt(dlogl[highextrapol])*ddens[highextrapol]
        attr(condvaldens,"lowdens") <- negddens
      }
    }
  } else {
    condvaldens <- jointvaldens/parvaldens
    if (which=="safe") condvaldens <- condvaldens*pmin(1,parvaldens/exp(object$thr_info$thr_dpar))
  }
  if (any(is.infinite(condvaldens))) warning("any(is.infinite(condvaldens)) is TRUE: expect an error.") 
  return(condvaldens)
}

.calc_outsideness <- function(v, constr_crits) {
  outsideness <- eval(constr_crits, envir = as.list(v))
  outsideness <- outsideness[outsideness>0]
  sqrt(sum(outsideness*outsideness))
}

predict.SLik_j <- function(
    object, 
    newdata, ## requests new fittedPars values! 
    log=TRUE, 
    which="safe", # tested v2.1.127 with no default value to check that fn is internally called with explicit values.
    tstat= t(get_from(object,"stat.obs")), # 1-row matrix...
    solve_t_chol_sigma_lists=object$clu_params$solve_t_chol_sigma_lists, # added to allow optimization via iterative smmothing (which does not work)
    constr_tuning=FALSE, 
    ...) {
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (is.null(colnames(newdata))) colnames(newdata) <- object$colTypes$fittedPars
  # not useful bc predict.MixmodResults does not handle it:
  #newdata <- cbind(newdata, t(replicate(nrow(newdata),get_from(object,"stat.obs"))))
  logl <- .get_densv(X=newdata, object=object, tstat.obs=tstat, log=log, which=which, 
                    solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)
  if ( which %in% c("lik","safe") && ! is.null(constr_crits <- object$constr_crits) ) {
    if (is.infinite(constr_tuning)) {
      is_outside <- apply(as.data.frame(newdata),1L, 
                          FUN=function(v) any(eval(object$constr_crits, envir = as.list(v))>0))
      logl[is_outside] <- -Inf
    } else if (is.na(constr_tuning)) {
      is_outside <- apply(as.data.frame(newdata),1L, 
                          FUN=function(v) any(eval(object$constr_crits, envir = as.list(v))>0))
      logl[is_outside] <- -.Machine$double.xmax
    } else if (constr_tuning) {
      outsidenesses <- apply(as.data.frame(newdata),1L, FUN=.calc_outsideness, 
                             constr_crits=object$constr_crits)
      logl <- logl - constr_tuning*outsidenesses
    }
  }
  if ( ! is.null(object$prior_logL)) logl <- logl + object$prior_logL(newdata)
  return(logl)
}

.wrap_calc_ecdf_t <- function(object, parm, nsim, cluster_args, nb_cores, ...) {
  newprojStats <- simulate(object, nsim=nsim) 
  if (inherits(object$jointdens,"MAF") && 
      .set_cluster_type(cluster_args, nb_cores)$spec >1L) {
    ## multiple processes may then use gpu... hmmm. 
    ## Provide the MAF objects:
    if (is.null(object$load_MAFs_info$pwd)) {
      stop("save_MAFs() should have been run so that child processes can load_MAFs().")
    } 
    torch_device <- .get_py_MAF_handle()$device$type
    simfun <- function(it) {
      config_mafR(torch_device=torch_device)
      load_MAFs_info <- object$load_MAFs_info
      object <- load_MAFs(object,prefix=load_MAFs_info$prefix, ext=load_MAFs_info$ext)
      summlik <- summLik(object, parm=NULL, data=newprojStats[it,, drop=FALSE], which="safe")
      attr(summlik,"profpt")[parm]
    }
    fit_env <- list(torch_device=torch_device)
  } else  {
    simfun <- function(it) {
      summlik <- summLik(object, parm=NULL, data=newprojStats[it,, drop=FALSE], which="safe")
      attr(summlik,"profpt")[parm]
    }
    fit_env <- list()
  }
  iterator <- seq_len(nsim)
  dim(iterator) <- c(1L,nsim)
  object <- .shrink(object, ...) 
  ecdf_t <- spaMM::dopar(newresp = iterator, fn = simfun, fit_env=fit_env,
                         cluster_args=cluster_args, nb_cores=nb_cores, ... )
  dim(ecdf_t) <- c(nsim,1L)
  ecdf_t
}

confint.SLik_j <- function(object, parm, ## parm is the parameter which CI is sought 
                         level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE), 
                         nsim=0L, # type="perc",
                         reset=TRUE,
                         cluster_args=NULL, nb_cores=NULL,
                         type="perc", # NB "raw" bootstrap CI does not exist...
                         ...) {
  ## 'Bartlett' requires à distribution of the LR stat (ecdf_2lr). 
  ## Some bootstrap computations only return an ecdf_t.
  ## confint.SLik_j -> confintAll will use ecdf_2lr if present; (but mind default reset=TRUE)
  ## confint.SLik_j -> boot.ci will use if ecdf_t is present.
  ## ecdf_t mais be available from a previous confint.SLik_j(., nsim) computation 
  ##   that stored its result in the CI object;
  ## But ecdf_t may also be available from a previous .ecdf_2lr computation
  ##   that stored its results in the bootreps object. ecdf_t may then be copied 
  ##   in the CI object.
  force(parm) # otherwise missing parm => cryptic error from .confintAll() (!)
  fittedPars <- object$colTypes$fittedPars

  # tryVM <- FALSE && # .VenzonM stops workflows (try H_7from17) ____F I X M E____ retry...
  #   .is_devel_session() && 
  #   length(fittedPars)>4L && 
  #   requireNamespace("proliks",quietly = TRUE) &&
  #   packageVersion("proliks")>="0.1.4" &&
  #   all(eigen(object$MSL$hessian,only.values = TRUE)$values>1e-3)
  # if (tryVM) {
  #   cat(cli::col_green("Trying VenzonM's method..."))
  #   optim_info <- list(est=object$MSL$MSLE, lower= object$LOWER[fittedPars], upper= object$UPPER[fittedPars])
  #   objectivefn <- function(v) {
  #     names(v) <- fittedPars
  #     return( - predict(object,newdata=v, which="safe", constr_tuning=FALSE)) ## removed log...   
  #   }
  #   .VenzonM <- get(".VenzonM", asNamespace("proliks"))
  #   ci <- try(.VenzonM(X=optim_info, 
  #                      nlogf=objectivefn,
  #                      parm=parm, prob=level, verbose=FALSE)$interval)
  # } 
  # 
  if (nsim>1L) { # provide  ecdf_t, and ecdf_2lr if "Bartlett" in type.
    
    # # next line may be useful only for devel hack calling confint.SLik_j on a non SLik_j method ?
    # if ( ! inherits(object,"SLik_j")) stop("Bootstrap CI computations are implemented only for objects of class 'SLik_j'.")
    
    if (is.null(object$bootCIenv)) {
      errmess <- paste("SLik_j object has no 'bootCIenv' element.\n",
                       "Was it created using an old version of Infusion?\n.",
                       "If so, running MSL() on it may fix this issue.")
      stop(errmess)
    }
    CIobject <- object$CIobject
    CIs <- CIobject$CIs
    if ("Bartlett" %in% type) { # then GET the necessary info. 
      bootreps_list <- object$bootCIenv$bootreps_list
      if (reset || is.null(ecdf_t <- bootreps_list[[parm]]$ecdf_t)) {
        MSLE <- object$MSL$MSLE
        # .ecdf_2lr() has MAF-handling code
        blob_2lr <- .ecdf_2lr(object, BGP=MSLE, h0_pars=parm,                       
                              nsim=nsim, nb_cores=max(1L,nb_cores), cluster_args=cluster_args, 
                              ... )
        bootreps_list[[parm]] <- blob_2lr
        assign("bootreps_list",bootreps_list,envir=object$bootCIenv)
        ecdf_t <- blob_2lr$ecdf_t
      } else if (nrow(ecdf_t)!=nsim) message("Saved object used as reset=FALSE: 'nsim' argument ignored.")
      ecdf_2lr <- bootreps_list[[parm]]$ecdf_2lr
    } else { # GET ecdf_t (only)
      if (reset || is.null(ecdf_t <- .get_ci_info(CIs, parm)$ecdf_t)) {
        ecdf_t <- .wrap_calc_ecdf_t(object=object, parm=parm, nsim=nsim, 
                                    cluster_args=cluster_args, nb_cores=nb_cores, ...)
      } else if (length(ecdf_t)!=nsim) message("Saved object used as reset=FALSE: 'nsim' argument ignored.")
      ecdf_2lr <- NULL
    }
    ## ecdf_t is a 1-col matrix
  } else ecdf_t <- ecdf_2lr <- NULL
  # if ( ( ! tryVM) || inherits(ci,"try-error")) {
    ci <- .confintAll(object=object, parm=parm, ## parm is the parameter whose CI is sought 
                givenmax = object$MSL$maxlogL,
                coverage=level,
                level= - qchisq(level,df=1)/2, ## df=1 for 1D profile; /2 to compare to logLR rather than 2 logLR  
                ecdf_2lr=ecdf_2lr, # => Bartlett correction if ecdf_2lr not NULL
                verbose=verbose,fixed=fixed,which=which, 
                cluster_args=cluster_args, nb_cores=nb_cores,...)
    #
    if ( ! is.null(ecdf_t) && length(bootCItype <- setdiff(type,"Bartlett"))  ) {
      ci$ecdf_t <- ecdf_t
      ci$bootCI <- boot.ci(
        boot.out = list(t=ecdf_t, t0=object$MSL$MSLE[parm], R=nrow(ecdf_t), sim="parametric"), 
        type=bootCItype, conf = level)
      CIs[[parm]] <- ci
      assign("CIs", CIs, envir=CIobject)
    }
  # }
  ci
}

.get_ci_info <- function(CIs, parm, level, prof_y=NULL) {
  ci4parm <- CIs[[parm]]
  if (inherits(ci4parm,c("ci_info","NULL"))) {
    ci_info <- ci4parm
  } else { # ci4parm expected to be a structured list of ci_info's with names corresponding to level
    # but at least in transitional stage, or in old fitobject, it might be only a list 
    if (inherits(ci4parm[[1]],"ci_info")) { # check that structure is OK
      if (length(ci4parm)>1L || # level should be specified in this case...
          ! missing(level) ) {
        ci_info <- ci4parm[[paste0(level)]] 
      } else ci_info <- ci4parm[[1]]
    } else if ("interval" %in% names(ci4parm)) ci_info <- ci4parm
  }
  if ( (! is.null(ci_info)) && 
       (! is.null(prof_y)) && 
       is.null(ci_info$prof_y)) ci_info$prof_y <- prof_y
  ci_info
}

.as_cis4parm <- function(ci_or_cis, level) {
  # ci4parm expected to be a structured list of ci_info's with names corresponding to level
  # but at least in transitional stage, or in old fitobject, it might be only a list 
  if (inherits(ci_or_cis[[1]],"ci_info")) { 
    # Nothing to change
  } else if ("interval" %in% names(ci_or_cis)) {
    class(ci_or_cis) <- c("ci_info", class(ci_or_cis))
    ci_or_cis <- list("NA"=ci_or_cis) 
  } else if (is.null(ci_or_cis)) {
    ci_or_cis <- list()
  } 
  ci_or_cis
}

.assign_ci_info <- function(ci_info, CIobject, parm, level) {
  CIs <- CIobject$CIs
  cis4parm <- .as_cis4parm(ci_or_cis=CIs[[parm]], level)
  cis4parm[[paste0(level)]] <- ci_info 
  CIs[[parm]] <- cis4parm
  assign("CIs", CIs, envir=CIobject)
}

allCIs <- function(object,level=0.95, verbose=TRUE, ...) {
  if (is.null(object)) {
    warning("allCIs() called on NULL 'object.' Check call.")
    return(NULL)
  }
  CIs <- cis4level <- list()
  for (st in object$colTypes$fittedPars) {
    cis4parm <- list() # .as_cis4parm(ci_or_cis=CIs[[st]], level)
    cis4parm[[paste0(level)]] <- 
      cis4level[[st]] <- confint(object,st,level=level, verbose=verbose,...) 
    CIs[[st]] <- cis4parm
    # .assign_ci_info(ci_info, object$CIobject, parm=st, level=level)
    # => does not make sense since CIobject is to be created by the present function
    # (and it has elements beyond CIs, preventing some changes in the code)
    # The 'bounds' etc won't match the structure of the 'CIs' as they are for one level only
  }
  lowers <- lapply(cis4level,function(li) {li$lowerpar})
  if ( ! is.null(lowers)) names(lowers) <-  paste("low.",names(lowers),sep="")
  uppers <- lapply(cis4level,function(li) {li$upperpar})
  if ( ! is.null(uppers)) names(uppers) <-  paste("up.",names(uppers),sep="")
  # the list elements are either numeric vectors or asingle NA... 
  ordre <- order(c(seq_len(length(lowers)),seq_len(length(uppers))))
  bounds <- c(lowers,uppers)[ordre]
  checkvec <- function(vec) {if (identical(vec,NA)) {return(NULL)} else {return(vec)} }
  whichNAs <- unlist(lapply(bounds,function(vec) {identical(vec,NA)}))
  missingBounds <- names(bounds[whichNAs])
  bounds <- do.call(rbind,bounds[ ! whichNAs])
  invisible(list(CIs=CIs, bounds=bounds, missingBounds=missingBounds, level=level, warn=NULL)) 
}

refine.SLik_j <- function(object, method=NULL, using=object$using, ...) {
  if (is.null(method)) method <- "mixmodCluster" ## no clear effect but distinct from SLik case
  if (is.null(using)) {
    using <- Infusion.getOption("mixturing") # in case an object from version < 2.0.4 is input...
  } else if (using=="mafR") {
    message('using="mafR" is interpreted as ="MAFmix".')
    using <- "MAFmix"
  }
  # if it's not NULL, any user's non-default for 'using' is passed to the next call and heeded. 
  resu <- refine.default(object, surfaceData=object$logLs, method=method, using=using, ...)
  # resu may simply be a table of new parameter values...:
  # if (inherits(resu,"SLik_j")) { # ... but then, MSL() has always been called by refine.default and this is not necess?
  #   assign("bootreps_list", list(), envir=resu$bootLRTenv) 
  #   assign("bootreps_list", list(), envir=resu$bootCIenv) 
  # }
    
  resu
}

summary.SLik_j <- summary.SLik

print.SLik_j <- function(x,...) {summary.SLik_j(x,...)} 

.rparam_dMixmod_around_focal <- function(object, # *dMixmod* object
                                 focal, 
                                 solve_t_chol_sigma_list,
                                 size) {
  varNames <- object@varNames
  if (is.null(dim(focal))) {
    datanames <- names(focal) # before they are erased by the next dim() assignment (!)
    dim(focal) <- c(1L, length(focal))
    if (is.null(datanames)) {
      colnames(focal) <- varNames
    } else colnames(focal) <- datanames
  } else if (inherits(focal,"data.frame")) {
    focal <- matrix(unlist(focal,recursive = FALSE, use.names = FALSE), nrow=1L,
                    dimnames=list(NULL,colnames(focal)) ) # newdata <- as.matrix(newdata)
  }
  if ( !is.null(Sobs_activeBoundaries <- object@Sobs_activeBoundaries)) { # if the Sobs matches some boundaries, identified in Sobs_activeBoundaries
    # then looks whether the newdata match all of the boundaries met by Sobs
    boundsdata <-  focal[,names(Sobs_activeBoundaries),drop=FALSE]
    atb <- apply(boundsdata,1L,`==`,y=Sobs_activeBoundaries)
    # If not (i.e. only partial matching) the object does not predict correctly the density => warning + heuristic patch:
    if (!all(atb)) {warning("'predict.dMixmod' cannot compute joint out-of-boundary density from conditional at-boundary density. ")}
    freq <- object@freq
    freqs <- atb*freq+(1-atb)*(1-freq) ## uses 1-freq instead of the unknown density of the variable(s) in boundaries 
    densitydata <- focal[,varNames,drop=FALSE]
  } else if ( !is.null(simuls_activeBoundaries <- object@simuls_activeBoundaries)) {
    ## do not reduce data in this case
    densitydata <- focal 
    ## only for the warning:
    boundsdata <-  focal[,names(simuls_activeBoundaries),drop=FALSE]
    atb <- apply(boundsdata,1L,`==`,y=simuls_activeBoundaries)
    if (any(atb)) {
      warning("'predict.dMixmod' cannot compute conditional at-boundary density from joint out-of-boundary density.")
      # return value is the same as for predict(,tcstat.obs=<newdata>) 
    }
  } else densitydata <- focal[,varNames,drop=FALSE]
  
  nbCluster <- object@nbCluster
  
  density <- numeric(nbCluster)
  for (k in 1:nbCluster) {
    density[k] <- # without the 'prior' cluster proba bc aim is not to predict the cluster given these 'prior' freqs but to modify this 'prior'  
      .fast_dmvnorm(densitydata, object@parameters["mean", k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=TRUE)
  }
  k <- which.max(density)
  trypoints <- rmvnorm(size, 
                       mean=object@parameters["mean", k], 
                       sigma= object@parameters["variance",k])
  colnames(trypoints) <- varNames
  trypoints
}

focal_refine <- function(object, focal, size, plotprof=TRUE, ...) {
  fittedPars <- object$colTypes$fittedPars
  prof <- profile(fitted=object,value=focal, return.optim=TRUE)
  fullfocal <- c(focal,prof$solution)[fittedPars]
  posteriordens <- .conditional_Rmixmod(object$jointdens, given=get_from(object,"stat.obs"), expansion=1) 
  solve_t_chol_sigma_list <- lapply(posteriordens@parameters["variance"], .solve_t_cholfn)
  trypoints <- .rparam_dMixmod_around_focal(posteriordens, focal = fullfocal, solve_t_chol_sigma_list=solve_t_chol_sigma_list, 
                                            size=size)
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],,drop=FALSE]
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],,drop=FALSE]
  trypoints <- data.frame(trypoints)
  trypoints <- cbind(trypoints,object$colTypes$fixedPars) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering and remove polluting things
  object <- refine(object, trypoints=trypoints, ...)
  if (plotprof) {
    pars <- names(focal)
    if (length(pars)==1L) {
      plot1Dprof(object, pars = pars, ...)
    } else if (length(pars)==2L) {
      plot2Dprof(object, pars = pars, ...)
    } 
  }
  object
}

# Handling box and 'constr_crits' constraints: 
.apply_par_constraints_LOWUP <- function(trypoints, object, fittedPars=object$colTypes$fittedPars) {
  inLOWUP <- TRUE
  for (vv in fittedPars) {
    par_vals <- trypoints[,vv]
    inLOWUP <- inLOWUP & par_vals>object$LOWER[vv] & par_vals<object$UPPER[vv]
  }
  trypoints <- trypoints[inLOWUP,,drop=FALSE]
  if (nrow(trypoints) && ! is.null(constr_crits <- object$constr_crits)) {
    obeys_par_constr <- apply(as.data.frame(trypoints),1L, function(v) all(eval(constr_crits, envir = as.list(v))<0))
    trypoints <- trypoints[obeys_par_constr,,drop=FALSE]
  }
  trypoints
}

.apply_par_constraints_lowup <- function(trypoints, object, fittedPars=object$colTypes$fittedPars,
                                         upper=object$upper, lower=object$lower, locrange=upper-lower) {
  for (st in fittedPars) {
    # Differs from .apply_par_constraints_LOWUP() here:
    locrange_st <- locrange[st]
    atlower <- which((trypoints[,st]-lower[st])/locrange_st<0.0005)
    atupper <- which((upper[st]-trypoints[,st])/locrange_st<0.0005)
    trypoints[atlower,st] <- trypoints[atlower,st] + runif(length(atlower),0,locrange_st*0.002)
    trypoints[atupper,st] <- trypoints[atupper,st] - runif(length(atupper),0,locrange_st*0.002)
  }
  if ( ! is.null(constr_crits <- object$constr_crits)) {
    obeys_par_constr <- apply(as.data.frame(trypoints),1L, 
                              function(v) all(eval(constr_crits, envir = as.list(v))<0))
    trypoints <- trypoints[obeys_par_constr,,drop=FALSE]
  }
  trypoints
}

.simVVV_rmvt <- function(parameters, n, seed = NULL, ...) {
  if (!is.null(seed)) 
    set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if (any(is.na(parameters[c("mean", "variance")])) || any(is.null(parameters[c("mean", 
                                                                                "variance")]))) {
    warning("parameters are missing")
    return(structure(matrix(as.double(NA), n, d + 1), modelName = "VVV"))
  }
  pro <- parameters$pro
  if (is.null(pro)) 
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins = G)
  x <- matrix(0, n, d)
  for (k in 1:G) {
    m <- ctabel[k]
    sigma <- parameters$variance$sigma[,,k] # for 1x1 matrix drop=FALSE => D array, vs default drops all dims, hence
    if (is.null(dim(sigma))) dim(sigma) <- c(1L,1L) # ugly...
    if (m>0) x[clabels == k, ] <- rmvt(n=m, 
                              delta= mu[, k], 
                              sigma= sigma)
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "VVV")
}

.sample_in_absol_constrs <- function(object, ceil_size, sampledPars,
                                     density, # should be sub-dimensional if there are givenpars
                                     givenpars=NULL,
                                     norm_or_t, # .wrap_rmvnorm or .wrap_rmvt 
                                     seed=NULL
                                     ) {
  if (inherits(density,"dMclust")) {
    if (identical(norm_or_t,.wrap_rmvt)) {
      trypoints <- .simVVV_rmvt(parameters=density$parameters,n=ceil_size) #matrix
    } else if (identical(norm_or_t,.wrap_rmvnorm)) {
      trypoints <- do.call("sim", c(density[c("modelName", "parameters")],list(n=ceil_size))) #matrix
    } else stop("Unhandled 'norm_or_t' value for dMclust object.")
    trypoints <- trypoints[,-1,drop=FALSE]  # remove group column
    if ( ! is.null(givenpars)) trypoints <- cbind(trypoints, t(givenpars)[rep(1,nrow(trypoints)),,drop=FALSE]) # matrix
    colnames(trypoints) <- c(sampledPars,names(givenpars)) # sim() does not keep track of param names
  } else if (inherits(density,"gllim")) {
    trypoints <- .simulate.gllimX(density, size=ceil_size, drop=TRUE, parNames=sampledPars)
    if ( ! is.null(givenpars)) trypoints <- cbind(trypoints, t(givenpars)[rep(1,nrow(trypoints)),,drop=FALSE]) # matrix
  } else if (inherits(density,"MAF")) { # "mafR" variants
    if (attr(density,"which")=="I_postdens") { 
      givens <- givenpars[,object$colTypes$statNames, drop=FALSE] # expected to match t(get_from(object,"stat.obs"))
    } else givens <- givenpars 
    trypoints <- .simulate.MAF(density, nsim=ceil_size, 
                               given=givens, 
                               object=object) # ___F I X M E___? provide non-default batchsize, 
                                              # fn of object's reftable dim and of py_handle's $gpu_memory as for the fit?
    fixed_here <- intersect(colnames(trypoints),names(givenpars))
    if (length(fixed_here)) {
      # Occurs in plot2Dprof in basic tests... 
      # profile.SLik_j -> .safe_init(NULL base_inits) -> .inits_from_postdens() -> here
      # That might theoretically happen in other case where no good base_inits are found. 
      # Cf also if .get_instr_postdens() was called with a non-stat.obs 'given'
      # hence if .inits_from_postdens() is called with 'givenpars'
      # a rare subcase of .safe_init()?
      # warning("Subperfect code for MAF in .sample_in_absol_constrs() with 'givenpars'.") # ____F I X M E____ it's still so despite the '#'
      for(st in fixed_here) trypoints[,st] <- givenpars[st] # QUICK PATCH
    } 
  } else { # using="Rmixmod", or "devel" for pardens, instr_postdens or a conditional distrib in instr_postdens
    trypoints <- .simulate.MixmodResults(density, size=ceil_size, drop=TRUE, seed=seed,
                                         norm_or_t=norm_or_t) 
    if ( ! is.null(givenpars)) trypoints <- cbind(trypoints, t(givenpars)[rep(1,nrow(trypoints)),,drop=FALSE]) # matrix
  }
  # For some time I've more or less assumed that parameter order would not matter;
  # but then in .safe_init I rbind() points that would then have different parameter orders:
  if ( ! is.null(givenpars)) trypoints <- trypoints[,object$colTypes$fittedPars,drop=FALSE] 
  trypoints <- .apply_par_constraints_LOWUP(trypoints, object=object, fittedPars=colnames(trypoints))
  trypoints
}

# Test cases: (use verbose$cloud_parm)
# weibull_noproj (true logL known) replicate 89 (originally also an issue of coordinate transfo...)
# N_7from17 replicate 2 for log1p.t23.; 
# N_7from17 replicate 46 for log1p.t34. maybe worst param-replicate combination of the 200 replicates)
## Points with high weights will be preferentially retained
.calc_filltop_weights <- function(object, trypoints, 
                                  log_rel_weight_penal, 
                                  exploration_fac=.Infusion.data$options$exploration_fac, 
                                  target_LR, 
                                  which # currently which="lik" in both calls
                                  ) {
  logLs <- predict(object,trypoints,which=which, constr_tuning=FALSE) 
  
  flatnd_max <- max(logLs)-target_LR # lower logL value determined by target LR
  upperslice <- ( logLs> flatnd_max ) ## defines an upper slice of logLs values.
  
  # log_rel_weight_penal by default penalizes for pardens, AND a bit by instr_postdens
  w_u <- - log_rel_weight_penal[upperslice] ## weights as opposite of penalizations

  # It suffice that one point has quite high weight (perhaps spuriously) to increase largely rejection proba of all other points,
  # => truncates high weights
  if ( (n_up <- length(which(upperslice)))>101L) {
    qs <- quantile(ecdf(w_u),probs=c(0.01,0.99))
    w_u <- pmin(w_u,qs[2]) 
    w_u <- pmax(w_u,qs[1]) # not sur this one is useful
  }
  
  w_u <- exp(w_u-max(w_u)) # w_u is now a probability and the max _realized_ w_u is 1

  if (any( ! upperslice)) {
    # Weights in ! upperslice decrease with decreasing (logLs-penal). 
    w_l <- logLs[ ! upperslice] - log_rel_weight_penal[ ! upperslice] 
    w_l <- exp(w_l-max(w_l)) # w_l a proba with realized max=1
    if ( ! is.null(exploration_fac)) { # allows some control of fraction of upperslice points
      # Quick test on N_7from17 #2 07/2024 shows exploration extremely efficient when exploration_fac is NULL
      # higher explo fac -> corr_fac more likely to be >1 -> more likely correction of w_l 
      corr_fac <- exploration_fac*sum(w_u)/sum(w_l)
      if (corr_fac>1) w_l <- w_l*corr_fac ## 'if sum(w_upper) too high, increase w_lower'
      ## so that sum(w_l)=exploration_fac*sum(w_u).
      ## * The opposite correction w_u <- w_u/corr_fac gave bad results.
      ## * This is not perfect: some (possibly many) w_l's become >1, which 
      ## no longer allows proportional sampling, and the lower slice remains 
      ## under-explored relative to the target defined by exploration fac. 
      ## * Nevertheless, this seems better than nothing, and still allowed the kind of 
      ## lateral expansion needed for N_7from17 replicate 2 
      ## [ Speculations: does lateral expansion
      ## result from indirect effect on importance of moreweight computation?.
      ## <=> high w_l's are important for good exploration? (so NO correction when corr_fac<1)
    }
    weights <- numeric(length(logLs))
    weights[upperslice] <- w_u
    weights[ ! upperslice] <- w_l
  } else weights <- w_u 
  
  ### Add a few top points
  # if (which=="safe") {
    local_logLs <- logLs
  # } else {
  #   safe_logLs <- predict(object,trypoints,which="safe", constr_tuning=FALSE)
  # }
  ## Adding 'top_pts' drives initial upwards push of the cloud of points.
  ## Some sort of balance between upward and lateral expansions may be required.
  top_pts <- head(order(local_logLs[upperslice], decreasing=TRUE), 
                  min(5, ceiling(length(upperslice)/10))) # beware number of trypoints >> target size
  logLs_u <- local_logLs[upperslice][top_pts]
  weights[upperslice][top_pts] <- pmax(weights[upperslice][top_pts], exp(logLs_u-logLs_u[1]))
  
  list(weights=weights, freq_upperslice=sum(upperslice)/length(upperslice))
}

.get_size_first_iter <- function(object) {
  size_first_iter <- attr(object$logLs,"n_first_iter")
  if (is.null(size_first_iter)) size_first_iter <- length(which(object$logLs$cumul_iter==1L)) 
  size_first_iter
}

.get_workflow_design_2.1.75 <- function(npar,cumn_over_maxit=NULL) {
  final_reft_size <- max(2000L*npar, 80L*npar*(2*npar+1L))
  refine_blocksize <- max(1000L, # for npar=1...
                         (1000L*( final_reft_size %/% 3000)))
  divisors_of_1000 <- c(2L,4L,5L,8L,10L,20L)
  subblock_nbr <- max(divisors_of_1000[((npar+5L)%/%3L) >= divisors_of_1000]) # O(npar/3)
  subblock_size <- refine_blocksize/subblock_nbr 
  init_reft_size <- min(refine_blocksize/5L,1000L)
  if (is.null(cumn_over_maxit)) stop("Explicit value of 'cumn_over_maxit' argument needed.")
  rm("divisors_of_1000")
  mget(ls(), environment())
}

.get_workflow_design_2.1.87 <- function(npar, 
                                        final_reft_size=NULL,
                                        refine_blocksize=NULL,
                                        subblock_nbr=NULL,
                                        init_reft_size=NULL,
                                        cumn_over_maxit=NULL
) {
  if (is.null(final_reft_size)) final_reft_size <- 1000L * round( pmax(2000L*npar, 240L*npar^2)/1000 )
  # Save 3 or 4 sizes:
  if (is.null(refine_blocksize)) {
    refine_blocksize <- max(1000L, # for npar=1...
                            (1000L*( final_reft_size %/% 3000)))
  }
  if (is.null(subblock_nbr)) {
    divisors <- seq(2,max(4L,sqrt(refine_blocksize/100)))
    subblock_nbr <- max(divisors[(refine_blocksize %% divisors)==0L])
    rm("divisors")
  }
  subblock_size <- refine_blocksize/subblock_nbr 
  if (is.null(init_reft_size)) init_reft_size <- min(refine_blocksize/5L,1000L)
  reftable_sizes <- seq(refine_blocksize,final_reft_size, refine_blocksize) 
  if (is.null(cumn_over_maxit)) cumn_over_maxit <- TRUE
  first_refine_ntot <- 2L*subblock_size
  mget(ls(), environment())
}

.get_workflow_design_2.1.112 <- function(npar, 
                                 final_reft_size=NULL,
                                 refine_blocksize=NULL,
                                 subblock_nbr=NULL,
                                 init_reft_size=NULL,
                                 cumn_over_maxit=NULL
) {
  #  if (is.null(final_reft_size)) final_reft_size <- 1000L * round( pmax(2000L*npar, 160L*npar^2)/1000 )
  if (is.null(final_reft_size)) final_reft_size <- 1000L * pmax(3L*npar-1L, round( pmax(2000L*npar, 200L*npar^2)/1000))
  # maybe 1000L * round( pmax(2500L*npar, 160L*npar^2)/1000) ?
  # Save 3 or 4 sizes:
  if (is.null(refine_blocksize)) {
    refine_blocksize <- max(1000L, # for npar=1...
                            (1000L*( final_reft_size %/% 3000)))
  }
  if (is.null(subblock_nbr)) {
    divisors <- seq(2,max(4L,refine_blocksize/(100*log(refine_blocksize))))
    subblock_nbr <- max(divisors[(refine_blocksize %% divisors)==0L])
    rm("divisors")
  }
  subblock_size <- refine_blocksize/subblock_nbr 
  if (is.null(init_reft_size)) init_reft_size <- min(refine_blocksize/5L,1000L)
  reftable_sizes <- seq(refine_blocksize,final_reft_size, refine_blocksize) 
  if (is.null(cumn_over_maxit)) cumn_over_maxit <- TRUE
  first_refine_ntot <- 2L*subblock_size
  mget(ls(), environment())
}

.get_workflow_design_2.1.162 <- function(npar, 
                                 n_proj_stats=npar,
                                 final_reft_size=NULL,
                                 refine_blocksize=NULL,
                                 subblock_nbr=NULL,
                                 init_reft_size=NULL,
                                 cumn_over_maxit=NULL,
                                 test_fac=NULL
) {
  nvar <- npar+n_proj_stats
  if (is.null(final_reft_size)) final_reft_size <- 1000L * pmax(round(nvar*3/2-1), 
                                                                round( pmax(1000L*nvar, 50L*nvar^2)/1000))
  if ( ! is.null(test_fac)) {
    final_reft_size <- final_reft_size*test_fac
    if (is.null(refine_blocksize)) {
      refine_blocksize <- max(1000L*test_fac, # for npar=1...
                              (1000L*test_fac*( final_reft_size %/% (3000*test_fac))))
    }
  } else {
    # maybe 1000L * round( pmax(2500L*npar, 160L*npar^2)/1000) ?
    # Save 3 or 4 sizes:
    if (is.null(refine_blocksize)) {
      refine_blocksize <- max(1000L, # for npar=1...
                              (1000L*( final_reft_size %/% 3000)))
    }
  }
  if (is.null(subblock_nbr)) {
    divisors <- seq(2,max(4L,refine_blocksize/(100*log(refine_blocksize))))
    subblock_nbr <- max(divisors[(refine_blocksize %% divisors)==0L])
    rm("divisors")
  }
  subblock_size <- refine_blocksize/subblock_nbr 
  if (is.null(init_reft_size)) init_reft_size <- min(refine_blocksize/5L,1000L)
  reftable_sizes <- seq(refine_blocksize,final_reft_size, refine_blocksize) 
  if (is.null(cumn_over_maxit)) cumn_over_maxit <- TRUE
  first_refine_ntot <- 2L*subblock_size
  # !!! : smaller steps:
  subblock_nbr <- 2L*subblock_nbr
  subblock_size <- subblock_size/2L
  mget(ls(), environment())
}

.get_workflow_design <- function(npar, 
                                 n_proj_stats=npar,
                                 n_latent=0L,
                                 final_reft_size=NULL,
                                 refine_blocksize=NULL,
                                 subblock_nbr=NULL,
                                 init_reft_size=NULL,
                                 cumn_over_maxit=NULL,
                                 test_fac=NULL
) {
  nvar <- npar+n_proj_stats+n_latent # dim of completedens
  if (is.null(final_reft_size)) final_reft_size <- 
      1000L * round(nvar*3/2-1)
    ## more simuls from nvar>=30 (certainly preferable for precision but less simple) 
    # 1000L* pmax(round(nvar*3/2-1), round( pmax(1000L*nvar, 50L*nvar^2)/1000))
  if ( ! is.null(test_fac)) {
    final_reft_size <- final_reft_size*test_fac
    if (is.null(refine_blocksize)) {
      refine_blocksize <- max(1000L*test_fac, # for npar=1...
                              (1000L*test_fac*( final_reft_size %/% (3000*test_fac))))
      # => ~ refine_blocksize= 1000L * (2npar-1)
    }
  } else {
    # maybe 1000L * round( pmax(2500L*npar, 160L*npar^2)/1000) ?
    # Save 3 or 4 sizes:
    if (is.null(refine_blocksize)) {
      refine_blocksize <- max(1000L, # for npar=1...
                              (1000L*( final_reft_size %/% 3000)))
    }
  }
  if (is.null(subblock_nbr)) {
    divisors <- seq(2,max(4L,refine_blocksize/(100*log(refine_blocksize))))
    subblock_nbr <- max(divisors[(refine_blocksize %% divisors)==0L])
    rm("divisors")
  }
  subblock_size <- refine_blocksize/subblock_nbr 
  if (is.null(init_reft_size)) init_reft_size <- min(refine_blocksize/5L,1000L)
  reftable_sizes <- seq(refine_blocksize,final_reft_size, refine_blocksize) 
  if (is.null(cumn_over_maxit)) cumn_over_maxit <- TRUE
  first_refine_ntot <- 2L*subblock_size
  if (first_refine_ntot != 
      reftable_sizes[1L]) reftable_sizes <- c(first_refine_ntot, reftable_sizes)
  if (final_reft_size != 
      tail(reftable_sizes,1L)) reftable_sizes <- c(reftable_sizes, final_reft_size)
  # !!! : smaller steps:
  subblock_nbr <- 2L*subblock_nbr
  # sapply(seq(20)*2L, function(npar) get_workflow_design(npar)$subblock_nbr) :
  # 8  8 10 14 18 22 26 30 34 38 42 40 40 50 50 50 60 56 50 60
  subblock_size <- subblock_size/2L
  # sapply(seq(20)*2L, function(npar) get_workflow_design(npar)$subblock_size) :
  # 125 375 500 500 500 500 500 500 500 500 500 575 625 540 580 620 550 625 740 650
  mget(ls(), environment())
}

get_workflow_design <- function(npar, 
                                n_proj_stats=npar,
                                n_latent=0L,
                                final_reft_size=NULL,
                                refine_blocksize=NULL,
                                subblock_nbr=NULL,
                                version=packageVersion("Infusion"),
                                cumn_over_maxit=NULL,
                                test_fac=NULL) {
  version <- as.package_version(version)
  if (version < "2.1.76") {
    wf <- .get_workflow_design_2.1.75(npar, cumn_over_maxit=cumn_over_maxit)
  } else if (version < "2.1.88") {
    wf <- .get_workflow_design_2.1.87(npar, cumn_over_maxit=cumn_over_maxit)
  } else if (version < "2.1.163") {
    wf <- .get_workflow_design_2.1.162(npar, cumn_over_maxit=cumn_over_maxit)
  } else wf <- .get_workflow_design(npar, n_proj_stats=n_proj_stats, 
                                    n_latent=n_latent,
                              final_reft_size=final_reft_size,
                              refine_blocksize=refine_blocksize,
                              subblock_nbr=subblock_nbr,
                              cumn_over_maxit=cumn_over_maxit,
                              test_fac=test_fac)
  if (version < "2.1.184.2") { ## older version of maxnbCluster():
    ## # of df for the different models: CeleuxG95
    ## Full MGM model has P = (nc*(nc+3)/2 +1)G-1) params but this has used 
    ## Q = (nc*(nc+1)/2 +1)G-1) (quite minor confusion). 
    ## A practically saturated model can be defined by setting max G so that Q ~ nr,
    ## but this will give very poor fits and failures.
    ## If we set max G so that 2 Q ~ nr
    ## 4 clusters, 6 cols => 2Q=174, but nr=200, G=4 still sometimes fails. 
    ## If we set max G so that 3 Q ~ nr
    ## G * ("1+"nc*(nc +1L))"*(3/2)" - 1 ~ 258. 
    ## ie G=4 is first attained for maxnbCluster(nr = [4*(1+6*7)*3/2-1 = 258],nc=6)   = 4
    message("Infusion's 'maxnbCluster' options reset to pre-v2.1.184.2 definition.")
    Infusion.options(maxnbCluster= function(projdata, nr= nrow(projdata), nc=ncol(projdata)) {
      nr_corr <- (nr*2L) %/% (3L) 
      floor((nr_corr+1L)/(nc*(nc+1L)+1L)) 
    })
  } else if (version < "2.1.186.2") { # transient version: quite good on MVNcovmat_identif,
    # but not so convincing on B_13from17
    message("Infusion's 'maxnbCluster' options set to transient ~ v2.1.185-186 definition.")
    Infusion.options(maxnbCluster= function(projdata, nr= nrow(projdata), nc=ncol(projdata)) {
      (nr+6L)%/%(nc*(nc+1L)*3L+6L) 
    })
  }
  wf$version <- version
  wf
}

.calc_ceil_size <- function(object, target_size) {
  rparam_info <- object$rparam_info # typically non-NULL after first iteration of first refine;
  # but it is theoretically possible that .sample_filltop_pts() generated no valid points, 
  # in which case this is NULL after first refine [in practice this has occurred only in a simulation 
  # workflow where the initial reftable did not respect its manually set LOWER, UPPER attrs
  #  => added check in add_reftable. ]
                                     
  if ( ! is.null(rparam_info)) { # then estimate the ceiling size (=> large if old_freq_full_rejection low)
    old_freq_full_rejection <- rparam_info$freq_full_rejection
    # Use smaller target sizes if $freq_runif_weights_thr*freq_runif_rejection falls beyond 1
    #and even smaller ones when in the initial subblocks or the initial 'refine_blocksize'
    adj_runif_thr <- .Infusion.data$options$freq_runif_weights_thr
    nr1 <- nrow(object$logLs)
    npar <- length(object$colTypes$fittedPars)
    refine_blocksize <- get_workflow_design(npar)$refine_blocksize 
    nr2 <- refine_blocksize * ((nr1 %/% refine_blocksize)+1L) # reconstruct the likely final size of a refine
    mxnbClufn <- .Infusion.data$options$maxnbCluster
    maxnb1 <- do.call(mxnbClufn, list(nr=nr1, nc=2L*npar))
    maxnb2 <- max(1L, do.call(mxnbClufn, list(nr=nr2, nc=2L*npar)))
    adj_runif_thr <- adj_runif_thr * min(1,maxnb1/maxnb2) # this reduces adj_target_size when 
    # maxnb1<maxnb2 ie when nr1<nr2 ie reftable size < refine_blocksize
    # instr_post sampling tends to be activated in that case.
    adj_target_size <- as.integer(ceiling( # adjust for lower target when poor success of runif() rejection
      target_size *min(1, 
                       max(0.1, 
                           adj_runif_thr * rparam_info$freq_runif_rejection)))) 
    
    # The option 'freq_full_rejection_bnd' controls bnded_old_freq_full_rejection -> ceil_size
    # which controls both steps of sampling: from logL and from postdens.
    bnded_old_freq_full_rejection <- max(.Infusion.data$options$freq_full_rejection_bnd,
                                         old_freq_full_rejection) # max() to avoid explosion of # of sampled points
    target_ceil_size <- ini <- adj_target_size/bnded_old_freq_full_rejection # a sort of expectation of the number of point to try to obtain enough good ones
    # A sort of predicted bound of an interval for the same number, computed recursively: 
    locfn <- function(ceil_size) {
      qnorm(0.001, mean=ceil_size*bnded_old_freq_full_rejection, 
            sd=sqrt(ceil_size*bnded_old_freq_full_rejection*(1-bnded_old_freq_full_rejection))
      ) -adj_target_size }
    while(locfn(ini)<0) ini <- ini*2
    ceil_size <- uniroot(f = locfn, interval=c(ini/2,ini))$root
    ceil_size <- max(ceiling(100/bnded_old_freq_full_rejection), ceil_size) # useful if low target_size 
  } else {
    adj_target_size <- target_size
    target_ceil_size <- ceil_size <-  ceiling(.Infusion.data$options$default_ceil_size_fac*target_size)
    rparam_info <- list(freq_in_par_constrs=1)
  }
  rparam_info$target_size <- target_size # pure info? Yes at the time of writing
  rparam_info$adj_target_size <- adj_target_size
  rparam_info$target_ceil_size <- target_ceil_size
  rparam_info$ceil_size <- as.integer(ceil_size) # 
  rparam_info
}
# The return value typically contains:
# freq_upperslice=0.39, # pure info, does not control the algo
# freq_in_par_constrs=0.547, # PREVIOUS length(weights)/ceil_size_blob$ceil_size
# freq_runif_rejection=0.0184, # PREVIOUS length(instr_good)/length(weights) # success of runif() rejection step
# freq_full_rejection=0.0101, # PREVIOUS length(instr_good)/ceil_size_blob$ceil_size
# target_size=1378, # input target size
# adj_target_size=., # ... 
# target_ceil_size=137061, #  target_size/bnded_old_freq_full_rejection
# ceil_size=148893, # safe upper bound for target_ceil_size: actual sampling size in MV mixture model

.sample_profile_in_absol_constrs <- function(
    object, refill, cluster_args, target_LR, 
    fittedPars=object$colTypes$fittedPars,
    safe #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
) {
  prof_info <- plot1Dprof(object, do_plot=FALSE, cluster_args=cluster_args, safe=safe, verbose=FALSE)
  profpts <- do.call(rbind, lapply(prof_info$profiles,`[[`, "profpts"))
  profLR <- do.call(c, lapply(prof_info$profiles,`[[`, "y"))
  profpts <- profpts[,,drop=FALSE]
  profgood <- which(profLR > - target_LR) # keep not-too-large negative values
  profLR[order(profLR[profgood])[1:5]]
  if (length(profgood)>refill) profgood <- profgood[sample(length(profgood), refill)]
  profpts <- profpts[profgood,,drop=FALSE]
  profpts <- .apply_par_constraints_lowup(profpts, object, fittedPars=fittedPars) # why do I need that?
  profpts
}


.get_instr_postdens <- function(object, given=get_from(object,"stat.obs")) {
  # instrumental dens for sampling will now be the 'instr_postdens'
  if (inherits(object$jointdens,"dMixmod")) { # using="Rmixmod" or "devel"
    instr_postdens <- .conditional_Rmixmod(object$jointdens, given=given, expansion=1)
  } else if (inherits(object$jointdens,"MAF")) { # using="mafR"
    if (length(grep("MAFmix",object$using))) {
      instr_postdens <- .conditional_Rmixmod(object$MGMjointdens, given=given, expansion=1)
    } else instr_postdens <- object$postdens
  } else if (inherits(object$jointdens,"dMclust")) {
    instr_postdens <- .conditional_mclust(object$jointdens, given=given, expansion=1)
  } else if (inherits(object$gllimobj,"gllim")) {
    # warning("No implementation of instr_postdens for gllim clustering: using pardens.")
    fittedPars <- object$colTypes$fittedPars 
    givenParNames <- intersect(names(given),fittedPars)
    if (length(givenParNames)) { # that's super ugly...
      instr_postdens <- .conditional_gllimobj(object$gllimobj, fittedPars=fittedPars, 
                                             given=given[givenParNames], expansion=1)
    } else instr_postdens <- object$gllimobj
  }
  instr_postdens
}

.calc_instr_densities <- function(object, 
                                  instr_dens, # typically the parameters' instrumental posterior 
                                  moretrypoints) {
  ## Selection of points ~uniformly in top of logl hill ('good' points)
  if (inherits(instr_dens,"dMixmod")) { # using="Rmixmod" or "devel"
    instr_densv <- predict.dMixmod( 
      instr_dens,
      newdata=moretrypoints, 
      solve_t_chol_sigma_list=lapply(instr_dens@parameters["variance"], .solve_t_cholfn),
      logproportions=log(instr_dens@parameters["proportions",]),
      clu_means=t(instr_dens@parameters["mean",]),
      log=TRUE) 
  } else if (inherits(instr_dens,"dMclust")) {
    nbCluster <- instr_dens$G
    solve_t_chol_sigma_list <- vector("list", nbCluster)
    sigma <- instr_dens$parameters$variance$sigma
    for (it in seq_len(nbCluster)) solve_t_chol_sigma_list[[it]] <- .solve_t_cholfn(sigma[,,it])
    
    instr_densv <- predict.dMclust( 
      instr_dens, newdata=moretrypoints, 
      solve_t_chol_sigma_list=solve_t_chol_sigma_list,
      logproportions = log(instr_dens$parameters$pro),
      log=TRUE) 
  } else if (inherits(object$pardens,"MAF")) { # using="mafR"
    # This does not handle zero-row 'moretrypoints'.
    instr_densv <- .predict_MAF(X=moretrypoints, # parameters only 
                                            object=object, 
                                            tstat.obs=NULL, 
                                            which="parvaldens") 
  } else if (inherits(object$gllimobj,"gllim")) {
    instr_densv <- predict(object,moretrypoints,which="parvaldens") # just a quick fix (bis)
  }
instr_densv
}

.sample_envelope_pts <- function(object, 
                                 instr_dens, # .get_instr_postdens(object)
                                 ceil_size_blob, fittedPars, target_LR, 
                                 fill, 
                                 norm_or_t=.Infusion.data$options$norm_or_t) {
  sample_instr_dens <- .sample_in_absol_constrs(object, ceil_size_blob$ceil_size, fittedPars, 
                                            density=instr_dens, 
                                            norm_or_t=norm_or_t)
  if ( nrow(sample_instr_dens)) {
    instr_densv <- .calc_instr_densities(object, instr_dens, sample_instr_dens)
    pardensv <- predict(object,newdata = sample_instr_dens, which="parvaldens", constr_tuning=FALSE) # constr_tuning useless as pred_data satisfy constraints
    penal_facs <- .Infusion.data$options$penal_facs
    log_rel_weight_penal <- penal_facs["pardens"] * pardensv +
      penal_facs["instr_dens"] * instr_densv 
    moreweights_blob <- .calc_filltop_weights(object=object, trypoints=sample_instr_dens, 
                                              log_rel_weight_penal = log_rel_weight_penal,
                                              which="lik", # one case where "lik" is clearly useful. # ___F I X M E___ another topic for repeated rethink
                                              target_LR=target_LR)
    moreweights <- moreweights_blob$weights
    # When 'morepoints' are needed, we allow more exploration at the boundary of the high logL region 
    # from the 'instrumental posterior':
    morewei_runif <- runif(n=length(moreweights))<moreweights ## allows some exploration by sampling some 'bad' points
    
    moregood <- which(morewei_runif)
    
    len_mg <- length(moregood)
    if (len_mg > fill) moregood <- moregood[sample(len_mg,size = fill)]
  } else { # it oocurred... with mafR
    moregood <- c()
  }
  
  sample_instr_dens[moregood,,drop=FALSE]
}


.sample_filltop_pts <- function(object, ceil_size_blob, fittedPars, 
                                target_size, 
                                verbose, target_LR) {
  if (inherits(object$gllimobj,"gllim")) {
    trypoints <- .sample_in_absol_constrs(object, ceil_size_blob$ceil_size, fittedPars, 
                                          density=object$gllimobj, 
                                          norm_or_t=.Infusion.data$options$norm_or_t)
  } else { # 
    trypoints <- .sample_in_absol_constrs(object, 
                                          ceil_size_blob$ceil_size, # tried .../2: v2.1.113.1, N_7from17. Not better 
                                          fittedPars, 
                                          density=object$pardens, 
                                          norm_or_t=.Infusion.data$options$norm_or_t)
  }
  NROW_trypoints <- NROW(trypoints)
  ok <- (NROW_trypoints >= target_size/ceil_size_blob$freq_in_par_constrs)
  if ( verbose && ! ok ) { # but note that new points may be added below
    if ( ! is.null(constr_crits <- object$constr_crits)) {
      locmess <- paste0("Only ",NROW(trypoints)," candidate points generated satisfying parameter constraints.")  
    } else  locmess <- paste0("Only ",NROW(trypoints)," candidate points generated satisfying bound constraints.")  
    message(locmess)
  }
  
  if (NROW_trypoints) {
    ## Selection of points ~uniformly in top of logl hill ('good' points)
    pardensv <- predict(object,trypoints,which="parvaldens") 
    weights_blob <- .calc_filltop_weights(object=object, trypoints=trypoints, 
                                          log_rel_weight_penal=pardensv, 
                                          exploration_fac = NULL,
                                          which="lik", target_LR=target_LR)
    weights <- weights_blob$weights
    # low prior => high weight => preferentially retained 
    wei_runif <- runif(n=length(weights))<weights ## allows some exploration by sampling some 'bad' points
    
    instr_good <- which(wei_runif)
    len_ig <- length(instr_good)
    rparam_info <- list(
      ## INFO:
      freq_upperslice=weights_blob$freq_upperslice, # info NOT USED currently;
      freq_in_par_constrs=NROW_trypoints/ceil_size_blob$ceil_size, # info NOT USED currently;
      # might be used by removing a 'TRUE || ...':
      freq_runif_rejection=len_ig/NROW_trypoints, # potentially used to control sample_best_clu
      # => The weights may be very non-uniform so which(wei_runif) may be much smaller than that of uniform sampling.
      # Ideally sampling from the instr posterior|the 'best' cluster will bring more high-weight candidates (but that may not be so simple)
      #
      # USED: product of two previous ones
      freq_full_rejection=len_ig/ceil_size_blob$ceil_size # used to control ceil_size in next iteration.
    )
    # First iteration, high dim -> one cluster -> many predicted 
    
    if (len_ig>target_size) instr_good <- instr_good[sample(len_ig, size = target_size)]
    
    trypoints <- trypoints[instr_good,,drop=FALSE] # these ones ideall fill the top
  } else rparam_info <- NULL
  list(trypoints=trypoints, rparam_info=rparam_info)
}


.rparam_SLik_j_in_out <- function(object, 
                                target_size=NULL,
                                fittedPars=object$colTypes$fittedPars,level,
                                target_LR,
                                maxit=10L,
                                size_first_iter=.get_size_first_iter(object),
                                verbose=FALSE,
                                # cluster_args=NULL,
                                safe #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
) {
  
  ## Defaults if not useful target_size info in input:
  prev_n_iter <- max(object$logLs$cumul_iter)
  if (is.null(target_size)) target_size <- max(size_first_iter/2, 
                                               size_first_iter*(prev_n_iter+1L)-nrow(object$logLs))
  
  ceil_size_blob <- .calc_ceil_size(object, target_size) # O(target_size/bnded_old_freq_full_rejection)
  if (verbose) {
    if (.is_devel_session()) {
      cat(cli::col_green(paste0("Target subblock size: ", target_size,"; target_LR: ", target_LR,"\n")))
      ublob <- unlist(ceil_size_blob)      
      cat(cli::col_green(
        paste0(names(ublob),"=",.prettysignif(ublob),c("","","","\n"),collapse=", ")))
    } else cat("Target subblock size: ", target_size,"\n")
  }
  
  ## Generate candidates points respecting parameter constraints par drawing in pardens
  ## loop to enforce absolute constraints by rejection
  ######################## FIRST SAMPLING STEP
  RESU <- .sample_filltop_pts( # -> first .calc_filltop_weights
    object=object, ceil_size_blob=ceil_size_blob, fittedPars=fittedPars, 
    target_size=target_size*0.9, # (___F I X M E___ constant rethink) 0.8 was not better in weibull_noproj test
    verbose=verbose, target_LR=target_LR)
  trypoints <- RESU$trypoints
  NROW_selected <- NROW(trypoints)
 
  incise <- NULL
  fill <- target_size - NROW_selected
  if (fill) { # presumably TRUE since we set target_size=target_size*0.9 in the above call.
    if (.is_devel_session()) incise <-  paste0(NROW_selected)
    ######################## INSTR POSTERIOR SAMPLING STEP
    envelope_pts <- .sample_envelope_pts( # -> second .calc_filltop_weights
      object=object,       instr_dens=.get_instr_postdens(object),
      ceil_size_blob=ceil_size_blob, fittedPars=fittedPars, 
      target_LR=target_LR, fill=fill) 
    trypoints <- rbind(trypoints, envelope_pts)
    if (.is_devel_session()) incise <- 
      paste0(incise, ", + ", nrow(envelope_pts)," using instr.posterior")
  } 
  
  if (verbose) {
    if (length(incise)) incise <-  cli::col_green(paste0(" (", incise,")"))
    nr <-  nrow(trypoints) # total after possible addition of 'moretrypoints'
    if (nr < target_size) {
      cat(paste0("Only ", nr, " points generated",incise,". "))
    } else cat(paste0(nr, " points generated",incise,". "))
  }
  
  trypoints <- cbind(trypoints,object$colTypes$fixedPars) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  RESU$trypoints <- trypoints
  RESU$fill_info <- list(NROW_selected=NROW_selected, 
                         fill=if (fill) {nrow(envelope_pts)} else {0L})
  # RESU already has $rparam_info
  RESU
}

.rparam_SLik_j_B_postdens <- function( # hack ASSUMING UNIFORM PRIOR
  # tested with Rmixmod, likely pb with MAF
    object, 
    target_size=NULL,
    fittedPars=object$colTypes$fittedPars,
    verbose=FALSE,
    ...
    # level, # not used
    # target_LR,
    # maxit=10L,
    # size_first_iter=.get_size_first_iter(object),
    # safe #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
    ### cluster_args=NULL
) {
  ceil_size_blob <- .calc_ceil_size(object, target_size) # O(target_size/bnded_old_freq_full_rejection)
  if (verbose) cat("Target subblock size: ", target_size,"\n")
  
  trypoints <- .sample_in_absol_constrs(object, ceil_size_blob$ceil_size, fittedPars, 
                                        density=.get_instr_postdens(object), 
                                        norm_or_t=.wrap_rmvnorm)
  NROW_trypoints <- NROW(trypoints) 
  if (NROW_trypoints) { # Penalize by instrumental prior density:
    pardensv <- predict(object,newdata = trypoints, which="parvaldens", constr_tuning=FALSE) # constr_tuning useless as pred_data satisfy constraints
    wei <- - pardensv
    # It suffice that one point has quite high weight (perhaps spuriously) to increase largely rejection proba of all other points,
    # => truncates high weights
    # NB this may results from pardens and postdens not deduced from a single joint distrib, 
    # so possible artefact to corretc, but introduces other artefacts
    if ( length(wei)>101L) {
      qs <- quantile(ecdf(wei),probs=c(0.01,0.99))
      wei <- pmin(wei,qs[2]) 
      wei <- pmax(wei,qs[1]) # not sur this one is useful
    }
    
    wei <- exp(wei-max(wei))
    wei_runif <- runif(n=length(wei))<wei 
    good <- which(wei_runif)
    len_good <- length(good)
    if (len_good > target_size) good <- good[sample(len_good,size = target_size)]
    trypoints <- trypoints[good,,drop=FALSE]
  }
  NROW_final <- NROW(trypoints)
  rparam_info <- list(
    freq_upperslice=1, # info NOT USED currently;
    freq_in_par_constrs=NROW_final/ceil_size_blob$ceil_size, # info NOT USED currently;
    freq_runif_rejection=len_good/NROW_trypoints, # potentially used to control sample_best_clu
    freq_full_rejection=len_good/ceil_size_blob$ceil_size # used to control ceil_size in next iteration.
  )
  trypoints <- cbind(trypoints,object$colTypes$fixedPars) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  RESU <- list(trypoints=trypoints,
               rparam_info=rparam_info,
               fill_info=list(NROW_selected=NROW_final, 
                              fill=0L) )
  
  RESU
}


# generic from stats package: profile(fitted,...), returns a # *log*L
profile.SLik_j <- function(fitted, value, fixed=NULL, return.optim=FALSE, 
                           init="default", which="safe", 
                           constr_crits=fitted$constr_crits, 
                           eq_constr=NULL, ...) {
  fixedPars <- names(fixed)   
  fittedPars <- fitted$colTypes$fittedPars
  if (! is.null(fixed)) { # warning introduced in v. 2.1.17 (2023/09/01). Doc for this argument 
    # suggested it was an overlooked copy from a plot-profile function.
    warning("Use of 'fixed' argument is redundant with 'value' argument, will be deprecated, and should be avoided",
            immediate. = TRUE)
    fittedPars <- setdiff(fittedPars,fixedPars)
  }
  fittedparamnbr <- length(fittedPars) 
  parm <- names(value)
  MLval <- fitted$MSL$MSLE[parm]
  if (anyNA(MLval)) {stop(paste0("'",paste0(parm[which(is.na(MLval))],collapse="','"),
                                 "' appear(s) to be incorrect parameter name(s). Check 'parm' argument."))}
  lowval <- fitted$lower[parm]
  hival <- fitted$upper[parm]
  v <- fitted$MSL$MSLE; v[names(v)] <- NA ## create template 
  lowval <- lowval + 0.002 * (MLval - lowval)
  hival <- hival - 0.002 * (hival - MLval)
  if (fittedparamnbr == 1L) {
    v[parm] <- value
    if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
    resu <- predict(fitted,newdata=v, which=which, constr_tuning=Inf)
  } else { 
    # ___F I X M E___ ? Possible trivial improvement: in the (not well-though) case where a full param vector is provided to profile(), 
    #a simple predict() call should be enough. 
    # Currently an expensive .safe_init appears to be performed (+ other inefficiencies?)  
    profiledNames <- names(fitted$lower)
    profiledNames <- setdiff(profiledNames, c(parm,fixedPars)) # [which( ! (profiledNames %in% parm))] 
    plower <- fitted$lower[profiledNames]
    pupper <- fitted$upper[profiledNames]
    v[parm] <- value 
    constr_tuning <- FALSE
    plogL <- function(pparv) {
      v[profiledNames] <- pparv
      if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
      return(( - predict(fitted,newdata=v, which=which, constr_tuning=constr_tuning))) ## removed log... (log=TRUE is the default)  
    }
    if ( ! is.null(constr_crits)) {
      neg_ineq_constrfn <- function(pparv) {
        v[profiledNames] <- pparv
        if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
        as.numeric(eval(constr_crits, envir = as.list(v)))
      }
    } else neg_ineq_constrfn <- NULL
    
    if ( ! is.null(eq_constr)) {
      eq_constrfn <- function(pparv) {
        v[profiledNames] <- pparv
        if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
        eval(eq_constrfn, envir = as.list(v))
      }
    } else eq_constrfn <- NULL
    
    if (identical(init,"default")) {
      if (inherits(fitted,"SLik_j")) {
        init <- .safe_init(object = fitted, given=c(value, fixed), plower, pupper,
                           more_inits=fitted$MSL$MSLE,
                           constr_crits=constr_crits) # does not account for eq_constrfn=eq_constrfn
        optr <- .safe_optim(init[profiledNames], plogL, lower=plower, upper=pupper, LowUp=list(), 
                            verbose=FALSE, object=fitted, neg_ineq_constrfn=neg_ineq_constrfn,
                            eq_constrfn=eq_constrfn, ...)
      } else optr <- .safe_optim(fitted$MSL$MSLE[profiledNames], plogL, lower=plower, upper=pupper, 
                                 LowUp=list(), verbose=FALSE, object=fitted, 
                                 neg_ineq_constrfn=neg_ineq_constrfn,
                                 eq_constrfn=eq_constrfn, ...)
    } else if (is.numeric(init)) { # else 'init' input argument must be user-provided vector including profiled pars 
      #Example is p <- profile(object, value=templateh0, init=init, which="safe") in SLRT() 
      optr <- .safe_optim(init[profiledNames], plogL, lower=plower, upper=pupper, LowUp=list(), 
                          verbose=FALSE, object=fitted, neg_ineq_constrfn=neg_ineq_constrfn,
                          eq_constrfn=eq_constrfn,
                          ...)
    } else if (identical(init,"clu_means")) { # multiple optims, one from each clu mean. Could probably be removed.
      newdata <- .newdata_from_pardens_clu_means(fitted,value)
      OPT <- function(init) {.safe_optim(init[profiledNames], plogL, lower=plower, 
                                         upper=pupper, LowUp=list(), verbose=FALSE, object=fitted, ...)}
      optrS <- apply(newdata,1L, OPT)  
      optr <- optrS[[which.min(sapply(optrS,`[[`, i="objective"))]]
    }
    
    if(return.optim) {
      optr$value <- - optr$objective 
      #optr$par <- optr$solution 
      return(optr)
    } else {
      resu <- - optr$objective
      attr(resu,"solution") <- optr$solution
      return(resu)
    }
  }
} # *log*L



.boo_SLik_joint <- function(simul, debug_level=0L, boo, stat.obs= attr(old_object$logLs,"stat.obs"), 
                            nbCluster_SVs, using=old_object$using, old_object,
                            constr_crits=old_object$constr_crits) {
  boo[] <- simul
  if (debug_level<2L) { # return valid, NULL or try-error
    densv <- suppressWarnings( ## suppress warnings for clustering failure
      try(infer_SLik_joint(boo,stat.obs=stat.obs, nbCluster=nbCluster_SVs, using=using,
                           verbose=list(most=FALSE,final=FALSE, constr_crits=constr_crits,
                                        latentVars=old_object$colTypes$latentVars)),silent=TRUE))
    if (inherits(densv,"try-error") && debug_level==0L) {
      return(NULL) ## used in which valid <- which( ! sapply(resu,is.null)) below
    } 
  } else { # return error -> useful in serial mode only
    densv <- suppressWarnings( ## suppress warnings for clustering failure
      infer_SLik_joint(boo,stat.obs=stat.obs, nbCluster=nbCluster_SVs, using=using,
                       verbose=list(most=FALSE,final=FALSE, constr_crits=constr_crits,
                                    latentVars=old_object$colTypes$latentVars)))
  }
  densv
}


# Called by .RMSEwrapper.SLik_j();
# this does not call the process-simulating function. Instead, it performs a sort of parametric bootstrap from the Gaussian mixture model,
# i.e., it jointly simulates parameters and projected statistics from the jointdens...
# An artefact is that it can draw meaningless values from the jointdens, e.g. negative values for variance parameters
.boot.SLik_j <- function(object,boot_nsim=2L, force=FALSE, verbose = TRUE, seed=NULL, 
                         parent_cores_info=NULL, # defined cluster by parent *.boot.SLik_j* calls
                         packages=get_from(object,"packages"), env=get_from(object,"env"), # not required if parent_cores_info is provided.
                         cluster_args=list(),
                         cl_seed=.update_seed(object)
                         ) {
  if (!force && boot_nsim<2L) stop("'boot_nsim' must be > 1")
  simuland <- object$completedens
  if (is.null(simuland)) simuland <- object$jointdens
  if (inherits(simuland,"dMixmod")) { # using="Rmixmod" or "devel"
    bootrepls <- .simulate.MixmodResults(simuland, n_tables=boot_nsim, size=nrow(object$logLs),
                                         drop=FALSE,
                                         norm_or_t=.wrap_rmvnorm) # simulates 10 (projected) reftables in parametric bootstrap spirit.
    nbCluster_SVs <- list(jointdens=object$jointdens@nbCluster,
                          pardens=object$pardens@nbCluster) # single values so that this exposes to clustering failure.
  } else if (inherits(simuland,"dMclust")) {
    bootrepls <- replicate(boot_nsim, do.call("sim", c(simuland[c("modelName", "parameters")],
                                                       list(n=nrow(object$logLs))))[,-1L],simplify=FALSE)
    nbCluster_SVs <- list(jointdens=object$jointdens$G,
                          pardens=object$pardens$G)
  } else if (inherits(simuland,"MAF")) { # using="mafR"
    bootrepls <- replicate(boot_nsim, expr = {
      .simulate.MAF(simuland, nsim= nrow(object$logLs), given=NULL) # ____F I X M E____ incomplete call ?
    },
    simplify=FALSE)
    nbCluster_SVs <- list()
  } else if (inherits(object$gllimobj,"gllim")) {
    bootrepls <- .simulate.gllimX(object$gllimobj, n_tables=boot_nsim, size=nrow(object$logLs),
                                  parNames=object$colTypes$fittedPars) 
    nclu <- length(object$gllimobj$pi)
    nbCluster_SVs <- list(jointdens=nclu, pardens=nclu)
    
    bootrepls <- .gllim.condsimul.stats(object$gllimobj, RGPpars=bootrepls, size=1L, drop=FALSE, cbind.=TRUE, 
                                 colTypes=object$colTypes) 
  } 
  
  boo <- object$logLs[,with(object$colTypes,c(inferredVars,statNames))]
  attr(boo, "allPars") <- object$colTypes$fittedPars # not the original $allPars as otherwise .boo_SLik_joint -> infer_SLik_joint() tries to use them, including fixedPars,
                                                     # while the bootstrap replicate does not include fixedPars.
  prevmsglength <- 0L
  resu <- vector("list")
  #
  if (is.null(parent_cores_info)) {
    cluster_args <- .set_cluster_type(cluster_args, nb_cores=NULL) # PSOCK vs FORK
    cores_info <- .init_cores(cluster_args=cluster_args)
    cl <- cores_info$cl
    if ( ! is.null(cl)) {
      parallel::clusterExport(cl, "packages",envir=environment()) ## passes the list of packages to load
      # Infusion not loaded in child process !
      #parallel::clusterExport(cl, list("nRealizations"),envir=environment()) 
      #parallel::clusterCall(cl,fun=Infusion.options,nRealizations=nRealizations)
      if ( ! is.null(env)) parallel::clusterExport(cl=cl, ls(env),envir=env)
    }
  } else cores_info <- parent_cores_info
  cl <- cores_info$cl
  if ( ! is.null(cl) && ! is.null(cl_seed) ) parallel::clusterSetRNGStream(cl = cl, cl_seed)
  resu <- pblapply(bootrepls, .boo_SLik_joint, cl = cl, boo=boo, #stat.obs=stat.obs, 
                   nbCluster_SVs=nbCluster_SVs, old_object=object)
  whichvalid <- which( ! sapply(resu,is.null))
  resu <- resu[whichvalid]
  if (length(resu)==0L) {
    message("All bootstrap replicates failed; this suggests a problem that cannot be solved by computing more replicates.\n Trying to diagnose the problem...")
    if (cores_info$nb_cores>1L) {
      message("\nIn parallel mode, first replicate gives:")
      resu <- pblapply(bootrepls[1],.boo_SLik_joint,cl = cores_info$cl, debug_level=1, boo=boo, 
                       #stat.obs=stat.obs, 
                       nbCluster_SVs=nbCluster_SVs, old_object=object)
      print(resu[[1]][1])
      if (identical(cluster_args$debug_info,TRUE)) { # ready-to-use debug code:
        sessioninfo <- utils::sessionInfo()
        infoptions <- Infusion.options()
        using <- object$using
        globoptions <- options()
        tmpname <-  .generateFileName("debug_info",ext=".rda")
        print(paste0("Saving info in file ", tmpname,":"))
        save(bootrepls, sessioninfo, infoptions, globoptions, object, nbCluster_SVs, using, verbose, prevmsglength,
             file=tmpname)
        # loading the .rda allows to run .boo_SLik_joint(bootrepls[[1]]) (and controlling the options, if needed)
      } else message(paste("\nSet cluster_args$debug_info=TRUE to save some debugging info in a file\n",
                           "and see the source of Infusion:::.boot.SLik_j() for further details about it.")) 
      message("\nTesting whether this generates an error in serial mode:") # might not fail if pb only in parallel mode
    } else  message("\nIn serial mode, first replicate gives:") # should fail again
    abyss <- lapply(bootrepls[1],.boo_SLik_joint, debug_level=2, boo=boo, #stat.obs=stat.obs, 
                    nbCluster_SVs=nbCluster_SVs, old_object=object)
    return(NULL) # Only if parallel failed and serial did not. 
  }
  while (length(resu)< boot_nsim) {
    message(paste("Mixture modelling with given nbCluster failed for",boot_nsim-length(resu),"replicate(s); drawing sample(s) again..."))
    moreresu <- .boot.SLik_j(object, boot_nsim=boot_nsim-length(resu), force=TRUE, verbose=verbose, parent_cores_info=cores_info) # recursive call uses cores_info but not cluster_args.
    resu <- c(resu,moreresu)
  }
  if (is.null(parent_cores_info)) .close_cores(cores_info)
  invisible(resu)
}

# Called by MSL()
.RMSEwrapper.SLik_j <- function(object, CIpoints=object$CIobject$bounds, useCI=TRUE, 
                                boot_nsim,
                                verbose=interactive(),
                                cluster_args=list(), level) {
  if (preds_is_mat <- ( useCI && ! is.null(CIpoints) )) {
    locdata <- data.frame(rbind(MSLE=object$MSL$MSLE,CIpoints))
  } else {
    if (is.null(CIpoints)) {
      message("No CI points available => RMSEs limited to the variance of logLik at the MSLE.", immediate. = TRUE)
    } else message("CI points available, but useCI set to FALSE => RMSEs limited to the variance of logLik at the MSLE.", immediate. = TRUE)
    locdata <- object$MSL$MSLE # vector
  }
  #
  bootrepls <- .boot.SLik_j(object,boot_nsim=boot_nsim,verbose=verbose,cluster_args=cluster_args) # returns gaussian mixture models for each resample; no MSL as not needed
  pred_data <- object$logLs[,object$colTypes$fittedPars, drop=FALSE]
  logLpreds <- vector("list", length(bootrepls)) 
  for (it in seq_along(bootrepls)) {
    bootrepl <- bootrepls[[it]]
    bootrepl$thr_info <- .calc_pardens_thr_info(bootrepl, level=level, 
                                                pred_data=pred_data)
    logLpreds[[it]] <- predict(bootrepl,newdata=locdata, which="safe", constr_tuning=FALSE)
  }
  if (preds_is_mat) {
    covmat <- cov(do.call(rbind, logLpreds)) 
    MSEs <- c(MSL=covmat[1,1],diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  } else MSEs <- c(MSL=var( unlist(logLpreds)))
  RMSEs <- sqrt(MSEs)
  if (length(MSEs)>1L) {
    headline <-  "* RMSEs of summary log-L maximum and of its ratio at CI bounds: *\n"
  } else {
    headline <-  "* RMSEs of summary log-L maximum: *\n"
  }
  if(verbose) {
    cat(headline)
    print(RMSEs)
  }
  return(RMSEs)
}

# overcome problems with plot Rmixmod which may need more than the info i na best result
# chercher 
# setMethod(
#   f="plot",
#   signature=c("Mixmod"),
# dans les sources pour savoir ce que Rmixmod fait.
# cf also plotCluster(lik_j$jointdens,data=slik_j$logLs)
# .plotClusters <- function(object,which,...) {
#   plotCluster(object[[which]],data=object$logLs,...)
#}

plot.MixmodResults <- function(x,...) {message('Use plotCluster(.,data=<slik object>$logLs,variable1=.,variable2=.)\n to plot a "MixModResults" object')}
