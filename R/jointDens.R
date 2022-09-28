.check_data_post_boundaries <- function(data) {
  data <- na.omit(data)
  if (nNAlines <- length(attr(data,"na.action"))) {
    message(paste0(nNAlines," lines of reference table still contained NA's after handling 'boundaries', and were removed. see ?handling_NAs for details."))
  }
  checkDegenerate <- cov(data)
  #if ( ! anyNA(checkDegenerate)) { # cf  use="na.or.complete" => NA if no complete cases
  eigvals <- eigen(checkDegenerate, only.values = TRUE)$values
  if (any(abs(eigvals<1e-14))) warning(paste(
    c("The covariance matrix of the (parameters,statistics) table seems singular,\n",
      "implying linear relationships between the variables. Problems will likely happen.\n",
      "Redundant variables should be eliminated."
    )),immediate. = TRUE)
  #}
  data
}


## fast version of dmvnorm that uses precomputed version of chol(sigma) and minimal checking
.fast_dmvnorm <- function (x, mean, solve_t_chol_sigma, log = FALSE) {
  if (is.vector(x)) dim(x) <- c(1,length(x))
  if (nrow(x)==1L) { 
    tmp <- tcrossprod(solve_t_chol_sigma, x-mean)
    rss <- sum(tmp^2)
  } else {
    tmp <- solve_t_chol_sigma %*% (t(x)-mean) # backsolve(chol_sigma, t(x) - mean, transpose = TRUE)
    rss <- .colSums(tmp^2,m=nrow(tmp), n=ncol(tmp))
  }
  logdet <- attr(solve_t_chol_sigma,"logdet")
  if (is.null(logdet)) logdet <- sum(log(diag(solve_t_chol_sigma))) # is (2017) method, precomputation is not interesting 
  logretval <- logdet - 0.5 * ncol(x) * log(2 * pi) - 0.5 * rss
  names(logretval) <- rownames(x)
  if (log) { logretval } else exp(logretval)
}

.marginalize_Rmixmod <- function(jointdens, colNames, For, over) {
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

.conditional_Rmixmod <- function(jointdens, #fittedPars, 
                                 given, expansion=Infusion.getOption("expansion")) {  # expansion=1 to get the conditional distribution
  nbCluster <- jointdens@nbCluster
  conddens <- jointdens
  MEAN <- conddens@parameters@mean
  givenNames <- names(given)
  colnames(MEAN) <- colNames <- jointdens@varNames
  For <- setdiff(colNames,givenNames) 
  conddens@parameters@mean <- MEAN[,For,drop=FALSE] # resizes, but will be modified
  condprop <- conddens@parameters@proportions
  for (clu_it in seq_len(nbCluster)) {
    COV <- conddens@parameters@variance[[clu_it]]
    colnames(COV) <- rownames(COV) <- colNames
    sig22 <-  COV[givenNames,givenNames,drop=FALSE]
    sig12 <-  COV[For,givenNames,drop=FALSE]
    mean2 <- MEAN[clu_it,givenNames]
    conddens@parameters@mean[clu_it,] <- MEAN[clu_it,For] + sig12 %*% solve(sig22,given-mean2)
    sig11 <- COV[For,For,drop=FALSE]
    conddens@parameters@variance[[clu_it]] <- expansion* (sig11 - sig12 %*% solve(sig22,t(sig12))) 
    condprop[clu_it] <- log(condprop[clu_it])+dmvnorm(t(given), # dmvnorm() tests is.vector(x) which returns FALSE if x has attributes other than names.
                                                 mean2, sigma= sig22, log=TRUE)
  }
  maxlog <- max(condprop)
  condprop <- condprop-maxlog
  condprop <- exp(condprop)
  condprop <- condprop/sum(condprop)
  conddens@parameters@proportions <- condprop
  conddens@varNames <- For
  conddens@simuls_activeBoundaries <- NULL
  #margdens_fP@criterionValue <- margdens_fP@likelihood <- "removed for safety"  
  return(conddens)
}

.conditional_mclust <- function(jointdens, fittedPars, given, 
                                expansion=Infusion.getOption("expansion")) { # expansion=1 to get the conditional distribution
  nbCluster <- jointdens$G
  conddens <- jointdens
  For <- fittedPars 
  conddens$data <- conddens$data[,For]
  conddens$parameters$variance$d <- ncol(conddens$data)
  MEAN <- conddens$parameters$mean
  givenNames <- names(given)
  colNames <- c(fittedPars, givenNames)
  conddens$parameters$mean <- MEAN[For,,drop=FALSE] # resizes, but will be modified
  sigma11 <- conddens$parameters$variance$sigma[For,For,,drop=FALSE] # resizes, but will be modified
  condprop <- conddens$parameters$pro
  for (clu_it in seq_len(nbCluster)) {
    sigma_it <- conddens$parameters$variance$sigma[,,clu_it] # from single array for all clusters # drops the clu_it dimension
    sig22 <-  sigma_it[givenNames,givenNames]
    sig12 <-  sigma_it[For,givenNames]
    mean2 <- MEAN[givenNames,clu_it]
    conddens$parameters$mean[,clu_it] <- MEAN[For,clu_it] + sig12 %*% solve(sig22,given-mean2)
    sigma11[,,clu_it] <- expansion* (sigma_it[For,For] - sig12 %*% solve(sig22,t(sig12))) 
    condprop[clu_it] <- condprop[clu_it]*dmvnorm(t(given), # dmvnorm() tests is.vector(x) which returns FALSE if x has attributes other than names.
                                                 mean2, sigma= sig22, log=FALSE)
  }
  maxlog <- max(condprop)
  condprop <- condprop-maxlog
  condprop <- exp(condprop)
  condprop <- condprop/sum(condprop) # normalisation (even in non-safe-log-exp version) was missing for a long time!
  conddens$parameters$pro <- condprop
  conddens$parameters$variance <- .do_call_wrap("sigma2decomp", 
                                                list(sigma=sigma11,tol=-1), # tol to avoid a simplification in the structure of $orientation that may hit a but in mclust...
                                                pack="mclust") 
  # sigma2decomp determines it is of model "VVV". We must copy this info so that the right simulation function is called by Mclust:
  conddens$modelName <- conddens$parameters$variance$modelName
  return(conddens)
}


.marginalize_mclust <- function(jointdens, colNames, For, over) {
  nbCluster <- jointdens$G
  margdens <- jointdens
  margdens$data <- margdens$data[,For, drop=FALSE]
  margdens$parameters$variance$d <- ncol(margdens$data)
  margdens$parameters$mean <- margdens$parameters$mean[For,,drop=FALSE]
  COV <- margdens$parameters$variance$sigma # single array for all clusters 
  COV <- COV[For,For,]
  if (length(For)==1L) {
    if (length(COV)==1L) {
      margdens$parameters$variance <- list(modelName="X", d=1, G=1L, sigmasq=COV)
    } else {
      # "V" may not be implied by the original fit, but we don't count dfs on this object
      margdens$parameters$variance <- list(modelName="V", d=1, G=length(COV), sigmasq=COV, scale=COV) 
    }
  } else {
    margdens$parameters$variance <- .do_call_wrap("sigma2decomp", 
                                                  list(sigma=COV,tol=-1), # tol to avoid a simplification in the structure of $orientation that may hit a but in mclust...
                                                  pack="mclust") 
    if (margdens$modelName=="VVV") {
      cholsigma <- array(0,dim=dim(COV))
      for (it in seq_len(dim(cholsigma)[3L])) cholsigma[,,it] <- chol(COV[,,it])
      margdens$parameters$variance$cholsigma <- cholsigma
    }
  }
  # Not necess useful here, but following thesamelogic as in .conditional_mclust():
  margdens$modelName <- margdens$parameters$variance$modelName
  return(margdens)
}

.get_best_mixmod_by_IC <- function(cluObject, which=Infusion.getOption("criterion"))  {
  if (inherits(cluObject,"try-error")) return(cluObject) ## passes original error info rather than hiding it under another error
  results <- cluObject@results
  if (length(results)==1L) return(results[[1L]])
  anyNaN <- logical(length(results))
  # patch for odd bug of Rmixmod: NaN in parameters but no error reported, and lik & BIC are real.
  for (it in seq_along(results)) anyNaN[it] <- anyNA(results[[it]]@proba)
  results <- results[ ! anyNaN]
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

.solve_t_cholfn <- function(mat, smoothing_mat=NULL) { # smoothing_mat arg added to allow optimization via iterative smmothing (which does not work)
  if (!is.null(smoothing_mat)) mat <- mat+smoothing_mat
  resu <- solve(t(chol(mat)))
  return(structure(resu,logdet=sum(log(diag(resu)))))
}

seq_nbCluster <- function(nr) {do.call(Infusion.getOption("seq_nbCluster"), list(nr=nr))}

get_nbCluster_range <- function(projdata, nr= nrow(projdata), nc=ncol(projdata), 
                                nbCluster=seq_nbCluster(nr)) { # data must be the mixmod data
  maxnb <- do.call(.Infusion.data$options$maxnbCluster, list(nr=nr,nc=nc))
  if (maxnb<1L) {
    warning("Too few simulations: Gaussian mixture modelling may fail; single cluster inferred", immediate. = TRUE)
    return(1L)
  } else if (any(nbCluster> maxnb)) {
    message(paste0("Gaussian mixture modelling constrained to at most ",maxnb," clusters by the number of columns."))
    return(nbCluster[nbCluster<=maxnb])
  } else nbCluster 
}

infer_SLik_joint <- function(data, ## reference table ~ abc
                            stat.obs,
                            logLname=Infusion.getOption("logLname"), ## not immed useful
                            Simulate=attr(data,"Simulate"), ## may be NULL
                            nbCluster=seq_nbCluster(nr=nrow(data)),
                            using=Infusion.getOption("mixturing"),
                            verbose=list(most=interactive(), ## must be explicitly set to FALSE in knitr examples
                                         pedantic=FALSE,
                                         final=FALSE),
                            marginalize=TRUE
) {
  if ( ! is.data.frame(data)) {
    stop(paste("'object' is not a data.frame.\n Did you mean to call infer_logLs() rather than infer_Slik_joint() ?"))
  }
  if (is.null(raw_data <- attr(data,"raw_data"))) {
    if (is.null(attr(data,"LOWER"))) { # Temporary patch? stop() instead of warning() ?
      warning('Some required attributes seem to be missing from the "data" [see "Value" in help("add_reftable")].\n Further execution could fail.')
    }
  } else if (is.null(attr(raw_data,"LOWER"))) { # Temporary patch? stop() instead of warning() ?
    warning('Some required attributes seem to be missing from attr(data,"raw_data" [see "Value" in help("add_reftable")].\n Further execution could fail.')
  }
  if ( ! is.null( cn <- colnames(stat.obs))) {
    message("Note: 'stat.obs' should be a numeric vector, not a matrix or data.frame. Converting...")
    stat.obs <- drop(stat.obs)
    names(stat.obs) <- cn ## minimal patch so that names() can be used, not colnames()
  }
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- c("most","final")[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- interactive()
  if (is.null(verbose$pedantic)) verbose$pedantic <- FALSE
  if (is.null(verbose$final)) verbose$final <- FALSE
  if (length(unique(colnames(data)))<ncol(data)) {
    stop("Some columns of the 'data' have identical names. Please correct.")
  }
  statNames <- names(stat.obs)
  if (nrow(data)<20000) check_raw_stats(na.omit(data),statNames = statNames,remove = FALSE, verbose=FALSE) # silent if not problem detected, but still verbose otherwise
  for (st in statNames) {
    tabstat <- table(data[,st])
    #which_mass <- as.numeric(names(which(tabstat>1L)))
    if (any(tabstat>Infusion.getOption("nodesize"))) { # this is a projector-specific value yet we use a single value fro the slik object
      # There is <ranger object>$min.node.size but no such info in the RandomForest results.
      massvalues <- as.numeric(names(which.max(tabstat))) 
      names(massvalues) <- rep(st,length(massvalues))
      attr(stat.obs,"boundaries") <- c(attr(stat.obs,"boundaries"),massvalues)
    }
  }
  if (is.null(allPars <- attr(data,"allPars"))) {  
    allPars <- setdiff(colnames(data),statNames) ## first iteration
  } ## else do not add the cumul_iter variable to the parameters !
  whichvar <- apply(data[,allPars,drop=FALSE],2L,function(v) length(unique(range(v)))>1L)
  fittedPars <- allPars[whichvar]
  fixedPars <- data[1L,!whichvar,drop=FALSE] ## has both names and values
  rownames(fixedPars) <- NULL## avoids warning on rownames when later used in cbind()
  if (! is.list(nbCluster)) nbCluster <- list(jointdens=nbCluster, pardens=nbCluster)
  nbCluster <- lapply(nbCluster,eval.parent)
  if (using=="mclust") {
    if ( ! "package:mclust" %in% search()) stop("'mclust' should be loaded first.")
    if (verbose$most) cat(paste("Densities modeling:",nrow(data),"points;\n"))
    
    if (length(nbCluster$jointdens)==1L) {
      jointdens <- .do_call_wrap("densityMclust",
                                 list(data=data[,c(fittedPars,statNames)],modelNames=Infusion.getOption("mclustModel"), 
                                      G=nbCluster$jointdens,verbose=FALSE),
                                 pack="mclust")
    } else { 
      models <- vector("list",length(nbCluster$jointdens))
      for (it in seq_along(nbCluster$jointdens)) {
        models[[it]] <- .do_call_wrap("densityMclust",
                                      list(data=data[,c(fittedPars,statNames)],modelNames=Infusion.getOption("mclustModel"), 
                                           G=nbCluster$jointdens[it],verbose=FALSE),
                                      pack="mclust")
      }
      jointdens <- .get_best_mclust_by_IC(models) 
    }
    if (verbose$most) cat(paste(" joint density modeling:",jointdens$G,"clusters;\n"))
    if (verbose$pedantic) if (jointdens$G==max(nbCluster$jointdens)) message("Inferred # of clusters= max of allowed range.")
    if (marginalize) {
      pardens <- .marginalize_mclust(jointdens, colNames= c(fittedPars,statNames),
                              For=fittedPars, over=statNames)
      # defective but should be sufficient for prediction. Some garde-fou:
      #pardens@proba <- matrix(NA_real_)
    } else {
      if (length(nbCluster$pardens)==1L) {
        pardens <- .do_call_wrap("densityMclust",
                                 list(data=data[,fittedPars,drop=FALSE],modelNames=Infusion.getOption("mclustModel"), G=nbCluster$pardens,verbose=FALSE),
                                 pack="mclust")
      } else { 
        models <- vector("list",length(nbCluster$pardens))
        for (it in seq_along(nbCluster$pardens)) {
          models[[it]] <- .do_call_wrap("densityMclust",
                                        list(data=data[,fittedPars,drop=FALSE],modelNames=Infusion.getOption("mclustModel"), 
                                             G=nbCluster$pardens[it], verbose=FALSE),
                                        pack="mclust")
        }
        pardens <- .get_best_mclust_by_IC(models) 
      }
      if (verbose$most) cat(paste(" parameter modeling:",pardens$G,"clusters.\n"))
    }
    resu <- list(jointdens=jointdens,pardens=pardens)
    # pdl <- vector("list", pardens$G)
    # sigma <- resu$pardens$parameters$variance$sigma
    # for (it in seq_len(pardens$G)) pdl[[it]] <- solve(t(chol(sigma[,,it])))
    # jdl <- vector("list", jointdens$G)
    # sigma <- resu$jointdens$parameters$variance$sigma
    # for (it in seq_len(jointdens$G)) jdl[[it]] <- solve(t(chol(sigma[,,it])))
    # resu$solve_t_chol_sigma_lists <- list(pardens=pdl, jointdens=jdl) 
  } else {
    models <- .do_call_wrap("mixmodGaussianModel",list(listModels=Infusion.getOption("mixmodGaussianModel")))
    if (verbose$most) cat(paste("Densities modeling:",nrow(data),"points;\n"))
    locarglist <- list(data=data[,c(fittedPars,statNames)],nbCluster=nbCluster$jointdens, 
                       seed=Infusion.getOption("mixmodSeed") , models=models)
    if (FALSE) { ## older code; still useful for devel
      jointdens <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
      jointdens <- .get_best_mixmod_by_IC(jointdens)
    }
    jointdens <- try(do.call(".densityMixmod",c(locarglist,list(stat.obs=stat.obs))),silent=TRUE) # using seed in locarglist ie that from geoOption
    while ((inherits(jointdens,"try-error") || anyNA(jointdens@parameters@mean)) && max(locarglist$nbCluster)>1L) {
      # Generate better (smaller) values with two cases whether a range or a single value were tried.
      if (length(locarglist$nbCluster)==1L) { # typical of RMSE bootstrap sample. Reduce by 1 or more
        nbCluster <- min(locarglist$nbCluster-1L, max(get_nbCluster_range(projdata=data, nbCluster=seq(locarglist$nbCluster))))
      } else  nbCluster <- get_nbCluster_range(projdata=data, nbCluster=seq(max(locarglist$nbCluster))) # generate a more standard range
      if ( length(setdiff(nbCluster, locarglist$nbCluster))) { # if there are new values to try
        locarglist$nbCluster <- nbCluster
        jointdens <- try(do.call(".densityMixmod",c(locarglist,list(stat.obs=stat.obs))),silent=TRUE) # using seed in locarglist ie that from geoOption
      } else break # = no more ideas to fix problem
    } ## else .densityMixmod has already run .get_best_mixmod_by_IC()
    if (inherits(jointdens,"try-error"))  {
      stop(jointdens) # i.e stop(<error object>) # (But in RMSE bootstrap, new bootstrap samples may be analyzed if a few failed) 
    } else if (anyNA(jointdens@parameters@mean)) stop("Rmixmod returned NaN's in parameter estimates.")
    # plotCluster(jointdens,data=locarglist$data,variable1="theta_p",variable2="theta") # var1: stat (prediction of projection); var2: actual param
    if (verbose$most) cat(paste(" joint density modeling:",jointdens@nbCluster,"clusters;\n"))
    if (verbose$pedantic) if (jointdens@nbCluster==max(nbCluster$jointdens)) message("Inferred # of clusters= max of allowed range.")
    if (marginalize) {
      pardens <- .marginalize_Rmixmod(jointdens, colNames= c(fittedPars,statNames),
                              For=fittedPars, over=statNames)
      # defective but should be sufficient for prediction. Some garde-fou:
      pardens@proba <- matrix(NA_real_)
      pardens@parameters@nbFreeParam <- NA_integer_
    } else {
      locarglist <- list(data=data[,fittedPars,drop=FALSE], nbCluster=nbCluster$pardens, 
                         seed=123 , models=models)
      pardens <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
      if (inherits(pardens,"try-error")) {
        # ignore any non-default infer_SLik_joint(.,nbCluster) argument
        locarglist$nbCluster <- get_nbCluster_range(projdata=data)
        pardens <- do.call("mixmodCluster",locarglist)
      }                               
      pardens <- .get_best_mixmod_by_IC(pardens) 
      if (verbose$most) cat(paste(" parameter modeling:",pardens@nbCluster,"clusters.\n"))
    }
    #plotCluster(pardens,data=data[,fittedPars]) to plot @results[[1L]] which is a 'MixModResults', not a 'mixmodCluster' object.
    # Rmixmod::plot(<mixmodCluster object>) with Rmixmod:: to plot from any envir, not only the global one
    resu <- list(jointdens=jointdens,pardens=pardens)
    resu$solve_t_chol_sigma_lists <- list(
      pardens=lapply(resu$pardens@parameters["variance"], .solve_t_cholfn),
      jointdens=lapply(resu$jointdens@parameters["variance"], .solve_t_cholfn)
    ) 
  }
  attr(resu,"EDFstat") <- "[see this string in infer_SLik_joint()]" ## 
  resu$logLs <- structure(data,stat.obs=stat.obs,Simulate=Simulate) ## as in infer_surface.logLs
  # attr(resu,"callfn") <- "infer_SLik_joint"
  resu$projectors <- attr(data,"projectors")
  resu$`Infusion.version` <- packageVersion("Infusion")
  resu$colTypes <- list(allPars=allPars, ## keeps the order of the columns
                                fittedPars=fittedPars,
                                fixedPars=fixedPars,
                                logLname=logLname,
                                statNames=statNames)
  resu$lower <- apply(data[,fittedPars,drop=FALSE],2L,min) ## really this
  resu$upper <- apply(data[,fittedPars,drop=FALSE],2L,max) # ... used by MSL -> optim
  if (any(is.na(c(resu$lower,resu$upper)))) stop("NA(s) in c(lower,upper))") 
  resu$LOWER <- attr(data,"LOWER") # used in many places
  if (is.null(resu$LOWER)) resu$LOWER <- resu$lower
  resu$UPPER <- attr(data,"UPPER") # ... used by confintAll
  if (is.null(resu$UPPER)) resu$UPPER <- resu$upper
  attr(resu,"Infusion.version") <- packageVersion("Infusion")
  class(resu) <- c("SLik_j",class(resu))
  resu$raw_data <- attr(data,"raw_data")
  return(resu)
}


## was infer_SLik_joint -> predict -> predict.SLik_j -> predict.MixmodResults
# except that now                                    -> .pointpredict.Rmixmod -> predict.dMixmod
predict.MixmodResults <- function(object, newdata,log=TRUE, solve_t_chol_sigma_list,...) {
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
        density[,k] <- log(object@parameters["proportions", k]) + 
          .fast_dmvnorm(newdata, object@parameters["mean", k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=log)
      }
      maxlogs <- apply(density,1,max)
      normalizedlogs <- apply(density,2L,`-`,maxlogs) ## highest value per row is 0
      ## apply return mess again: 
      if (is.vector(normalizedlogs)) dim(normalizedlogs) <- c(1L, length(normalizedlogs))
      mixture <- .rowSums(exp(normalizedlogs),m=nrow(normalizedlogs), n=ncol(normalizedlogs)) ## exp(normalizedlogs) in (0,1); sum(exp(logLi-maxlog))= exp(-maxlog) sum(exp(logLi))= exp(-maxlog) sum(Li)
      mixture <- log(mixture) + maxlogs ## log(sum(Li))
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

.pointpredict.Rmixmod <- function(X, # parameters only 
                                  object, 
                                  tstat.obs, # 1-row matrix as otherwise more cases should be considered for cbind'ing
                                  log, 
                                  which, # "safe" version ignores, by correcting, spuriously high logL in area of low parameter density.
                                  solve_t_chol_sigma_lists=object$solve_t_chol_sigma_lists, # added to allow optimization via iterative smmothing (which does not work)
                                  ...)  {
  if (is.null(dim(X))) {
    newjointX <- cbind(t(X),tstat.obs) # cbind two 1-row matrices
  } else {
    if (length(intersect(colnames(X),colnames(tstat.obs)))) stop("'X' should contain only parameters, not summary statistics")
    newjointX <- cbind(X,tstat.obs[rep(1,nrow(X)),,drop=FALSE]) # cbind two nrow(X)-row matrices
  }
  if (which!="parvaldens") {
    jointvaldens <- predict(object$jointdens,
                            newdata=newjointX,
                            solve_t_chol_sigma_list=solve_t_chol_sigma_lists$jointdens,
                            log=log,...)
  }
  if (which=="jointvaldens") return(jointvaldens)
  parvaldens <- predict(object$pardens,
                        newdata=newjointX, # statistics will be ignored
                        solve_t_chol_sigma_list=solve_t_chol_sigma_lists$pardens,
                        log=log,...) 
  if (which=="parvaldens") return(parvaldens)
  # ELSE: "lik", or "safe" for safe version of "lik" using object$thr_dpar
  if (log) {
    condvaldens <- jointvaldens - parvaldens
    if (which=="safe") condvaldens <- condvaldens + pmin(0,parvaldens-object$thr_dpar) # decreasing returned logL when parvaldens<object$thr_dpar
  } else {
    condvaldens <- jointvaldens/parvaldens
    if (which=="safe") condvaldens <- condvaldens*pmin(1,parvaldens/object$thr_dpar)
  }
  return(condvaldens) # vector if X was a matrix
}


.predict_SLik_j_mclust <- function(object, newdata, tstat.obs, log, which, solve_t_chol_sigma_lists=object$solve_t_chol_sigma_lists)  {
  if (which!="parvaldens") jointvaldens <- predict(object$jointdens, newdata=cbind(data.frame(newdata),tstat.obs), 
                                                   solve_t_chol_sigma_lists=solve_t_chol_sigma_lists) # ! order of elements in newdata must be that of fittedPars as in object$jointdens$data
  if (which=="jointvaldens") return(jointvaldens)
  parvaldens <- predict(object$pardens,newdata=newdata, 
                        solve_t_chol_sigma_lists=solve_t_chol_sigma_lists) 
  if (which=="parvaldens") return(parvaldens)
  if (log) {
    jointvaldens <- pmax(jointvaldens, .Machine$double.xmin) # F I X M E for the value
    condvaldens <- log(jointvaldens) - log(parvaldens)
  } else condvaldens <- jointvaldens/parvaldens
  if (any(is.infinite(condvaldens))) warning("any(is.infinite(condvaldens)) is TRUE: expect an error.") 
  return(condvaldens)
}

predict.SLik_j <- function (object, 
                            newdata, ## requests new fittedPars values! 
                            log=TRUE, 
                            which="lik", # may be still preferable, to "safe", for drawing new points
                            tstat= t(attr(object$logLs,"stat.obs")), # 1-row matrix...
                            solve_t_chol_sigma_lists=object$solve_t_chol_sigma_lists, # added to allow optimization via iterative smmothing (which does not work)
                            ...) {
  if (is.null(nrow(newdata)) ) newdata <- t(as.matrix(newdata)) # as.matrix keeps names
  if (is.null(colnames(newdata))) colnames(newdata) <- object$colTypes$fittedPars
  # not useful bc predict.MixmodResults does not handle it:
  #newdata <- cbind(newdata, t(replicate(nrow(newdata),attr(object$logLs,"stat.obs"))))
  if (inherits(object$jointdens,"densityMclust")) {
    logl <- .predict_SLik_j_mclust(object=object, newdata=newdata, tstat.obs=tstat, log=log, which=which, solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)
  } else {
    logl <- .pointpredict.Rmixmod(X=newdata, object=object, tstat.obs=tstat, log=log, which=which, solve_t_chol_sigma_lists=solve_t_chol_sigma_lists)
    names(logl) <- rownames(newdata) ## predict.glm (say) keeps names; 
                                    ## RMSEs computation uses names of CI bounds, to be retained here, to name its results
    # apply(newdata,1L, .pointpredict.Rmixmod, object=object, tstat.obs=tstat, log=log, which=which)
  }
  if ( ! is.null(object$prior_logL)) logl <- logl + object$prior_logL(newdata)
  return(logl)
}

confint.SLik_j <- function(object, parm, ## parm is the parameter which CI is sought 
                         level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE), 
                         #Bartlett=0L, 
                         ...) {
  .confintAll(object=object, parm=parm, ## parm is the parameter which CI is sought 
             givenmax = object$MSL$maxlogL,
             level= - qchisq(level,df=1)/2, ## df=1 for 1D profile; /2 to compare to logLR rather than 2 logLR  
             #Bartlett=Bartlett,
             verbose=verbose,fixed=fixed,which=which,...)
}

refine.SLik_j <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "mixmodCluster" ## no clear effect but distinct from SLik case
  resu <- refine.default(object, surfaceData=object$logLs, method=method, ...)
  assign("bootreps_list", list(), envir=resu$bootLRTenv) # or reset the env by resu$bootLRTenv <- list2env(list(bootreps_list=list()))?
  resu
}

summary.SLik_j <- summary.SLik

print.SLik_j <- function(x,...) {summary.SLik_j(x,...)} 

.rparam_dMixmod_around_focal <- function(object, # *dMixmod* object
                                 focal, ## single point # recall that there exists .rparam_from_SLik_j_post(), 
                                        # which is similar is 'object' is the posterior density, EXCEPT that it does not handle this additional constraint
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
  fullfocal <- c(focal,prof$par)[fittedPars]
  posteriordens <- .conditional_Rmixmod(object$jointdens, given=attr(object$logLs,"stat.obs"), expansion=1) 
  solve_t_chol_sigma_list <- lapply(posteriordens@parameters["variance"], .solve_t_cholfn)
  trypoints <- .rparam_dMixmod_around_focal(posteriordens, focal = fullfocal, solve_t_chol_sigma_list=solve_t_chol_sigma_list, size=size)
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

.rparam_from_SLik_j <- function(object, 
                                frac, 
                                target_size=NULL,
                                fittedPars=NULL,level,tol=0) {
  if (is.null(fittedPars)) fittedPars <- object$colTypes$fittedPars
  CI_LRstat <- qchisq(level,df=length(fittedPars))/2 
  size_first_iter <- attr(object$logLs,"n_first_iter")
  if (is.null(size_first_iter)) size_first_iter <- length(which(object$logLs$cumul_iter==1)) ## FR->FR temporary back compat
  prev_n_iter <- max(object$logLs$cumul_iter)
  if (is.null(target_size)) target_size <- frac*max(size_first_iter/2, size_first_iter*(prev_n_iter+1L)-nrow(object$logLs))
  #
  ceil_size <- target_size
  if ( ! is.null(freq_good <- attr(object$logLs,"freq_good")$default)) ceil_size <- ceiling(ceil_size/freq_good)
  #
  if (inherits(object$jointdens,"densityMclust")) {
    trypoints <- do.call("sim", c(object$jointdens[c("modelName", "parameters")],list(n=ceil_size)))
    trypoints <- trypoints[,1L+seq_len(length(fittedPars)),drop=FALSE]  
  } else if (FALSE && # inhibits block for generation of specific param points when profiling has found a new maximum (__F I X M E__)
             (identical(object$MSL$updated_from_prof,TRUE) || # object$MSL envir recreated by user after plot1Dprof, plot2Dprof or plot() found a new max
              ! is.null(object$MSL$init_from_prof) # only .MSL_update() run -> object$MSL$MSLE modified but envir not recreated, $init_from_prof still there.
             )
            ) {
      # if we sample always in the posterior = prior*lik, we will end sampling according to the lik, so too much concentrated. 
      size_from_pardens <- ceil_size* (1-1/object$jointdens@nbCluster)
      trypoints <- .simulate.MixmodResults(object$jointdens, nsim=1L, size=size_from_pardens, drop=TRUE) ## simulates joint
      trypoints <- trypoints[,seq_len(length(fittedPars)),drop=FALSE] # sample in marginal distrib of params from sample in jointdens
      if (FALSE) {
        postpoints <- .rparam_dMixmod_around_focal(object$jointdens, # dMixmod object
                                                   focal=c(object$MSL$MSLE,attr(object$logLs,"stat.obs")), 
                                                   solve_t_chol_sigma_list=object$solve_t_chol_sigma_lists$jointdens,
                                                   size= ceil_size - nrow(trypoints))
        # this samples the joint density around the focal point so +/- conditionally on stat.obs. 
        # Then next line +/- samples in posterior distrib of params given stat.obs  
        postpoints <- postpoints[,seq_len(length(fittedPars)),drop=FALSE] 
      } else {
        posteriordens <- .conditional_Rmixmod(object$jointdens, given=attr(object$logLs,"stat.obs"), expansion=1) 
        solve_t_chol_sigma_list <- lapply(posteriordens@parameters["variance"], .solve_t_cholfn)
        postpoints <- .rparam_dMixmod_around_focal(posteriordens, # dMixmod object
                                                   focal=object$MSL$MSLE, 
                                                   solve_t_chol_sigma_list=solve_t_chol_sigma_list,
                                                   size= ceil_size - nrow(trypoints))
      }
      trypoints <- rbind(trypoints,postpoints)
    } else {
      trypoints <- .simulate.MixmodResults(object$jointdens, nsim=1L, size=ceil_size, drop=TRUE) ## simulates joint
      trypoints <- trypoints[,seq_len(length(fittedPars)),drop=FALSE] # sample in marginal distrib of params from sample in jointdens
      ## : samples the joint and keep only the params => distrib of trypoints=current distr of param (pardens)
    }
  
  colnames(trypoints) <- fittedPars
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],,drop=FALSE]
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],,drop=FALSE]
  #
  prior <- predict(object,trypoints,which="parvaldens") 
  logLs <- predict(object,trypoints,which="lik") ## likelihood fn  
  flatnd_max <- max(logLs)-CI_LRstat ## defines an upper slice of logLs values.
  upperslice <- ( logLs> flatnd_max )
  weights <- numeric(length(logLs))
  weights[ upperslice] <- max(prior[upperslice])-prior[upperslice] ## in log density units
  max_upper <- max(weights[ upperslice]) ## in log density units
  weights[ ! upperslice] <- logLs[ ! upperslice] -prior[ ! upperslice] 
  max_lower <- max(weights[ ! upperslice])
  weights[ ! upperslice] <- weights[ ! upperslice] + max_upper-max_lower # hence same maximum 'max_upper' in the two regions.
  # Refinement of the same idea but with lesser impact of logLs variation outside the upperslice:
  weights[ ! upperslice] <- ((1-tol)*weights[ ! upperslice]+tol*max_upper) ## still in log density units
  # we could correct  ! upperslice acording to 'prior' but this not be worth the effort. 
  #
  weights <- exp(weights-max(weights)) 
  wei_runif <- runif(n=length(weights))<weights ## allows some exploration by sampling some 'bad' points
  good <- which(wei_runif) ## samples uniformly in the current fitted distrib of params with logL>threshold,
                                                  ## with tapering around
  freq_good <- length(good)/ceil_size 
  if (length(good)>target_size) {
    good <- good[sample(length(good),size = target_size)]
  } else if (length(good)>target_size-10L) { # just missing a few points... 
    bad <- which(!wei_runif)
    order_bad <- order(weights[bad],decreasing=TRUE)
    supplement <- bad[order_bad][seq(min(length(bad),target_size-length(good)))]
    good <- c(good,supplement)
  } else if (freq_good < 1/5) {
    bad <- which(!wei_runif)
    order_bad <- order(weights[bad],decreasing=TRUE)
    supplement <- bad[order_bad][seq(ceil_size/5-length(good))]
    good <- c(good,supplement)
  } # else may stay a few points below target, rather than add points that were not retained bc of low weights (it's better to wait next iteration)
  trypoints <- trypoints[good,,drop=FALSE]
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  return(list(trypoints=trypoints, 
              freq_good=max(1/5,freq_good) # max() to avoid explosion of # of sampled points
             ))
}


.rparam_from_SLik_j_post <- function(object , frac) {
  if (frac==0) return(list(trypoints=NULL, freq_good=NA))
  size_first_iter <- attr(object$logLs,"n_first_iter")
  if (is.null(size_first_iter)) size_first_iter <- length(which(object$logLs$cumul_iter==1)) ## FR->FR temporary back compat
  prev_n_iter <- max(object$logLs$cumul_iter)
  ceil_size <- target_size <- frac*max(size_first_iter/2, size_first_iter*(prev_n_iter+1L)-nrow(object$logLs))
  if ( ! is.null(freq_good <- attr(object$logLs,"freq_good")$posterior)) ceil_size <- ceiling(ceil_size/freq_good)
  fittedPars <- object$colTypes$fittedPars
  if (inherits(object$postdens,"densityMclust")) {
    trypoints <- do.call("sim", c(object$postdens[c("modelName", "parameters")],list(n=ceil_size)))
    trypoints <- trypoints[,1L+seq_len(length(fittedPars)),drop=FALSE]
  } else {
    trypoints <- .simulate.MixmodResults(object$postdens, nsim=1L, size=ceil_size, drop=TRUE)
  } 
  colnames(trypoints) <- fittedPars
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  trypoints <- data.frame(trypoints)
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],,drop=FALSE]
  for (vv in fittedPars) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],,drop=FALSE]
  freq_good <- nrow(trypoints)/ceil_size
  return(list(trypoints=trypoints, freq_good=freq_good))
}




## FR->FR define rparam as generic with methods SLik_j, etc ?
.bad_rparam_from_SLik_j <- function(object, size=1L,CIweight=Infusion.getOption("CIweight")) {
  locpredict <- function(x) {predict.SLik_j(object,newdata=x)}
  MSLE <- object$MSL$MSLE
  fittedNames <- names(MSLE)
  RMSEs <- get_from(object,"RMSEs")
  #
  nvec <- .make_n(RMSEs=RMSEs, fittedPars=fittedNames, n=size, CIweight=CIweight)
  #
  trypoints <- matrix(NA,ncol=length(MSLE),nrow=0L)
  if (nvec["MSL"]>0L) {
    hess <- hessian(locpredict,x=MSLE)
    sigma <- - solve(hess) # *(RMSEs["MSL"])
    rmv <- rmvnorm(nvec["MSL"],mean=MSLE,sigma=sigma) # 1st sampling
    trypoints <- rbind(trypoints,rmv) ## centers points on the MSLE 
  }
  CIpoints <- object$CIobject$bounds
  for (kt in seq_len(NROW(CIpoints))) {
    focalpt <- CIpoints[kt,] 
    locst <- rownames(CIpoints)[kt]
    if (nvec[locst]>0L) {
      dotpos <- regexpr(".",locst,fixed=TRUE)[1L]
      parm <- substring(locst,dotpos+1L)
      parpos <- which(fittedNames==parm)
      curv <- hessian(locpredict,x=focalpt)
      curv[parpos,] <- 0
      curv[,parpos] <- 0
      curv[parpos,parpos] <- -grad(locpredict,x=focalpt)[parpos]^2
      sigma <- - solve(curv) # *(RMSEs[locst])
      rmv <- rmvnorm(nvec[locst],mean=focalpt,sigma=sigma)
      trypoints <- rbind(trypoints,rmv)
    }
  }
  # those for which no CI bound is available:
  for (st in object$CIobject$missingBounds) {
    if (nvec[st]>0L) {
      dotpos <- regexpr(".",st,fixed=TRUE)[1L]
      parm <- substring(st,dotpos+1L)
      side <- substring(st,0,dotpos-1L)
      focalpt <- object$MSL$MSLE
      if (side=="low") {
        bnd <- object$LOWER[parm]
      } else bnd <- object$UPPER[parm]
      if (abs(bnd-focalpt[parm])>(object$UPPER[parm]-object$LOWER[parm])/1e4) { ## if MSLE not at boundary
        focalpt[parm] <- bnd*0.999+0.001*focalpt[parm] ## create refpoint near boundary 
        parpos <- which(fittedNames==parm)
        curv <- hessian(locpredict,x=focalpt)
        curv[parpos,] <- 0
        curv[,parpos] <- 0
        curv[parpos,parpos] <- -grad(locpredict,x=focalpt)[parpos]^2
        sigma <- - solve(curv) ## proceeds as if RMSE = 1
        rmv <- rmvnorm(nvec[st],mean=focalpt,sigma=sigma)
        trypoints <- rbind(trypoints,rmv)
      }
    }
  }
  colnames(trypoints) <- fittedNames
  for (vv in fittedNames) trypoints <- trypoints[trypoints[,vv]>object$LOWER[vv],]
  for (vv in fittedNames) trypoints <- trypoints[trypoints[,vv]<object$UPPER[vv],]
  trypoints <- (cbind(trypoints,object$colTypes$fixedPars)) ## add fixedPars for simulation
  trypoints <- trypoints[,object$colTypes$allPars,drop=FALSE] ## column reordering
  return(trypoints)
}
#debug(rparam_from_SLik_j)
#rparam_from_SLik_j(mjd,7)

profile.SLik_j <- function(fitted,...) profile.SLik(fitted=fitted,...) ## no need for two distinct methods here


.boo_SLik_joint <- function(simul, debug_level=0L, boo, stat.obs, nbCluster_SVs) {
  boo[] <- simul
  if (debug_level<2L) { # return valid, NULL or try-error
    densv <- suppressWarnings( ## suppress warnings for clustering failure
      try(infer_SLik_joint(boo,stat.obs=stat.obs, nbCluster=nbCluster_SVs,
                           verbose=list(most=FALSE,final=FALSE)),silent=TRUE))
    if (inherits(densv,"try-error") && debug_level==0L) {
      return(NULL) ## used in whichvalid <- which( ! sapply(resu,is.null)) below
    } else return(densv)
  } else { # return error -> useful in serial mode only
    densv <- suppressWarnings( ## suppress warnings for clustering failure
      infer_SLik_joint(boo,stat.obs=stat.obs, nbCluster=nbCluster_SVs,
                       verbose=list(most=FALSE,final=FALSE)))
    return(densv)
  }
}


# Called by .RMSEwrapper.SLik_j();
# this does not call the process-simulating function. Instead, it performs a sort of parametric bootstrap from the Gaussian mixture model,
# i.e., it jointly simulates parameters and projected statistics from the jointdens...
# An artefact is that it can draw meaningless values from the jointdens, e.g. negative values for variance parameters
.boot.SLik_j <- function(object,nsim=2L, force=FALSE, verbose = TRUE, seed=NULL, 
                         parent_cores_info=NULL, # defined cluster by parent *.boot.SLik_j* calls
                         packages=attr(object$logLs,"packages"), env=attr(object$logLs,"env"), # not required if parent_cores_info is provided.
                         cluster_args=list(),
                         cl_seed=.update_seed(object)
                         ) {
  if (!force && nsim<2L) stop("'nsim' must be > 1")
    if (inherits(object$jointdens,"densityMclust")) {
    bootrepls <- replicate(nsim, do.call("sim", c(object$jointdens[c("modelName", "parameters")],
                                                  list(n=nrow(object$logLs))))[,-1L],simplify=FALSE)
    nbCluster_SVs <- list(jointdens=object$jointdens$G,
                      pardens=object$pardens$G)
  } else {
    bootrepls <- .simulate.MixmodResults(object$jointdens, nsim=nsim, size=nrow(object$logLs),
                                         drop=FALSE) 
    nbCluster_SVs <- list(jointdens=object$jointdens@nbCluster,
                      pardens=object$pardens@nbCluster) # single values so that this exposes to clustering failure.
  } ##  marginal=current distr of param
  
  boo <- object$logLs[,with(object$colTypes,c(fittedPars,statNames))]
  attr(boo, "allPars") <- object$colTypes$fittedPars # late 2021/01/25 correction as otherwise infer_SLik_joint() tries to use the original allPars that include fixed ones
  stat.obs <- attr(object$logLs,"stat.obs")
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
  resu <- pblapply(bootrepls, .boo_SLik_joint, cl = cl, boo=boo, stat.obs=stat.obs, nbCluster_SVs=nbCluster_SVs)
  whichvalid <- which( ! sapply(resu,is.null))
  resu <- resu[whichvalid]
  if (length(resu)==0L) {
    message("All bootstrap replicates failed; this suggests a problem that cannot be solved by computing more replicates.\n Trying to diagnose the problem...")
    if (cores_info$nb_cores>1L) {
      message("\nIn parallel mode, first replicate gives:")
      resu <- pblapply(bootrepls[1],.boo_SLik_joint,cl = cores_info$cl, debug_level=1, boo=boo, stat.obs=stat.obs, nbCluster_SVs=nbCluster_SVs)
      print(resu[[1]][1])
      if (identical(cluster_args$debug_info,TRUE)) { # ready-to-use debug code:
        sessioninfo <- utils::sessionInfo()
        infoptions <- Infusion.options()
        globoptions <- options()
        tmpname <-  .generateFileName("debug_info",ext=".rda")
        print(paste0("Saving info in file ", tmpname,":"))
        save(bootrepls, sessioninfo, infoptions, globoptions, object, nbCluster_SVs, verbose, prevmsglength,
             file=tmpname)
        # loading the .rda allows to run .boo_SLik_joint(bootrepls[[1]]) (and controlling the options, if needed)
      } else message(paste("\nSet cluster_args$debug_info=TRUE to save some debugging info in a file\n",
                           "and see the source of Infusion:::.boot.SLik_j() for further details about it.")) 
      message("\nTesting whether this generates an error in serial mode:") # might not fail if pb only in parallel mode
    } else  message("\nIn serial mode, first replicate gives:") # should fail again
    abyss <- lapply(bootrepls[1],.boo_SLik_joint, debug_level=2, boo=boo, stat.obs=stat.obs, nbCluster_SVs=nbCluster_SVs)
    return(NULL) # Only if parallel failed and serial did not. 
  }
  while (length(resu)< nsim) {
    message(paste("Mixture modelling with given nbCluster failed for",nsim-length(resu),"replicate(s); drawing sample(s) again..."))
    moreresu <- .boot.SLik_j(object, nsim=nsim-length(resu), force=TRUE, verbose=verbose, parent_cores_info=cores_info) # recursive call uses cores_info but not cluster_args.
    resu <- c(resu,moreresu)
  }
  if (is.null(parent_cores_info)) .close_cores(cores_info)
  invisible(resu)
}

# Called by MSL()
.RMSEwrapper.SLik_j <- function(object, CIpoints=object$CIobject$bounds, useCI=TRUE, nsim=10L,verbose=interactive(),
                                cluster_args=list()) {
  bootrepls <- .boot.SLik_j(object,nsim=nsim,verbose=verbose,cluster_args=cluster_args) # returns gaussian mixture models for each resample; no MSL as not needed
  if( useCI && ! is.null(CIpoints) ) {
    locdata <- data.frame(rbind(MSLE=object$MSL$MSLE,CIpoints))
    covmat <- cov(do.call(rbind,lapply(bootrepls,predict,newdata=locdata, which="lik"))) ## covmat of predict=of logL
    MSEs <- c(MSL=covmat[1,1],diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  } else {
    locdata <- object$MSL$MSLE
    MSEs <- structure( var( unlist(lapply(bootrepls,predict,newdata=locdata, which="lik")))  ,names="MSL")
  }
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
