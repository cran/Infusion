.RMSEwrapper <- function(object, CIpoints=object$CIobject$bounds, useCI=TRUE, verbose=interactive()) {
  if( useCI && ! is.null(CIpoints) ) {
    locdata <- data.frame(rbind(MSLE=object$MSL$MSLE,CIpoints))
    etas <- predict(object, newdata=locdata,variances=list(predVar=TRUE,cov=TRUE))
    covmat <- attr(etas,"predVar")
    #object$CIobject$cov <- covmat
    MSEs <- c(MSL=covmat[1,1],diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  } else {
    etas <- predict(object,newdata=object$MSL$MSLE,variances=list(predVar=TRUE))
    MSEs <- structure(attr(etas,"predVar")[1L],names="MSL") ## of predict= of logL
  }
  if ( any(MSEs<0) ) {
    message("Inaccurate MSE computation, presumably from nearly-singular covariance matrices.")
    # in particular in infer_surface.tailp 
    # infinite lambda => MSE < 0 => nvec[1]<0 => size<0 in rvolTriangulation crashes 
    MSEs[MSEs<0] <- NA ## quick patch
  }  
  if (inherits(object,"SLikp")) {
    dmudeta <- object$fit$family$mu.eta(as.numeric(etas)) ## also removes attributes, quite useful
    MSEs <- MSEs * (dmudeta^2) 
  }
  RMSEs <- sqrt(MSEs)
  if (length(MSEs)>1L) {
    if (inherits(object,"SLik")) {
      headline <-  "* RMSEs of summary log-L maximum and of its ratio at CI bounds: *\n"
    } else if (inherits(object,"SLikp")) {
      headline <-  "* RMSEs of MaxSumm-tail p and of MSQ for CIs: *\n"
    }
  } else {
    if (inherits(object,"SLik")) {
      headline <-  "* RMSEs of summary log-L maximum: *\n"
    } else if (inherits(object,"SLikp")) {
      headline <-  "* RMSE of MaxSumm-tail p: *\n"
    }      
  }
  if(verbose) {
    cat(headline)
    print(RMSEs)
  }
  return(RMSEs)
}

# creates <slik>$par_RMSEs 
.par_RMSEwrapper <- function(object,verbose=interactive()) {
  RMSEs <- object$RMSEs$RMSEs
  CIs <- object$CIobject$CIs
  if (( ! is.null(CIs)) && (lenRMSEs <- length(RMSEs))>1L)  { # lenRMSEs test to detect if more than MSL in object
    pars <- names(CIs)
    par_RMSEs <- rep(NA,lenRMSEs-1L) ## (over)size as the MSEs have no NAs
    par_ests <- rep(NA,lenRMSEs-1L) ## (over)size as the MSEs have no NAs
    names(par_RMSEs) <- names(par_ests) <- names(RMSEs)[-1L]
    map <- c(low=1L,up=2L)
    for (nam in names(par_RMSEs)) {
      stterms <- regmatches(nam, regexpr("\\.", nam), invert = TRUE)[[1]] ## regmatches(<...>) splits at the first '.'
      parm <- stterms[2L]
      interval <- .get_ci_info(CIs, parm)$interval
      th <- interval[map[stterms[1L]]]
      par_ests[nam] <- th
      if (! is.na(RMSEs[nam])) {
        objfn <- function(value) {
          names(value) <- parm
          profile(object, value=value)
        } 
        dlogLdth <- grad(objfn,x=th)
        par_RMSEs[nam] <- RMSEs[nam]/abs(dlogLdth)
      } else par_RMSEs[nam] <- NA
    }
    resu <- rbind(par=par_ests,par_RMSE=par_RMSEs,LR_RMSE=(RMSEs)[-1L]) 
    if (verbose && length(par_RMSEs)>0L) {
      par_headline <- "*** Interval estimates and RMSEs ***\n"
      cat(par_headline)
      print(resu)
      return(invisible(resu))
    } else return(resu)
  } else return(NULL)
}

.check_identifiables <- function(hess, parnames, ad_hoc_mess="par"){
  npar <- ncol(hess)
  ev  <-  eigen(hess)
  zero.tol <- sqrt(.Machine$double.eps/7e-7) # the one sought for elements of the hessian matrix by hessian() default method. 
  which_lowev <-  which( ev$values < zero.tol*npar )
  if (length(which_lowev)) {
    which_lowest <- which.min(ev$values)
    absv <- abs(ev$vectors[,which_lowest])
    ord <- order(absv, decreasing = TRUE)
    if (npar >1L && absv[ord[2]]>0.001) { 
      bad <- parnames[ord[1:2]]
      if (ad_hoc_mess=="stat") {
        list(non=bad, message=paste0("statistic '",bad[1],"' appears linearly dependent with '",bad[2],"' and possibly some others."))
      } else list(non=bad, message=paste0("parameter '",bad[1],"' may be practically unidentifiable, possibly in combination with '",bad[2],"' and some others."))
    } else {
      bad <- parnames[ord[1]]
      list(non=bad, message=paste0("parameter '",bad[1],"' may be practically unidentifiable."))
    }
  } else return(list(non=NULL,message=NULL))
  # pb <- list()
  # for (jt in which_lowev) {
  #   absv <- abs(ev$vectors[, jt])
  #   pbj <- absv > 0.1/npar
  #   pb[[as.character(jt)]] <- parnames[ which(pbj)] # identifiable parameters should have 0 weight in these eigenvectors 
  # }
  # list(non=pb)
  
}

check_raw_stats <- 
  function(x, statNames, remove=FALSE, verbose=interactive()) {
    if (requireNamespace("caret",quietly=TRUE)) {
      findLinearCombos <- get("findLinearCombos", asNamespace("caret")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      if (inherits(x,"reftable") && missing(statNames)) {
        parnames <- names(attr(x,"LOWER"))
        statNames <- setdiff(colnames(x),parnames)
      } else parnames <- setdiff(colnames(x),statNames) # statNames required here
      rawstats <- x[,statNames,drop=FALSE]
      lindeps <- caret::findLinearCombos(rawstats) # but in fact gllim also produces cases where the variance of rawstats is 0
      lindeps$remove <- statNames[lindeps$remove]
      if (length(lindeps$remove)) {
        if (remove) {
          if (verbose) print(paste0("Variables ",paste(lindeps$remove,collapse=", "),
                                    ",\n which induce linear dependencies between variables, are removed."), quote=FALSE)
          ok <- setdiff(colnames(x),lindeps$remove)
          x[,ok, drop=FALSE]
        } else {
          print(paste0("Variables ",paste(lindeps$remove,collapse=", "),
                       " induce linear dependencies between variables, and would better be removed."), quote=FALSE)
          for (it in seq_along(lindeps$linearCombos)) lindeps$linearCombos[[it]] <- statNames[lindeps$linearCombos[[it]]]
          lindeps
        }
      } else {
        if (verbose >= 1L) print("No linear dependencies detected", quote=FALSE) # allows 0.5 ...
        if (remove) {
          x
        } else {
          lindeps
        }
      }
    } else {
      if ( ! .one_time_warnings$caret_warned ) {
        message(paste("  If the 'caret' package were installed, linear dependencies among raw statistics could be checked."))
        .one_time_warnings$caret_warned <- TRUE
      }
      if (remove) {
        x
      } else {
        "package 'caret' not available for checking linear dependencies"
      }
    }
  }

.inits_from_postdens <- function(object, max_base=.Infusion.data$options$max_base, 
                                 nsim= min(2L*nrow(object$logLs),max_base),
                                 newobs=NULL, 
                                 givenpars=NULL, # vector 
                                 profiledNames,
                                 seed=NULL) {
  if (is.null(newobs)) {
    givens <- c(givenpars, get_from(object,"stat.obs"))
  } else givens <- c(givenpars,drop(unlist(newobs))) # where unlist(<1-row data frame>)
  instr_postdens <- .get_instr_postdens(object, given=givens)
  # For MAF we need the conditioning value ie sobs passed asa 1-row matrix
  if (inherits(instr_postdens,"MAF") && 
      attr(instr_postdens,"which")=="I_postdens") givenpars <- t(givens)
  inits <- .sample_in_absol_constrs(object, ceil_size = nsim, profiledNames, 
                                    density=instr_postdens, givenpars=givenpars, 
                                    norm_or_t=.wrap_rmvnorm, seed=seed)
}

# For profiles:
# .get_init_from_hull <- function(object, given) {
#   X <- object$logLs[,object$colTypes$fittedPars]
#   parname <- names(given)
#   MSLE <- object$MSL$MSLE
#   ontheright <- given>MSLE[parname]
#   ### Find best point beyond the 'given' value:
#   if (ontheright) {
#     subX <- X[X[,parname]>given,]
#   } else subX <- X[X[,parname]<given,]
#   if ( ! nrow(subX)) return(NULL)
#   rownames(X) <- seq(nrow(X))
#   preds_safe <- predict(object, subX, which="safe")
#   outer <- which.max(preds_safe)
#   logL_outer <- preds_safe[outer]
#   outer <- subX[outer,] # parameter point
#   ###
#   ### Find best point within the MSLE-'given' range, close to 'given':
#   if (ontheright) {
#     subX <- X[X[,parname]>MSLE[parname] & X[,parname]<outer[,parname],]
#   } else subX <- X[X[,parname]<MSLE[parname] & X[,parname]>outer[,parname],]
#   if ( ! nrow(subX)) return(NULL)
#   preds_safe <- predict(object, subX, which="safe")
#   top <- which(preds_safe>logL_outer)
#   if (length(top)>10L) subX <- subX[top,] # else no subselection of subX
#   subX <- rbind(subX, outer)
#   xy <- cbind(x=subX[,parname],
#               logL=predict(object, subX, which="safe"))
#   if (TRUE) {
#     convhull_info <- convhulln(xy)
#     convhullids <- unique(as.vector(convhull_info))
#   } else {
#     convhull_info <- resetCHull(xy, quickndirty = TRUE)$vertices
#     convhullids <- rownames(convhull_info)
#   }
#   inn <- subX[names(which(predict(object, subX[convhullids,], which="safe") >
#                             logL_outer+1e-6 # => exclude 'outer' point
#   )),, drop=FALSE]
#   if ( ! nrow(inn)) return(NULL)
#   if (ontheright) {
#     inner <- inn[which.max(inn[,parname]),]
#   } else inner <- inn[which.min(inn[,parname]),]
#   # Must respect 'given':
#   init <- (abs(outer[,parname]-given)*inner+abs(given - inner[,parname])*outer)/ abs(outer[,parname]- inner[,parname])
#   # The level of extrapolation may affect both this 'init' (through the above predict(., which="safe") calls),
#   # and the predicted logL of this 'init'
#   init
# }

# For profiles:
.get_init_from_hull <- function(object, given, 
                                newobs=t(get_from(object,"proj_data")), 
                                MSLE,# =object$MSL$MSLE, # alternative newobs and alternative MSLE used in.ecdf_2lr
                                which="safe") { 
  X <- object$logLs[,object$colTypes$fittedPars]
  parname <- names(given)
  ontheright <- given>MSLE[parname]
  ### Find best point beyond the 'given' value:
  if (ontheright) {
    subX <- X[X[,parname]>given,]
  } else subX <- X[X[,parname]<given,]
  if ( ! nrow(subX)) return(NULL)
  rownames(X) <- seq(nrow(X))
  preds_safe <- summLik(object, parm=subX, data=newobs, which=which)
  outer <- which.max(preds_safe)
  logL_outer <- preds_safe[outer] # max logL for points outside the interval
  outer <- subX[outer,] # parameter point
  ###
  ### Find best point within the MSLE-'given' range, close to 'given':
  if (ontheright) {
    subX <- X[X[,parname]>MSLE[parname] & X[,parname]<given,]
  } else subX <- X[X[,parname]<MSLE[parname] & X[,parname]>given,]
  if ( ! nrow(subX)) return(NULL)
  if (nrow(subX)>1L) {
    preds_safe <- summLik(object, parm=subX, data=newobs, which=which)
    top <- which(preds_safe>logL_outer)
    if (length(top)>10L) {
      subX <- subX[top,] 
      preds_safe <- preds_safe[top]
    } # else no subselection of subX
    subX <- rbind(subX, outer)
    xy <- cbind(x=subX[,parname], 
                logL=c(preds_safe, logL_outer))
    if (TRUE) {
      convhull_info <- convhulln(xy)
      convhullids <- unique(as.vector(convhull_info))
    } else {
      convhull_info <- resetCHull(xy, quickndirty = TRUE)$vertices
      convhullids <- rownames(convhull_info)
    }
    convhullids <- setdiff(convhullids, nrow(xy)) # => exclude 'outer' point
    hullX <- subX[convhullids,, drop=FALSE]
    if (ontheright) {
      inner <- hullX[which.max(hullX[,parname]),]
    } else inner <- hullX[which.min(hullX[,parname]),]
  } else inner <- subX
  # Must respect 'given':
  init <- (abs(outer[,parname]-given)*inner+abs(given - inner[,parname])*outer)/ abs(outer[,parname]- inner[,parname])
  # The level of extrapolation may affect both this 'init' (through the above predict(., which="safe") calls),
  # and the predicted logL of this 'init'
  init
}


        
# .safe_init() -> predict(,"safe") assumes that object$thr_info has been computed, which is done by MSL()
.safe_init <- function(object, 
                       given=NULL, # in MSL() in particular: no 'given' nor 'newobs".
                       plower, pupper, 
                       newobs=NULL, 
                       more_inits=NULL, # *full-length*; to add to locally generated ones. It should be rbind()able -> dim attr
                       base_inits=NULL, # possible alternative to calling .inits_from_postdens().
                       profiledNames=names(plower),
                       max_base=.Infusion.data$options$max_base, # Inf by default, 
                       nsim=min(2L*nrow(object$logLs),max_base),
                       constr_crits=object$constr_crits,
                       MSLE, # for .get_init_from_hull
                       seed=NULL) { # 
  if (is.null(base_inits)) { # ideally always TRUE for clarity.
    if ( ! is.null(given) && length(given)==1L) { # 1D profile CIs => allow efficient procedure
      # .inits_from_postdens() appears INefficient to generate good points for non-NULL 'given',
      # so we first try this (major improvement)
      if ( ! is.null(newobs)) {
        if (missing(MSLE)) stop("The 'MSLE' argument is required when calling .safe_init(.,given=<one parm>, newobs).")
        if ( ! is.null(MSLE)) base_inits <- .get_init_from_hull(object=object, given=given,
                                                                newobs=newobs, 
                                                                MSLE=MSLE)
        # if MSLE is explicit NULL, no base_inits is yet computed
      } else base_inits <- .get_init_from_hull(object=object, given=given,
                                        newobs=t(get_from(object,"proj_data")), 
                                        MSLE=object$MSL$MSLE) # may be NULL
    } 
    if (is.null(base_inits) || # 1st cond: NOT 1D profile nor SLRT 1D; rather, plot2Dprof() or SLRT !1D
        # SLRT 1D -> init_h0 with given and newobs -> a base init has just been provided 
        # so 1st cond is FALSE. But we can force call of .inits_from_postdens 
        # by the second condition (SLRT 1D):
        ! missing(nsim) 
      ) { 
      base_inits <- .inits_from_postdens(
        object, max_base=max_base, 
        nsim= nsim, # default is large...
        newobs=newobs, seed=seed, givenpars=given, profiledNames=profiledNames)
      ## Does not seem useful:
      # base_inits <- rbind(object$logLs[,object$colTypes$fittedPars],
      #                     base_inits)
    }
  } else { # tried for MSL call in particular: cf also .ecdf_2lr()
    nr <- nrow(base_inits)
    if (nr>max_base) base_inits <- base_inits[(nr-max_base):nr ,, drop=FALSE]
  }
  inits <- rbind(more_inits, base_inits, object$MSL$MSLE)
  for (st in names(given)) inits[,st] <- given[st]
  if ( ! is.null(constr_crits)) { # selection of points that satisfy parameter constraints
    dx <- pupper-plower
    margin <- dx/10000
    safelow <- plower+margin
    safeup <- pupper-margin
    for(st in names(plower)) inits[,st] <- pmax(safelow[st],pmin(safeup[st],inits[,st]))
    constrs <- apply(inits,1L, function(v) all(eval(constr_crits, envir = as.list(v))<0))
    inits <- inits[constrs,, drop=FALSE]
  }
  if (nrow(inits)) {
    if (is.null(newobs)) { 
      inits_pred <- predict(object,inits, which="safe", constr_tuning=FALSE)
      # if ( ! is.null(given)) {
      #   init <- .init_params_from_pardens(object = object, given, profiledNames, plower, pupper)
      #   if ( (  ! is.null(constr_crits) && any(eval(constr_crits, envir = as.list(c(init,given)))>0)) ||
      #        max(inits_pred) > 
      #        predict(object,c(init[profiledNames],given), which="safe", constr_tuning=FALSE)
      #   ) init <- unlist(inits[which.max(inits_pred),])
      # } else init <- unlist(inits[which.max(inits_pred),, drop=FALSE]) 
      # in MSL() in particular: no 'given' nor 'newobs". drop=FALSE for 1-parameter case (otherwise parName is lost) 
    } else {
      inits_pred <- summLik(object,inits, data=newobs, which="safe")
      # if ( ! is.null(given)) {
      #   init <- .init_params_from_pardens(object = object, given, profiledNames, plower, pupper, newobs=newobs)
      #   if ( (  ! is.null(constr_crits) && any(eval(constr_crits, envir = as.list(c(init, given)))>0)) ||
      #        max(inits_pred) > summLik(object,c(init[profiledNames],given), data=newobs, which="safe")
      #   ) init <- unlist(inits[which.max(inits_pred),])
      # } else init <- unlist(inits[which.max(inits_pred),])
    }
    if (inherits(inits, "data.frame")) { #  was always the case prior to using .inits_from_postdens():
      # still the case for xLLiM clustering
      init <- unlist(inits[which.max(inits_pred),,drop=FALSE]) # drop=FALSE for 1-parameter case (otherwise parName is lost) 
    } else init <- inits[which.max(inits_pred),] 
  } else {
    candidates <- object$logLs[,names(object$LOWER)]
    for (st in names(given)) candidates[,st] <- given[st]
    constrs <- apply(candidates,1L, function(v) all(eval(constr_crits, envir = as.list(v))<0))
    candidates <- candidates[constrs,, drop=FALSE]
    nr <- nrow(candidates)
    if (nr>1L) {
      init <- candidates[which.max(predict(object, newdata=candidates, which = "safe")), ]
    } else if (nr) {
      init <- candidates[, ]
    } else return("No init point found to satisfy the constraints.") # if nloptr's init is 'outside' it returns a the solution which may be 'outside'
    init <- unlist(init)
  }
  init[profiledNames] # named vector ordered as plower
}

.get_gradfn <- function(object, init) { # grad of  - logL
  col_ids <- which(object$colTypes$fittedPars %in% names(init))
  nbCluster <- object$jointdens@nbCluster
  logproportions <- object$clu_params$logproportions
  solve_t_chol_sigma_lists <- object$clu_params$solve_t_chol_sigma_lists
  stat.obs <- get_from(object,"stat.obs")
  gradfn <- function(v) { 
    grad_log_pardens <- .grad.dMixmod(object=object$pardens,newdata=v,
                                          logproportions=logproportions, clu_means=object$clu_params$pardens_means,
                                          nbCluster=nbCluster, col_ids=col_ids, solve_t_chol_sigma_list=solve_t_chol_sigma_lists$pardens)
    grad_log_jointdens <- .grad.dMixmod(object=object$jointdens,newdata=c(v, stat.obs),
                    logproportions=logproportions, clu_means=object$clu_params$jointdens_means,
                    nbCluster=nbCluster, col_ids=col_ids, solve_t_chol_sigma_list=solve_t_chol_sigma_lists$jointdens) 
    # return value could be further penalized by - [parvaldens<object$thr]* grad_log_pardens 
    # so that it would be the grad of the "safe" logL. This would require parvaldens being available locally.
    grad_log_pardens - grad_log_jointdens #  of  - logL
  }  
  gradfn # returns the FUNCTION
}

# Will use gradient only if "safe' is not TRUE => gradient not used in practice
.ad_hoc_opt <- function(init, objfn, lower, upper, safe=TRUE, LowUp=list(), gradient, 
                        object, neg_ineq_constrfn, template=NULL, ...) {
  time1 <- Sys.time()
  if (safe) {
    msl <- .safe_optim(init, objfn=objfn, 
                     lower=lower,upper=upper, LowUp=LowUp, verbose=FALSE, object=object, 
                     neg_ineq_constrfn=neg_ineq_constrfn, template=template,...)
    # msl$par <- msl$solution
    if ("template" %in% names(formals(objfn))) {
      obj_init <- objfn(init, template=template)
    } else obj_init <- objfn(init)
    if (obj_init<msl$objective) { # bobyqa at a boundary...
      msl$solution <- init
      msl$objective <- obj_init
    }
  } else { 
    msl <- nlminb(init, objfn, lower=lower, upper=upper, gradient=gradient, template=template, ... )
    msl$solution <- msl$par
  }
  msl$optim_time <- round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1) ## spaMM:::.timerraw(time1)
  msl$value <- - msl$objective # $value: logL value; $objective: minimization objective. Trying to change this was confusing.
  msl
}

.calc_pardens_thr_info <- function(object, fittedPars=object$colTypes$fittedPars, level, 
                            pred_data=object$logLs[,fittedPars, drop=FALSE],
                            thr_info_ctrl=.Infusion.data$options$thr_info_ctrl) {
  logls <- predict(object,newdata = pred_data, which="lik", constr_tuning=FALSE)
  logls_order <- order(logls,decreasing = TRUE)
  whichmaxlogl <- logls_order[as.integer(1L+eval(thr_info_ctrl))]
  maxlogl <- logls[whichmaxlogl]
  if (is.infinite(maxlogl)) stop("is.infinite(maxlogl)") # gllim problem in bootstrap replicate...
  ## If I set thr_dpar <- predict(object,newdata = pred_data[whichmaxlogl,], which="parvaldens", ...)
  ## I see an artefact [see odp plot in devel/explanations]: profile values for several parameters tend to be the 'maxlogl' => identical p values for several SLRTs
  ## This occurs when pardens for the tested value is practically always below the dpar_thr. 
  ## * In that case any value with lik>maxlogl will be penalized, so the maximization is over
  ## parms with lik<maxlogl, and if lik>maxlogl easily at low densities, 
  ## we just get max(lik | lik<maxlogl) = maxlogl.
  ##
  ## * This artefact suggests we need more simul to increase parvaldens for tested value.
  ##
  ## Practical solution below it to allow lower thr_dpar by exploring the putative CI range
  # which_in_putative_CI <- which(logls> maxlogl-qchisq(level, df=1)/2)
  # thr_dpar <- min(predict(object,newdata = pred_data[which_in_putative_CI,], which="parvaldens", constr_tuning=FALSE))
  thr <- qchisq(level, df=1)/2 
  which_in_top_slice <- which(logls> maxlogl-thr)
  ## trivial protection, but not optimized in any sense:
  while (length(which_in_top_slice)<2L) {
    thr <- 2*thr
    which_in_top_slice <- which(logls> maxlogl-thr)
  }
  dpar_range <- range(predict(object,newdata = pred_data[which_in_top_slice,], which="parvaldens", constr_tuning=FALSE))
  #
  # If I instead penalized lik even when it is<maxlogl (how?), and assuming the MSLE has high pardens,
  # this will lower the profile values and the p-values.
  # Set the threshold for "safe" logL computation as function of parameter density:
  # list2env(
    list(thr_dpar=dpar_range[1], # sum(dpar_range*c(0.75,0.25)),        # threshold value of parameter density
         reft_maxlogl=maxlogl,
         # next ones maybe not used but useful in debug session
         par_thr=pred_data[logls_order[1],], # allowed 'max'
         par_lik=pred_data[whichmaxlogl,] # observed max 'lik'
         ) 
  #  , 
  #        parent=emptyenv()) # envir so that it can be updated similarly to the the MSL envir
  ## Result is NOT modified when a new maximum is found, as it depends only on the logls of the reftable points
  ## Attempts to hack it differently were not successful and rather suggested fixing other issues (such as using
  ## nsteps in SLRT, or more deeply, by better parameter space exploration:
  ## B_13from17, replicate 49, fit by v2.1.106 vs 2.1.112)
}

.get_thr_info <- function(object) {
  if (is.null(thr_info <- object$thr_info)) { # handles older objects with only object$thr_dpar
    stop("Run MSL() to update the object for use with recent version of Infusion")
  } 
  thr_info
}

# both SLik and SLikp, with different methods used in -> allCIs -> confint
MSL <- function (object,
                 CIs=prod(dim(object$logLs))<12000L, level=0.95, 
                 verbose=interactive(),
                 RMSE_n=Infusion.getOption("RMSE_nsim"), 
                 eval_RMSEs=(RMSE_n>1L) * prod(dim(object$logLs))<12000L, # RMSEs on likelihood (ratio) values
                 cluster_args=list(), nb_cores=NULL, init=NULL, prior_logL=NULL) { ##
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- (c("most","from_refine"))[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- (interactive())
  if (is.null(verbose$from_refine)) verbose$from_refine <- FALSE
  
  fittedPars <- object$colTypes$fittedPars
  if (inherits(object,"SLik")) {
    vertices <- object$fit$data[,fittedPars,drop=FALSE]
    lowup <- sapply(vertices,range)
    lower <- lowup[1L,]
    upper <- lowup[2L,]
  } else {
    if (verbose$most && ! verbose$from_refine) cat(paste(nrow(object$logLs),"simulated samples;\n"))
    lower <- object$lower
    upper <- object$upper
  }
  # (a constrOptim would be even better)
  parscale <- (upper-lower) 
  ## pred_data was previously used for .calc_thr_dpar() [replaced by .calc_pardens_thr_info] and for init
  ## It turns out that this may not contain the best point anyway.
  # nr <- nrow(object$logLs)
  # recent <- seq(max(1L,nr-1000L),nr)
  # pred_data <- object$logLs[recent,fittedPars, drop=FALSE] 
  object$thr_info <- .calc_pardens_thr_info(object=object, fittedPars=fittedPars, level=level)
  if (inherits(object,"SLik_j")) { # reset all envs that keep results of boot with BGP potentially = MSLE 
    object$bootLRTenv <- list2env(list(bootreps_list=list()),
                                  parent=emptyenv())
    object$bootCIenv <- list2env(list(bootreps_list=list()),
                                   parent=emptyenv())
    object$CIobject <- list2env(list(CIs=NULL),
                                parent=emptyenv())
  }
  if (is.null(init)) {
    if (inherits(object,"SLik_j")) {
      # Here .safe_init() will use a large draw from instr_postdens. 
      # This still misses local maxima; using the reftable does not help...
      init <- .safe_init(object = object, # base_inits=pred_data, 
                         profiledNames = fittedPars,
                         # could add non-NULL seed to control repeatability
                         plower=lower, pupper=upper)
    } else {
      init <- unlist(object$obspred[which.max(object$obspred[,attr(object$obspred,"fittedName")]),fittedPars]) 
      init <- init*0.999+ colMeans(vertices)*0.001
    }
  }
  prev_init <- object$MSL$init_from_prof ## provided by plot1Dprof() -> profil(). 
  ## prev_init is not NULL if MSL is here called by .MSL_update() or by [refine() after a plot], but then object$MSL$MSLE should not be used ;
  ## in other calls by refine(), the object is being reconstructed from scratch, and then object$MSL$MSLE is NULL too ;
  ## in neither case the next line is useful
  # if (is.null(prev_init)) prev_init <- object$MSL$MSLE 
  ##   which="safe" is effective only for SLik_j objects. Otherwise, it is ignored.
  constr_tuning <- FALSE
  objectivefn <- function(v) { 
    - predict(object,newdata=v, which="safe", control=list(fix_predVar=TRUE), constr_tuning=constr_tuning)
  } 
  if ( ! is.null(constr_crits <- object$constr_crits)) {
    neg_ineq_constrfn <- function(v) {
      names(v) <- names(init)
      as.numeric(eval(constr_crits, envir = as.list(v)))
    }
  } else neg_ineq_constrfn <- NULL
  
  
  # fix_predVar=TRUE inhibits messages from spaMM:::.get_invColdoldList. See comments there. (effective for some example in the vignette, of old workflow -> MSL.Slik -> ... predict.HLfit())
  if ( (! is.null(prev_init)) &&
       objectivefn(prev_init) <  objectivefn(init)) {init <- prev_init}
  
  # if (essai) {
  #   msl <- .try_optim_clu_it(object, proj_data=get_from(object,"proj_data")) # [,object$jointdens@statNames])
  # } else {
    #### optimization call of up-to-date workflow. Handles constraints defined by constr_crits (.ad_hoc_opt -> .safe_optim -> iterative algo )
    msl <- .ad_hoc_opt(init, objfn=objectivefn, lower=lower, upper=upper,
                       gradient=.get_gradfn(object, init), object=object, neg_ineq_constrfn=neg_ineq_constrfn)
  # }
  
  if (inherits(object,"SLik") && length(fittedPars)>1L) { # This should NOT be run for SLik_j objects
    locchull <- resetCHull(vertices, formats=c("constraints"))
    if( ! (isPointInCHull(msl$solution, constraints=locchull[c("a", "b")]))) { ## if simple optim result not in convex hull
      # this block can be tested by test-3par, with mixmodGaussianModel="Gaussian_pk_Lk_Ck"
      ui <- -locchull$a
      ci <- -locchull$b
      candttes_in_hull <- 0.99*vertices + 0.01*colMeans(vertices)
      best_candidate <- unlist(candttes_in_hull[which.max(predict(object,newdata=candttes_in_hull, which="safe")),])
      constr_tuning <- FALSE
      objectivefn <- function(v) { 
        - as.numeric(predict(object,newdata=v, which="safe", control=list(fix_predVar=TRUE), constr_tuning=constr_tuning))
      }
      objectivefn.grad <- function(x) {grad(func=objectivefn, x=x, constr_tuning=constr_tuning)} ## no need to specify fixedlist here
      if (.Infusion.data$options$constrOptim) {
        msl <- constrOptim(theta=best_candidate,f=objectivefn,grad=objectivefn.grad, ui=ui, ci=ci , 
                           mu=1e-08, ## a low mu appear important
                           method = "BFGS",control=list(parscale=parscale))
        msl$value <- - msl$value
      } 
    }
  }
  if (verbose$most) {
    if (inherits(object,"SLik")) {
      cat("* Summary ML: *\n")
      print(c(msl$solution,"logL"=msl$value))
    } else if (inherits(object,"SLik_j")) {
      cat("* Summary ML: *\n")
      print(c(msl$solution,"logL"=msl$value))
    } else if (inherits(object,"SLikp")) {
      cat("* Summary MQ: *\n")
      print(c(msl$solution,"tailp"=msl$value))
    }
  }
  if (length(fittedPars)==1L) names(msl$solution) <- fittedPars
  # hess ideally SPD (as that of - logL)
  hess <- hessian(objectivefn, msl$solution) # numDeriv::hessian can give bad results, even in seemingly inoccuous 1D model
  # comparison with pracma::hessian may be instructive
  colnames(hess) <- rownames(hess) <- fittedPars
  rnge <- diag(x=sqrt(upper-lower), ncol=length(upper))
  hess <- rnge %*% hess %*% rnge
  # chk_identif <- .check_identifiables(hess,parnames=fittedPars) 
  # if (length(chk_identif$non)) message(chk_identif$message)
  #
  sol <- .sanitize_optim_solution(msl$solution, lower, upper)
  object$MSL <- list2env(list(MSLE=sol, maxlogL=msl$value, hessian=hess, 
                              updated_from_prof=(! is.null(prev_init))),
                         parent=emptyenv())
  if (inherits(object$fit,"HLfit")) object$MSL$predVar <- 
    attr(predict(object,newdata=msl$solution,variances=list(predVar=TRUE)),"predVar")
  # CIs
  prevmsglength <- 0L
  if(CIs) {
    locverbose <- (verbose$most && ! inherits(object$fit,"HLfit")## for HLfit object, printing is later
                   && msl$optim_time*length(msl$solution)>3) # guess when it is useful to be verbose from the time to find the maximum 
    if (locverbose) {
      prevmsglength <- .overcat("Computing confidence intervals...\n", 0L) 
    } 
    object$CIobject <- list2env(allCIs(object,verbose=locverbose, level=level),
                                parent=emptyenv())  ## may be NULL
  } 
  # MSEs computation
  if (inherits(object$fit,"HLfit")) {
    object$RMSEs <- list2env(list(RMSEs=.RMSEwrapper(object,verbose=FALSE),
                                  warn=NULL),
                             parent=emptyenv())
  } else {
    if (eval_RMSEs) {
      if  (verbose$most) .overcat("Computing RMSEs... (may be slow)\n",prevmsglength)
      if (is.null(cluster_args$spec)) cluster_args$spec <- nb_cores # which means that cluster_args$spec overrides nb_cores
      object$RMSEs <- list2env(
        list(RMSEs=.RMSEwrapper.SLik_j(object,
                                       boot_nsim=RMSE_n,
                                       cluster_args=cluster_args,verbose=FALSE,
                                       level=level),
             warn=NULL),
        parent=emptyenv())
    } 
  }
  if (is.null(object$logLs$cumul_iter)) {
    object$logLs$cumul_iter <- 1L ## adds a column to the data frame; needed soon for .par_RMSEwrapper() -> . -> profile.SLik()
    attr(object$logLs,"n_first_iter") <- nrow(object$logLs) # but this attribute is not kept in later iterations...
  }
  # par_RMSEs computation requires prior $CIobject information and $RMSEs information:
  object$par_RMSEs <- list2env(list(par_RMSEs=.par_RMSEwrapper(object,verbose=FALSE)),
                               parent=emptyenv())
  if  (verbose$most) {
    if ( ! is.null(object$par_RMSEs$par_RMSEs)) {
      .overcat("*** Interval estimates and RMSEs ***\n",prevmsglength)
      print(object$par_RMSEs$par_RMSEs)
    } else {
      bounds <- .extract_intervals(object,verbose=FALSE) 
      if (length(bounds)) {
        .overcat("*** Interval estimates ***\n",prevmsglength)
        print(bounds)
      }
    }
  }
  if ( ! is.null(prior_logL)) object$prior_logL <- prior_logL
  if ( ! is.null(object$projectors) && is.null(object$timeStamps)) {
    # expected to occur only on the first MSL()
    tmstp <- Sys.time() # not the exact times, but sufficient for their use.
    object$timeStamps <- list2env(list(proj_stats=tmstp, projectors=tmstp), parent=emptyenv())
  }
  object$`Infusion.version` <- packageVersion("Infusion")
  invisible(object)
}

# Hack the object to avoid recomputation of GMM clusterint when new Data are analyzed: 
.update_obs <- function(object, # result of infr_SLik_joint() (or MSL() or refine()... )
                        new.obs, # named vector of summary statistics (ideally projected ones, but raw ones are handled)
                        CIs=FALSE, eval_RMSEs = FALSE, 
                        ... # further arguments passed to MSL(). e.g. CIs=FALSE, eval_RMSEs=FALSE may be useful to speed up simulations studies.
                        ) { 
  object$RMSEs <- NULL
  object$CIobject <- NULL
  if ( ! is.null( cn <- colnames(new.obs))) {
    if (nrow(new.obs)>1L) stop("'new.obs' has several rows, but it should instead be a vector.")
    message("Note: 'new.obs' should be a numeric vector, not a matrix or data.frame. Converting...")
    if (is.data.frame(new.obs)) {
      new.obs <- unlist(new.obs)
    } else new.obs <- drop(new.obs)
    # names(stat.obs) <- cn ## minimal patch so that names() can be used, not colnames()
  }
  if ( ! is.null(projectors <- object$projectors)) {
    statNames <- names(new.obs)
    if (length(setdiff(statNames, names(projectors)))) {
      new.obs <- .project_reftable_raw(new.obs, projectors=projectors, is_trainset=FALSE, use_oob=FALSE)
    }
  }
  attr(object$logLs,"stat.obs") <- new.obs
  MSL(object, CIs=CIs, eval_RMSEs=eval_RMSEs, ...)
}
  

  
