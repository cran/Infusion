.RMSEwrapper <- function(object, CIpoints=object$CIobject$bounds, useCI=TRUE, verbose=interactive()) {
  if( useCI && ! is.null(CIpoints) ) {
    locdata <- data.frame(rbind(MSLE=object$MSL$MSLE,CIpoints))
    etas <- predict(object, newdata=locdata,variances=list(linPred=TRUE,dispVar=TRUE,cov=TRUE))
    covmat <- attr(etas,"predVar")
    #object$CIobject$cov <- covmat
    MSEs <- c(MSL=covmat[1,1],diag(covmat[-1,-1,drop=FALSE])+covmat[1,1]-2*covmat[1,-1])
  } else {
    etas <- predict(object,newdata=object$MSL$MSLE,variances=list(linPred=TRUE,dispVar=TRUE))
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
  if (( ! is.null(CIs)) && (lenRMSEs <- length(RMSEs))>1L)  {
    pars <- names(CIs)
    par_RMSEs <- rep(NA,lenRMSEs-1L) ## (over)size as the MSEs have no NAs
    par_ests <- rep(NA,lenRMSEs-1L) ## (over)size as the MSEs have no NAs
    names(par_RMSEs) <- names(par_ests) <- names(RMSEs)[-1L]
    map <- c(low=1L,up=2L)
    for (nam in names(par_RMSEs)) {
      stterms <- regmatches(nam, regexpr("\\.", nam), invert = TRUE)[[1]] ## regmatches(<...>) splits at the first '.'
      parm <- stterms[2L]
      interval <- CIs[[parm]]$interval
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
    if (absv[ord[2]]>0.001) { 
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

check_raw_stats <- local({
  caret_warned <- FALSE
  function(x, statNames, remove=FALSE, verbose=interactive()) {
    if (requireNamespace("caret",quietly=TRUE)) {
      findLinearCombos <- get("findLinearCombos", asNamespace("caret")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      if (inherits(x,"reftable") && missing(statNames)) {
        parnames <- names(attr(x,"LOWER"))
        statNames <- setdiff(colnames(x),parnames)
      } else parnames <- setdiff(colnames(x),statNames) # statNames required here
      rawstats <- x[,statNames,drop=FALSE]
      lindeps <- caret::findLinearCombos(rawstats)
      lindeps$remove <- statNames[lindeps$remove]
      if (length(lindeps$remove)) {
        if (remove) {
          if (verbose) print(paste0("Variables ",paste(lindeps$remove,collapse=", "),
                                    ",\n which induce linear dependencies between variables, are removed."))
          ok <- setdiff(colnames(x),lindeps$remove)
          x[,ok, drop=FALSE]
        } else {
          print(paste0("Variables ",paste(lindeps$remove,collapse=", "),
                       " induce linear dependencies between variables, and would better be removed."))
          for (it in seq_along(lindeps$linearCombos)) lindeps$linearCombos[[it]] <- statNames[lindeps$linearCombos[[it]]]
          lindeps
        }
      } else {
        if (verbose) print("No linear dependencies detected")
        if (remove) {
          x
        } else {
          lindeps
        }
      }
    } else {
      if ( ! caret_warned ) {
        message(paste("  If the 'caret' package were installed, linear dependencies among raw statistics could be checked."))
        caret_warned <<- TRUE
      }
      if (remove) {
        x
      } else {
        "package 'caret' not available for checking linear dependencies"
      }
    }
  }
})
        
        


# both SLik and SLikp, with different methods used in -> allCIs -> confint
MSL <- function (object,CIs=TRUE,level=0.95, verbose=interactive(),
                 eval_RMSEs=TRUE, #inherits(object,"SLik"),    # RMSEs on likelihood (ratio) values
                 cluster_args=list(), init=NULL, prior_logL=NULL,
                 ...) { ##
  fittedPars <- object$colTypes$fittedPars
  if (inherits(object,"SLik")) {
    vertices <- object$fit$data[,fittedPars,drop=FALSE]
    lowup <- apply(vertices,2L,range)
    lower <- lowup[1L,]
    upper <- lowup[2L,]
  } else {
    if (verbose) cat(paste(nrow(object$logLs),"simulated samples;\n"))
    lower <- object$lower
    upper <- object$upper
  }
  # (a constrOptim would be even better)
  parscale <- (upper-lower) 
  if (is.null(init)) {
    if (inherits(object,"SLik_j")) {
      object$bootLRTenv <- list2env(list(bootreps_list=list()))
      pred_data <- object$logLs[,fittedPars, drop=FALSE] # checks likelihood of all points to initiate maximization. Should be efficient
      if (TRUE || is.null(object$thr_dpar)) {
        pardensv <- predict(object,newdata = pred_data, which="parvaldens")
        # Used by predict(., which="safe") to remove spurious high logL in parameters regions where it is poorly estimated (low parameter density)
        object$thr_dpar <- min(max(pardensv)- qchisq(1-(1-level)/2,df=1)/2 , ## (fixme rethink threshold? Actually, the quantile() is lower and more important) 
                               quantile(pardensv,probs=1/sqrt(length(pardensv)))
        )
        #print(object$thr_dpar)
      } 
      predsafe <- predict(object,newdata=pred_data, which="safe") 
      init <- unlist(pred_data[which.max(predsafe), ])
    } else {
      init <- unlist(object$obspred[which.max(object$obspred[,attr(object$obspred,"fittedName")]),fittedPars]) 
      init <- init*0.999+ colMeans(vertices)*0.001
    }
  }
  prev_init <- object$MSL$init_from_prof ## provided by plot1Dprof() -> profil(). 
  ## prev_init is not NULL if MSL called by .MSL_update() or by [refine() after a plot], but then object$MSL$MSLE should not be used ;
  ## in other calls by refine(), the object is being reconstructed from scratch, and then object$MSL$MSLE is NULL too ;
  ## in neither case the next line is useful
  # if (is.null(prev_init)) prev_init <- object$MSL$MSLE 
  ##   which="safe" is effective only for SLik_j objects. Otherwise, it is ignored.
  objectivefn <- function(v) { - predict(object,newdata=v, which="safe", control=list(fix_predVar=TRUE))} 
     # fix_predVar=TRUE inhibits messages from spaMM:::.get_invColdoldList. See comments there. (effective for some example in the vignette, of old workflow -> MSL.Slik -> ... predict.HLfit())
  if ( (! is.null(prev_init)) &&
       objectivefn(prev_init) <  objectivefn(init)) {init <- prev_init}
  time1 <- Sys.time()
  msl <- .safe_opt(init, objfn=objectivefn, 
                   lower=lower,upper=upper, LowUp=list(), verbose=FALSE)
  optim_time <- round(as.numeric(difftime(Sys.time(), time1, units = "secs")), 1) ## spaMM:::.timerraw(time1)
  msl$value <- - msl$objective 
  msl$par <- msl$solution
  if (inherits(object,"SLik") && length(fittedPars)>1L) { # for SLik_j objects, this should not be run
    locchull <- resetCHull(vertices, formats=c("constraints"))
    if( ! (isPointInCHull(msl$par, constraints=locchull[c("a", "b")]))) { ## if simple optim result not in convex hull
      # this block can be tested by test-3par, with mixmodGaussianModel="Gaussian_pk_Lk_Ck"
      ui <- -locchull$a
      ci <- -locchull$b
      candttes_in_hull <- 0.99*vertices + 0.01*colMeans(vertices)
      best_candidate <- unlist(candttes_in_hull[which.max(predict(object,newdata=candttes_in_hull, which="safe")),])
      objectivefn <- function(v) { - as.numeric(predict(object,newdata=v, which="safe", control=list(fix_predVar=TRUE)))}
      objectivefn.grad <- function(x) {grad(func=objectivefn, x=x)} ## no need to specify fixedlist here
      if (TRUE || .Infusion.data$options$constrOptim) {
        msl <- constrOptim(theta=best_candidate,f=objectivefn,grad=objectivefn.grad, ui=ui, ci=ci , 
                           mu=1e-08, ## a low mu appear important
                           method = "BFGS",control=list(parscale=parscale))
        msl$value <- - msl$value
      } else { ## more or less works, but the result is not as accurate as by constrOptim. Presumably related to the  method used, with bound constraints.
        msl <- .safe_constrOptim(theta=best_candidate,f=objectivefn,grad=objectivefn.grad, ui=ui, ci=ci , 
                                 lower=lower,upper=upper,
                           mu=1e-08)
        msl$value <- - msl$objective
        msl$par <- msl$solution
      }
      ## then we will optimize in the convex hull using constroptim(R). But we need a good starting point
    }
  }
  if (verbose) {
    if (inherits(object,"SLik")) {
      cat("* Summary ML: *\n")
      print(c(msl$par,"logL"=msl$value))
    } else if (inherits(object,"SLik_j")) {
      cat("* Summary ML: *\n")
      print(c(msl$par,"logL"=msl$value))
    } else if (inherits(object,"SLikp")) {
      cat("* Summary MQ: *\n")
      print(c(msl$par,"tailp"=msl$value))
    }
  }
  if (length(lower)==1) names(msl$par) <- names(lower)
  #
  hess <- hessian(objectivefn, msl$par)
  rnge <- diag(x=sqrt(upper-lower), ncol=length(upper))
  hess <- rnge %*% hess %*% rnge
  chk_identif <- .check_identifiables(hess,parnames=names(msl$par))
  if (length(chk_identif$non)) message(chk_identif$message)
  #
  object$MSL <- list2env(list(MSLE=msl$par, maxlogL=msl$value, hessian=hess, updated_from_prof=(! is.null(prev_init))))
  if (inherits(object$fit,"HLfit")) object$MSL$predVar <- attr(predict(object,newdata=msl$par,variances=list(linPred=TRUE,dispVar=TRUE)),"predVar")
  # CIs
  prevmsglength <- 0L
  if(CIs) {
    locverbose <- (verbose && ! inherits(object$fit,"HLfit")## for HLfit object, printing is later
                   && optim_time>3) # guess when it is useful to be verbose from the time to find the maximum 
    if (locverbose) {
      prevmsglength <- .overcat("Computing confidence intervals...\n", 0L) 
    } 
    object$CIobject <- list2env(.allCIs(object,verbose=locverbose, level=level))  ## may be NULL
  } 
  # MSEs computation
  if (inherits(object$fit,"HLfit")) {
    object$RMSEs <- list2env(list(RMSEs=.RMSEwrapper(object,verbose=FALSE),
                                  warn=NULL))
  } else {
    if (eval_RMSEs) {
      if  (verbose) .overcat("Computing RMSEs... (may be slow)\n",prevmsglength)
      object$RMSEs <- list2env(list(RMSEs=.RMSEwrapper.SLik_j(object,cluster_args=cluster_args,verbose=FALSE),
                                    warn=NULL))
    } 
  }
  if (is.null(object$logLs$cumul_iter)) {
    object$logLs$cumul_iter <- 1L ## adds a column to the data frame; needed soon for .par_RMSEwrapper() -> . -> profile.SLik()
    attr(object$logLs,"n_first_iter") <- nrow(object$logLs)
  }
  # par_RMSEs computation requires prior $CIobject information and $RMSEs information:
  object$par_RMSEs <- list2env(list(par_RMSEs=.par_RMSEwrapper(object,verbose=FALSE)))
  if  (verbose) {
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
  object$`Infusion.version` <- packageVersion("Infusion")
  invisible(object)
}

