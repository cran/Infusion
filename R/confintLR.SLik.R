# generic: confint(object, parm, level = 0.95, ...)

## same methods for SLikp and SLik because the code differs only on one line...
confint.SLik <- function(object, parm, ## parm is the parameter which CI is sought 
                         level=0.95, verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE),...) {
  .confintAll(object=object, parm=parm, ## parm is the parameter which CI is sought 
             givenmax = object$MSL$maxlogL,
             level= - qchisq(level,df=1)/2, ## df=1 for 1D profile; /2 to compare to logLR rather than 2 logLR  
             verbose=verbose,fixed=fixed,which=which,...)
}

# experimental; proved insuficient
.newdata_from_pardens_clu_means <- function(object, given) {
  if (inherits(object$jointdens,"Mclust")) { # code not tested (_F I X M E_)
    conddens <- .conditional_mclust(object$pardens,fittedPars=object$colTypes$fittedPars ,given=given, expansion=1)
    means <- conddens$parameters$mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  } else {
    conddens <- .conditional_Rmixmod(object$pardens, given=given, expansion=1) 
    solve_t_chol_sigma_list <- lapply(conddens@parameters["variance"], .solve_t_cholfn)
    means <- conddens@parameters@mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  }
  fittedPars <- object$colTypes$fittedPars
  newdata <- matrix(NA,nrow = 1L, ncol=length(fittedPars), dimnames=list(NULL,fittedPars))
  newdata <- newdata[rep(1,nrow(means)),, drop=FALSE]
  newdata[,colnames(means)] <- means
  newdata[,names(given)] <- given
  newdata
}

.init_params_from_pardens <- function(object, 
                                      given, # always conditional to some values (=> subset of parameter values for the profiles)
                                      profiledNames, # in the default method (max_conddens=FALSE) this only serves to select the returned parameters
                                      plower, pupper,
                                      newobs=NULL) {
  if (inherits(object$jointdens,"Mclust")) { # code not tested (_F I X M E_)
    condpardens <- .conditional_mclust(object$pardens,fittedPars=object$colTypes$fittedPars ,given=given, expansion=1)
    means <- condpardens$parameters$mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  } else if (inherits(object$gllimobj,"gllim")) { 
    # we want the means of the conddens given 'given' , a nclu * length(fittedPars) matrix
    means <- .condProfoutParMeans.gllim(object$gllimobj, fittedPars=object$colTypes$fittedPars, given=given, expansion=1) 
  } else {
    condpardens <- .conditional_Rmixmod(object$pardens, given=given, expansion=1) 
    means <- condpardens@parameters@mean # use predictions in mean of each Gaussian component to define initial value of search of maximum of the density
  }
  fittedPars <- object$colTypes$fittedPars
  DGPparms <- matrix(NA,nrow = 1L, ncol=length(fittedPars), dimnames=list(NULL,fittedPars))
  DGPparms <- DGPparms[rep(1,nrow(means)),, drop=FALSE]
  DGPparms[,colnames(means)] <- means
  DGPparms[,names(given)] <- given
  if (is.null(newobs)) {
    best_clu <- which.max(predict(object, newdata=DGPparms, which="safe"))
  } else best_clu <- which.max(summLik(object, parm=DGPparms, data=newobs, which="safe"))
  resu <- DGPparms[best_clu,profiledNames]
  resu
}

# generic: profile(fitted,...), returns a # *log*L
profile.SLik <- function(fitted, value, fixed=NULL, return.optim=FALSE, init="default", which="safe", ...) {
  fixedPars <- names(fixed)   
  fittedPars <- fitted$colTypes$fittedPars
  if (! is.null(fixed)) fittedPars <- setdiff(fittedPars,fixedPars)
  fittedparamnbr <- length(fittedPars) 
  parm <- names(value)
  MLval <- fitted$MSL$MSLE[parm]
  if (is.na(MLval)) {stop(paste("'",parm,"' appears to be an incorrect parameter name. Check 'parm'."))}
  lowval <- fitted$lower[parm]
  hival <- fitted$upper[parm]
  v <- fitted$MSL$MSLE; v[names(v)] <- NA ## create template 
  lowval <- lowval + 0.002 * (MLval - lowval)
  hival <- hival - 0.002 * (hival - MLval)
  if (fittedparamnbr == 1L) {
    v[parm] <- value
    if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
    resu <- predict(fitted,newdata=v, which=which)
  } else { 
    profiledNames <- names(fitted$lower)
    profiledNames <- setdiff(profiledNames, c(parm,fixedPars)) # [which( ! (profiledNames %in% parm))] 
    plower <- fitted$lower[profiledNames]
    pupper <- fitted$upper[profiledNames]
    v[parm] <- value 
    plogL <- function(pparv) {
      v[profiledNames] <- pparv
      if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
      return(( - predict(fitted,newdata=v, which=which))) ## removed log... (log=TRUE is the default)  
    }
    if (identical(init,"clu_means")) { # multiple optims, one from each clu mean
      newdata <- .newdata_from_pardens_clu_means(fitted,value)
      OPT <- function(init) {.safe_opt(init[profiledNames], plogL, lower=plower, upper=pupper, LowUp=list(), verbose=FALSE, ...)}
      optrS <- apply(newdata,1L, OPT)  
      optr <- optrS[[which.min(sapply(optrS,`[[`, i="objective"))]]
    } else { # single optim from best clu mean by logLik 
      if (identical(init,"default")) {
        if (inherits(fitted,"SLik_j")) { 
          init <- .safe_init(object = fitted, given=value, plower, pupper)
        } else init <- fitted$MSL$MSLE  
      } # else init must be user-provided vector including profiled pars 
      optr <- .safe_opt(init[profiledNames], plogL, lower=plower, upper=pupper, LowUp=list(), verbose=FALSE, ...)
    }
    if(return.optim) {
      optr$value <- - optr$objective 
      optr$par <- optr$solution
      return(optr)
    } else return(- optr$objective) 
  }
} # *log*L


.confintAll <- function(object, parm, ## parm is the parameter which CI is sought
                       givenmax=NULL,
                       #Bartlett=0L,
                       level, # value of 2 LR  
                       verbose=interactive(),fixed=NULL,which=c(TRUE,TRUE),...) {
  if (is.null(givenmax)) stop("The point estimates should be computed before using 'confint'.")
  fixedPars <- names(fixed)   
  if ( ! is.null(bootreps <- object$bootLRTenv$bootreps_list[[sort(parm)]])) {
    meanLR <- mean(bootreps)
    level <- level *meanLR / length(parm)  # correcting the threshold: inverse of correcting the LR stat
  }
  fittedPars <- object$colTypes$fittedPars
  if (! is.null(fixed)) fittedPars <- setdiff(fittedPars,fixedPars)
  fittedparamnbr <- length(fittedPars) 
  MLval <- object$MSL$MSLE[parm]
  if (is.na(MLval)) {stop(paste("'",parm,"' appears to be an incorrect parameter name. Check 'parm'."))}
  lowval <- object$LOWER[parm]
  hival <- object$UPPER[parm]
  ## FR->FR uniroot has no parscale agument 
  tol <- .Machine$double.eps^0.25* (hival-lowval)
  if (tol< .Machine$double.eps ) {
    warning(paste("UPPER-LOWER range for",parm,"is as narrow as",signif(hival-lowval,3),
                  "\n which is likely to cause problems in various numerical procedures.\n Consider using another parameter scale."))
  }
  v <- object$MSL$MSLE; v[names(v)] <- NA ## create template 
  lowval <- lowval + 0.002 * (MLval - lowval)
  hival <- hival - 0.002 * (hival - MLval)
  shift <- - givenmax - level ## opposite of predict value at the bound 
  if (fittedparamnbr == 1) {
    profiledNames <- c()
    objectivefn1D <- function(CIvarval) {
      v[parm] <- CIvarval
      if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
      predict(object,newdata=v, which="safe")+shift ## removed log... 
    }
    objectivefn <- objectivefn1D
  } else { 
    profiledNames <- names(object$lower)
    profiledNames <- setdiff(profiledNames, c(parm,fixedPars)) # [which( ! (profiledNames %in% parm))] 
    plower <- object$lower[profiledNames]
    pupper <- object$upper[profiledNames]
    objectivefnmultiD <- function(CIvarval,return.optim=FALSE) {
      #print(CIvarval)
      v[parm] <- CIvarval 
      plogL <- function(pparv) {
        v[profiledNames] <- pparv
        if (! is.null(fixed)) v[fixedPars] <- fixed[fixedPars]
        return(( - predict(object,newdata=v, which="safe"))) ## removed log...   
      }
      optr <- .safe_opt(object$MSL$MSLE[profiledNames], plogL, lower=plower, upper=pupper, LowUp=list(), verbose=FALSE)
      if(return.optim) {
        optr$value <- - optr$objective 
        optr$par <- optr$solution
        return(optr)
      } else return(- optr$objective+shift) ## then returns shifted value for uniroot
    }
    objectivefn <- objectivefnmultiD    
  }
  CIlo <- NA
  CIup <- NA
  fupper <- -level
  if (ML_at_bound <- (abs(lowval - MLval) < abs(MLval * 1e-08)) || !which[1]) {
    CIlo <- NA
    if (ML_at_bound && verbose) {
      mess <- paste0("Lower CI bound for ",parm," cannot be computed because point estimate is at lower bound of 'fitted' parameter range.")
      print(mess,quote=FALSE)
    }
  } else {
    flower <- objectivefn(lowval)
    if (is.na(flower)) {
      stop("From 'confint': 'flower' is NA")
    } else {
      if (flower < 0) {
        CIlo <- try((uniroot(objectivefn, interval = c(lowval, 
                                                       MLval), f.lower = flower, f.upper = fupper,tol=tol))$root, 
                    TRUE)
      } else if (verbose) {
        mess <- paste0("Lower CI bound for ",parm," cannot be computed as it appears to exceed the 'fitted' parameter range.")
        print(mess,quote=FALSE)
      }
    }
    if (inherits(CIlo,"try-error")) {
      CIlo <- NA
      errmsg <- paste("Lower CI bound for ", parm, " could not be computed (maybe out of sampled range)", 
                      sep = "")
      message(errmsg)
    } 
  }
  if (is.na(CIlo)) {
    lowfit <- NA
    lowerEstv <- NA
  } else {
    lowerEstv <- object$lower
    lowerEstv[parm] <- CIlo
    if (length(profiledNames)>0) {
      lowfit <- objectivefn(CIlo,return.optim=TRUE)
      lowerEstv[profiledNames] <- lowfit$par
    }      
  }
  ###
  flower <- fupper
  if ( ML_at_bound <- (abs(hival - MLval) < abs(MLval * 1e-08)) || !which[2]) {
    CIup <- NA
    if (ML_at_bound && verbose) {
      mess <- paste0("Upper CI bound for ",parm," cannot be computed because point estimate is at upper bound of 'fitted' parameter range.")
      print(mess,quote=FALSE)
    }
  } else {
    fupper <- objectivefn(hival)
    if (is.na(fupper)) {
      stop("From 'confint': 'fupper' is NA")
    } else {
      if (fupper < 0) {
        CIup <- try((uniroot(objectivefn, c(MLval, hival), 
                             f.lower = flower, f.upper = fupper,tol=tol))$root, 
                    TRUE)
      } else if (verbose) {
        mess <- paste0("Upper CI bound for ",parm," cannot be computed as it appears to exceed the 'fitted' parameter range.")
        print(mess,quote=FALSE)
      }
    }
    if (inherits(CIup,"try-error")) {
      CIup <- NA
      errmsg <- paste("Upper CI bound for ", parm, " could not be computed (maybe out of sampled range)", 
                      sep = "")
      message(errmsg)
    } 
  }
  if (is.na(CIup)) {
    upfit <- NA
    upperEstv <- NA
  } else {
    upperEstv <- object$upper
    upperEstv[parm] <- CIup
    if (length(profiledNames)>0) {
      upfit <- objectivefn(CIup,return.optim=TRUE)
      upperEstv[profiledNames] <- upfit$par
    } 
  }
  interval <- c(CIlo,CIup)
  names(interval) <- paste(c("low","up"),parm,sep=".")
  if (verbose) print(interval)
  invisible(list(lowerpar=lowerEstv,upperpar=upperEstv,interval=interval))
}


