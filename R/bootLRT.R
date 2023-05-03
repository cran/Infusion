.optim_mlogL_newobs <- function(object, 
                            givenpars=NULL, # optional profile parameters, but necessary only if their value differ from the MSLE
                            newobs, # 1-row matrix in projected summstats space 
                            init=NULL,
                            lower, upper, # in restricted space if profile
                            log=TRUE) {
  if (is.null(init)) init <- .safe_init(object=object, given=givenpars, newobs = newobs, plower=lower,pupper=upper)
  template <- object$MSL$MSLE
  if (!is.null(givenpars)) template[names(givenpars)] <- givenpars 
  if (inherits(object$jointdens,"Mclust")) { # code not tested (_F I X M E_)
    prof_mlogLfn_Dsim <- function(parmv) { 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .predict_SLik_j_mclust(object=object, newdata=fullpar, tstat.obs=newobs, log=log, which="")
    }
  } else {
    prof_mlogLfn_Dsim <- function(parmv) { 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .get_dens_from_GMM(object=object, X=fullpar, tstat.obs=newobs, log=log, which="lik")}
  }
  .safe_opt(init, objfn=prof_mlogLfn_Dsim, lower=lower,upper=upper, LowUp=list(), verbose=FALSE)
}

.ecdf_2lr <- function(object, 
                      BGP, # bootstrap-generating parameters
                      h0_pars=NULL, 
                      nsim=100L,
                      nb_cores=1L,
                      ... # parallelisation is possible, but PSOCK is not useful in low-dimensional cases
                      ) { # using gaussian-mixture approximation of distribution of 'data'
  # data simulated given 'BGP' for LRT of 'h0' consistent with BGP
  if (inherits(object$gllimobj,"gllim")) {
    statempdist_h0 <- .gllim.condsimul.stats(object$gllimobj, RGPpars=BGP, size=nsim, colTypes=object$colTypes, cbind.=FALSE)
  } else {
    statdens_h0 <- .conditional_Rmixmod(object$jointdens, given=BGP, expansion=1) # stat dens|ML parameter estimates
    statempdist_h0 <- .simulate.MixmodResults(statdens_h0, nsim=1L, size=nsim, drop=TRUE) # directly in projected space
  }
  if (length(h0_pars)) { # profile
    profiledOutpars <- setdiff(names(BGP), h0_pars)
    prof_init_pars <- BGP[profiledOutpars] # true DGP should be a good starting value
    prof_lower <- object$lower[profiledOutpars]
    prof_upper <- object$upper[profiledOutpars]
  }
  locfn <- function(it) { # other arguments are in the environment of the function, which seems sufficient.
    newobs <- statempdist_h0[it,, drop=FALSE] # in projected summstats space
    # Maximization of likelihood wrt full param (for the simulated data):
    #init <- BGP
    if (length(h0_pars)) {
      # Maximization of likelihood wrt profiledOutpars (for the simulated data):
      # it should be like running profile(object, value=BGP[h0_pars]) except that new data are used hence not the lik surf of the 'object'
      #init <- prof_init_pars
      opt_prof_mlogL_Dsim <- .optim_mlogL_newobs(object, 
                                                 givenpars=BGP[h0_pars], # would be useful only for non-default (object$MSL$MSLE) values
                                                 newobs=newobs,
                                                 init=NULL, # use .safe_init() internally; alternative to prof_init_pars 
                                                 lower=prof_lower, upper=prof_upper)
      logL_h0 <- - opt_prof_mlogL_Dsim$objective
      init_h1 <- c(opt_prof_mlogL_Dsim$solution, BGP[h0_pars])[object$colTypes$fittedPars]
      init_h1 <- .safe_init(object=object, newobs = newobs, plower=object$lower,pupper=object$upper, more_inits=init_h1)
    } else {
      logL_h0 <- .get_dens_from_GMM(object, X=BGP, tstat.obs=newobs, log=TRUE, which="lik")
      init_h1 <- NULL # .safe_init() will be called internally by .optim_mlogL_newobs() 
    }
    opt_mlogL_Dsim <- .optim_mlogL_newobs(object, newobs=newobs,
                                          init=init_h1, 
                                          lower=object$lower, upper=object$upper)
    resu <- - opt_mlogL_Dsim$objective - logL_h0
    # if (resu<0) browser()
    resu
  }
  if (nb_cores==1L) { # alternative works also in this case but with some (moderate) overhead.
    ecdf_lr <- numeric(nsim)
    for (it in seq(nsim)) ecdf_lr[it] <- locfn(it)
  } else {
    seq_it_as_1row_mat <- seq(nsim)
    dim(seq_it_as_1row_mat) <- c(1L,nsim)
    ecdf_lr <- spaMM::dopar(seq_it_as_1row_mat, fn = locfn, fit_env = NULL, nb_cores=nb_cores, ...)
  }
  2* ecdf_lr # 2* !
}

get_LRboot <- function(object, h0_pars=NULL, nsim=100L, reset=TRUE, 
                       BGP=object$MSL$MSLE,  # bootstrap-generating parameters
                       ...) {
  if ( ! inherits(object,"SLik_j")) stop("Bootstrap objects for LRTs are implemented only for objects of class 'SLik_j'.")
  if (is.null(object$bootLRTenv)) stop("SLik_j object has no 'bootLRTenv' element (was it created using an old version of Infusion?)")
  if (is.null(h0_pars)) {
    element <- "<full>"
  } else element <- paste0(sort(h0_pars),collapse=",")
  bootreps_list <- object$bootLRTenv$bootreps_list
  if (reset || is.null(resu <- bootreps_list[[element]])) {
    resu <- .ecdf_2lr(object, BGP=BGP, h0_pars=h0_pars, nsim=nsim, ...)
    bootreps_list[[element]] <- resu
    assign("bootreps_list",bootreps_list,envir=object$bootLRTenv)
  } else if (length(resu)!=nsim) message("Saved object used as reset=FALSE: 'nsim' argument ignored.")
  return(resu)
}

SLRT <- function(object, h0, nsim=0L, BGP=NULL, ...) {
  resu <- list()
  fittedPars <- object$colTypes$fittedPars
  chk <- setdiff(names(h0), fittedPars)
  if (length(chk)) stop("'h0' parameters not all in fitted parameter names: this case is not handled.")
  profiledOutPars <- setdiff(fittedPars,names(h0))
  if (length(profiledOutPars)) {
    LRTori <- 2*(object$MSL$maxlogL - profile(object, value=h0, which="safe"))
    # prof <-  profile(object, value=h0, which="safe", return.optim=TRUE)
    # LRTori <- 2*(object$MSL$maxlogL + prof$objective) # '+' bc objfn is mlogL
    # if (LRTori<0) {
    #   MSL(object, init= c(h0, prof$solution))
    # }
    df <- length(h0)
    h0_pars <- names(h0)
  } else {
    LRTori <- 2*(object$MSL$maxlogL - predict(object, newdata=h0, which="safe"))
    df <- length(fittedPars)
    h0_pars <- NULL
  }
  resu$basicLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=1-pchisq(LRTori,df=df))
  if (nsim) {
    if (is.null(BGP)) {
      if (length(profiledOutPars)) {
        profiledoutpars <- profile.SLik_j(object,value=h0, return.optim=TRUE)$solution
        BGP <- object$MSL$MSLE
        BGP[names(h0)] <- h0
        BGP[names(profiledoutpars)] <- profiledoutpars
      } else {
        BGP <- h0
      }
    }
    LRboot <- get_LRboot(object, h0_pars=h0_pars, nsim=nsim, BGP=BGP, ...) 
    meanLR <- mean(LRboot)
    LRTcorr <- LRTori *df / meanLR
    resu$BartBootLRT <- data.frame(chi2_LR=LRTcorr,df=df,p_value=1-pchisq(LRTcorr,df=df))
    if ("gamma" %in% .Infusion.data$options$SLRTopts) { # private option
      dat <- data.frame(y=LRboot)
      gfit <- fitme(y ~1, family=Gamma(log), data=dat )
      mu <- exp(fixef(gfit))
      gammap <- 1-pgamma(LRTori, shape= 1 / gfit$phi, scale=mu*gfit$phi)
      resu$gammaCorrLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=gammap)
    }
    rawPvalue <- (1+sum(LRboot>=LRTori))/(length(LRboot)+1) ## DavisonH, p.141
    resu$rawBootLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=rawPvalue)
    if (meanLR < df) {
      resu$safeBartBootLRT <- resu$rawBootLRT
    } else resu$safeBartBootLRT <- resu$BartBootLRT
  }
  resu
}
