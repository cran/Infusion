# called by SLRT with bootstrap -> .ecdf_2lr(),  and by latint()
.optim_mlogL_newobs <- function(object, 
                            givenpars=NULL, # optional profile parameters, but necessary only if their value differ from the MSLE
                            newobs, # 1-row matrix in projected summstats space 
                            init=NULL,
                            lower, upper, # in restricted space if profile
                            which,
                            log=TRUE) {
  if (is.null(init)) init <- .safe_init(object=object, given=givenpars, newobs = newobs, plower=lower,pupper=upper)
  template <- object$MSL$MSLE
  if (!is.null(givenpars)) template[names(givenpars)] <- givenpars 
  if ( inherits(object$jointdens,"dMixmod")) { 
    thr_info <- .get_thr_info(object)
    prof_mlogLfn_Dsim <- function(parmv) { 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .get_dens_from_GMM.Rmixmod(object=object, X=fullpar, tstat.obs=newobs, 
                                   log=log, which=which, thr_info=thr_info)
    }
  } else if (inherits(object$jointdens,"Mclust")) { # code presumably not tested (_F I X M E_)
    thr_info <- .get_thr_info(object)
    prof_mlogLfn_Dsim <- function(parmv) { 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .predict_SLik_j_mclust(object=object, X=fullpar, tstat.obs=newobs, 
                               log=log, which=which, thr_info=thr_info)
    }
  } else {
    prof_mlogLfn_Dsim <- function(parmv) { 
      fullpar <- template
      fullpar[names(init)] <- parmv
      - .get_densv(object=object, X=fullpar, tstat.obs=newobs, log=log, which=which)}
  }
  if ( ! is.null(constr_crits <- object$constr_crits)) {
    neg_ineq_constrfn <- function(parmv) {
      fullpar <- template
      fullpar[names(init)] <- parmv
      as.numeric(eval(constr_crits, envir = as.list(fullpar)))
    }
  } else neg_ineq_constrfn <- NULL
  
  .safe_optim(init, objfn=prof_mlogLfn_Dsim, lower=lower,upper=upper, 
              neg_ineq_constrfn=neg_ineq_constrfn, LowUp=list(), verbose=FALSE, object=object)
}

.ecdf_2lr <- function(object, 
                      BGP, # bootstrap-generating parameters
                      h0_pars=NULL, 
                      nsim=100L,
                      nb_cores=1L,
                      sim_object=TRUE, # simulate the fit, not the SGP 
                      #  To simulate from the SGP, set it to FALSE
                      #  [which may be passed as private option through the dots of SLRT() -> get_LRboot()]
                      cluster_args=NULL, # parallelisation is possible, but PSOCK is not useful in low-dimensional cases
                      ... 
) { # using gaussian-mixture approximation of distribution of 'data'
  # data simulated given 'BGP' for LRT of 'h0' consistent with BGP
  if (sim_object) {
    statempdist_h0 <- simulate(object, given=BGP, nsim=nsim)  
  } else {
    Simulate <- get_from(object,"Simulate")
    if (is.null(Simulate)) stop("No 'Simulate' function specified: only parameter points are returned.")
    if (is.null(dim(BGP))) BGP <- t(BGP)
    # reftable_cluster_args <- .lookup(cluster_args,try_in="reftable")
    
    statempdist_h0 <- add_reftable(Simulate=Simulate, 
                                   parsTable=as.data.frame(BGP[rep(1L, nsim),]),
                                   # verbose=verbose$progress_bars, 
                                   control.Simulate=get_from(object,"control.Simulate"),
                                   # cluster_args=reftable_cluster_args, 
                                   packages=get_from(object,"packages")$add_simulation, 
                                   env=get_from(object,"env"),
                                   cl_seed=.update_seed(object))
    if ( ! is.null(projectors <- object$projectors)) {
      # project_methodArgs <- .lookup(cluster_args, try_in="project")
      # project_methodArgs[names(methodArgs)] <- methodArgs # user-level methodArgs override cluster_args,
      missing_forests <- .missing_forests(object)
      if (any(missing_forests)) {
        # either forest is missing (as per ranger(., write.forest=FALSE))
        # or I removed one of its bulky elements.
        message("'forest' information not available from some projector(s). Regenerating it...")
        reftable_raw <- .get_reft_raw(object)
        sub_reftable <- reftable_raw[object$proj_trainset,]
        for (st in names(which(missing_forests))) {
          projectors[[st]] <- .update_projector.SLik_j(
            object, proj=st, reftable_raw=sub_reftable, methodArgs = list(write.forest=TRUE))
        }
      } 
      statempdist_h0 <- .project_reftable_raw(statempdist_h0, projectors=projectors, 
                                              methodArgs=list(), # need to control parall independently for simuls and for project. 
                                              use_oob=FALSE, # bc only newdata are projected and they are not in training set of re-used predictors
                                              ext_projdata=object$projdata
      )
      statempdist_h0 <- statempdist_h0[,object$colTypes$statNames, drop=FALSE]
    }
  }
  if (length(h0_pars)) { # profile
    profiledOutpars <- setdiff(names(BGP), h0_pars)
    prof_init_pars <- BGP[profiledOutpars] # true DGP should be a good starting value
    prof_lower <- object$lower[profiledOutpars]
    prof_upper <- object$upper[profiledOutpars]
  }
  
  given <- BGP[h0_pars]
  # locfn is like a local version of summLik since both call .optim_mlogL_newobs() internally.
  if (length(given)==1L && length(profiledOutpars)) { ## special case where 
    ## .safe_init([base_inits=NULL, length(given)=1]) => .get_init_from_hull() 
    ## is used as much as possible (.inits_from_postdens() called too).
    ## This special case is the most important, for bootstraps in SLRT in particular.
    ## The alternative code may be valid in this case too but has not been similarly tuned.
    
    init_nsim <- 1000L*length(object$lower) # explicit nsim to force call 
    # of .inits_from_postdens() despite newobs+given. # => (v2.1.184 drastic effect on bootstrap results)
    
    .locfn_MLsol_with_LR <- function(it, object) { # other arguments are in the environment of the function, 
      # which seems sufficient. For clarity 'object' arg receives locobject in MAF version. 
      newobs <- statempdist_h0[it,, drop=FALSE] # in projected summstats space
      
      ## "ML" for the simulated newobs (not prof lik) => solution for all params in object$lower
      #
      init_h1 <- .safe_init(object, # accounts for constr_crits...
                            newobs=newobs, plower=object$lower, pupper=object$upper,
                            nsim=init_nsim) 
      # Using the hessian of summLik to find init_h1 did not give good results.
      optr_h1 <- .optim_mlogL_newobs(object, newobs=newobs,
                                     init=unlist(init_h1), which="safe",
                                     lower=object$lower, upper=object$upper)
      
      ## proflik for tested hypo
      init_h0 <- .safe_init(object, given=given, profiledNames=profiledOutpars, 
                            newobs=newobs, plower=prof_lower, pupper=prof_upper,
                            MSLE=optr_h1$solution, nsim=init_nsim) 
      optr_h0 <- .optim_mlogL_newobs(object, 
                                     givenpars=given, # would be useful only for non-default (object$MSL$MSLE) values
                                     newobs=newobs, which="safe",
                                     init=unlist(init_h0)[profiledOutpars], 
                                     lower=prof_lower, upper=prof_upper)
      
      LRstat <- - optr_h1$objective + optr_h0$objective # logL_h1 - logL_h0
      if (LRstat<0) { # This case reappears in the 13-param demo-genetic simulation tests,
        # so addition v2.1.196.1 (sanitize in 196.3... but what about constr_crits?):
        init_h1 <- c(optr_h0$solution, given)[object$colTypes$fittedPars]
        init_h1 <- .sanitize_optim_solution(init_h1, lower=object$lower, upper=object$upper)
        optr_h1 <- .optim_mlogL_newobs(object, newobs=newobs,
                                       init=init_h1, which="safe",
                                       lower=object$lower, upper=object$upper)
        LRstat <- - optr_h1$objective + optr_h0$objective # logL_h1 - logL_h0  
      }
      MLsol_with_LR <- optr_h1["solution"]
      MLsol_with_LR$LRstat <- LRstat
      MLsol_with_LR
    }
  } else {
    ## alternative where there is no .get_init_from_hull() procedure to call. 
    ## .safe_init([base_inits=NULL, length(given)!=1]): inits_from_postdens() is still called.
    .locfn_MLsol_with_LR <- function(it, object) { # other arguments are in the environment of the function, which seems sufficient.
      newobs <- statempdist_h0[it,, drop=FALSE] # in projected summstats space
      
      if (length(profiledOutpars)) {
        # Maximization of likelihood wrt profiledOutpars (for the simulated data):
        # it should be like running profile(object, value=BGP[h0_pars]) except that new data are used hence not the lik surf of the 'object'
        optr_h0 <- .optim_mlogL_newobs(object, 
                                       givenpars=given, # would be useful only for non-default (object$MSL$MSLE) values
                                       newobs=newobs, which="safe",
                                       init=NULL, # use .safe_init() internally; alternative to prof_init_pars 
                                       lower=prof_lower, upper=prof_upper)
        logL_h0 <- - optr_h0$objective
        init_h1 <- c(optr_h0$solution, given)[object$colTypes$fittedPars]
        init_h1 <- .safe_init(object=object, newobs = newobs, plower=object$lower,
                              pupper=object$upper, more_inits=init_h1)
      } else {
        logL_h0 <- .get_densv(object, X=BGP, tstat.obs=newobs, log=TRUE, which="safe")
        init_h1 <- NULL # .safe_init() will be called internally by .optim_mlogL_newobs() 
      }
      
      ## "ML" for the simulated newobs (not prof lik) => solution for all params in object$lower
      optr_h1 <- .optim_mlogL_newobs(object, newobs=newobs,
                                     init=init_h1, which="safe",
                                     lower=object$lower, upper=object$upper)
      LRstat <- - optr_h1$objective - logL_h0 
      MLsol_with_LR <- optr_h1["solution"]
      MLsol_with_LR$LRstat <- LRstat
      MLsol_with_LR
    }
  }
  
  if (nb_cores==1L) { # alternative works also in this case but with some (moderate) overhead.
    resu <- vector("list", nsim)
    for (it in seq(nsim)) resu[[it]] <- .locfn_MLsol_with_LR(it, object=object)
  } else {
    if (inherits(object$jointdens,"MAF") && 
        .set_cluster_type(cluster_args, nb_cores)$spec >1L) {
      ## multiple processes may then use gpu... hmmm. 
      ## Provide the MAF objects:
      if (is.null(object$load_MAFs_info$pwd)) {
        stop("save_MAFs() should have been run so that child processes can load_MAFs().")
      } 
      torch_device <- .get_py_MAF_handle()$device$type
      # locfn() will use 'object' from its defining envir
      # .locfn() has an 'object' argument so that 'object' can be modified locally (1st case)
      locfn <- function(it) {
        config_mafR(torch_device=torch_device)
        load_MAFs_info <- object$load_MAFs_info
        locobject <- load_MAFs(object,prefix=load_MAFs_info$prefix, ext=load_MAFs_info$ext)
        .locfn_MLsol_with_LR(it,object=locobject)
      }
      fit_env <- list(torch_device=torch_device)
    } else  { 
      fit_env <- list()
      locfn_MLsol_with_LR <- function(it) .locfn_MLsol_with_LR(it,object=object) 
    }
    seq_it_as_1row_mat <- seq(nsim)
    dim(seq_it_as_1row_mat) <- c(1L,nsim)
    if (nb_cores>1L) object <- ..shrink(object)
    resu <- spaMM::dopar(seq_it_as_1row_mat, fn = locfn_MLsol_with_LR, 
                         fit_env = fit_env, nb_cores=nb_cores, ...)
  }
  resu <- do.call(rbind, resu)
  boot_ests <- do.call(rbind,resu[,"solution"])
  list(ecdf_2lr=2*unlist(resu[,"LRstat"], use.names = FALSE), 
       ecdf_t=boot_ests[,h0_pars, drop=FALSE])
}


.call4print <- function (x, max.print=10L) { 
  x <- as.list(x)
  xnames <- names(x) 
  for (it in seq_along(x)) {
    xx <- x[[it]]
    if (xnames[it]=="boot.out") {
      for (jt in seq_along(xx)) {
        xxx <- xx[[jt]]
        if (is.numeric(xxx) & length(xxx)> max.print) {
          xx[[jt]] <- capture.output(str(xxx))[[1]]
          names(xx)[jt] <- paste0("str(long numeric '",names(xx[jt]),"')")
        }
      }
      x[[it]] <- as.call(c(quote(list),xx)) # hack so that boot:::print.bootci 
                    #  does not drop the names of the modified 'boot.out' list.
    } else {
      if (is.numeric(xx) & length(xx)> max.print) {
        x[[it]] <- capture.output(str(xx))
        names(x)[it] <- paste0("str(long numeric '",names(x[it]),"')")
      }
    }
  }
  as.call(x)
}


get_LRboot <- function(object, h0_pars=NULL, nsim=100L, reset=TRUE, 
                       BGP=object$MSL$MSLE,  # bootstrap-generating parameters
                       which="ecdf_2lr", bootCI.args=list(type="perc", conf = 0.95),
                       ...) {
  # no local control of the seed. No use for a seed in the \dots. A standard set.seed should provide control. 
  if ( ! inherits(object,"SLik_j")) stop("Bootstrap objects for LRTs are implemented only for objects of class 'SLik_j'.")
  
  if (is.null(object$bootLRTenv)) {
    errmess <- paste("SLik_j object has no 'bootLRTenv' element.\n",
                     "Was it created using an old version of Infusion?\n.",
                     "If so, running MSL() on it may fix this issue.")
    stop(errmess)
  }
  # : created := last run of MSL(), but run possibly not at user-level so I cannot mention it
  
  if (is.null(h0_pars)) {
    element <- "<full>"
  } else element <- paste0(sort(h0_pars),collapse=",")
  bootreps_list <- object$bootLRTenv$bootreps_list
  #
  if (reset || is.null(resu <- bootreps_list[[element]])) {
    resu <- .ecdf_2lr(object, BGP=BGP, h0_pars=h0_pars, nsim=nsim, ...)
    if (length(h0_pars)==1L  && nsim>2L ) {
      bootwarn <- NULL
      withCallingHandlers(
        bootCI <- do.call("boot.ci", c(
          list(boot.out = list(t=resu$ecdf_t, t0=BGP[h0_pars], R=nsim, sim="parametric")), 
          bootCI.args)),
        warning=function(w) {
          bootwarn <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        })
      if (length(bootwarn)) 
        warning(paste("In get_LRboot() -> boot.ci():", bootwarn), call.=FALSE)
      bootCI$call <- .call4print(bootCI$call)
      resu$bootCI <- bootCI
    }
    bootreps_list[[element]] <- resu
    assign("bootreps_list",bootreps_list,envir=object$bootLRTenv)
  } else if (length(resu)!=nsim) message("Saved object used as reset=FALSE: 'nsim' argument ignored.")
  if ( ! is.null(which)) {
    resu <- resu[["which"]]
  }
  resu
}

.fit_2lr <- function(chi2_obs, ecdf_2lr, df) {
  if (df>1L) stop(".fit_2lr() does not handle df>1") # we're working on correcting 
  #   the root of the sum of squares of 'df' nominally gaussian RVs...
  # so the 'chi" (not chi2) distrib... CRAN task view -> package Runuran if we really need it.
  if (chi2_obs <0.01) return(chi2_obs)
  r_obs <- sqrt(chi2_obs)
  r <- sqrt(pmax(0.01,ecdf_2lr))
  objfn <- function(u.r) {
    rn <- r+log.abs.u.r/r
    - sum(dnorm(rn,mean = 0,sd = 1,log = TRUE))
  }
  rfit <- .safe_opt(init = 1, objfn = objfn,lower = -Inf, upper=Inf, 
                            LowUp=list(), verbose=FALSE)
  log.abs.u.r <- rfit$solution
  rs <- r_obs+log.abs.u.r/r_obs # adjusted r
  rs^2
}


SLRT <- function(object, h0, nsim=0L, BGP=NULL, type="perc", level=0.95, nsteps=10L, 
                 variants=NULL,
                 ...) {
  # no local control of the seed. No use for a seed in the \dots. A standard set.seed should provide control. 
  resu <- list()
  fittedPars <- object$colTypes$fittedPars
  chk <- setdiff(names(h0), fittedPars)
  if (length(chk)) stop("'h0' parameters not all in fitted parameter names: this case is not handled.")
  profiledOutPars <- setdiff(fittedPars,names(h0))
  MSL <- object$MSL
  profNames <- names(h0)
  if ( any(c(h0-object$lower[profNames],object$upper[profNames]-h0)<0)) {
    warning("'h0' is out of the effectively sampled parameter ranges.")
  }
  if (length(profiledOutPars)) {
    profLik <- profile(object, value=h0, which="safe") # without init, so that the the convex-hull method will be used. 
    if (nsteps>1L) {  ## not tested for any length(h0)>1L
      # Not as efficient individually as the convex-hull method, but still helpful in some cases:
      # cf comparison of SLRTs B_13from17 in summaries_list.1_200.v2.1.112.rda  vs ... 2.1.111.rda
      ## not tested for any length(h0)>1L
      # Not as efficient individually as the convex-hull method, but still helpful in some cases:
      # cf comparison of SLRTs B_13from17 in summaries_list.1_200.v2.1.112.rda  vs ... 2.1.111.rda
      #       #prof_xs <- seq(MSL$MSLE[profNames], h0, length.out=nsteps)
      prof_xs <- do.call(cbind,
                         lapply(profNames, function(v) seq(MSL$MSLE[v], h0[v], length.out=nsteps)) )
      templateh0 <- h0
      prof_ys <- numeric(nsteps)
      prof_ys[1] <- MSL$maxlogL
      if (TRUE) {
        init <- MSL$MSLE
        for (it in 2L:(nsteps-1L)) {
          templateh0[] <- prof_xs[it,]
          p <- profile(object, value=templateh0, init=init, which="safe")
          prof_ys[it] <- p[]
          init <- .sanitize_optim_solution(attr(p,"solution"), object$lower, object$upper) 
        }
      } else { # more intensive try. Not convincing. 
        template_init <- MSL$MSLE
        plower <- object$LOWER[profiledOutPars]
        pupper <- object$UPPER[profiledOutPars]
        init_nsim <- (1000L*length(object$lower))%/% nsteps
        for (it in 2L:(nsteps-1L)) {
          templateh0[] <- prof_xs[it,]
          init <- .safe_init(object,given=templateh0,plower = plower,pupper = pupper,
                             more_inits = template_init, nsim=init_nsim)
          p <- profile(object, value=templateh0, init=init, which="safe")
          prof_ys[it] <- p[]
          template_init[profiledOutPars] <- attr(p,"solution")
        }
        init <- .safe_init(object,given=templateh0,plower = plower,pupper = pupper,
                           more_inits = template_init, nsim=init_nsim)
      }
      prof_ys[nsteps] <- proflast <- profile(object, value=h0, init=init, which="safe")
      if (proflast> profLik) profLik <- proflast # cannot use max() bc it would lose the attribute
    }  
  } else profLik <- predict(object, newdata=h0, which="safe", constr_tuning = Inf)
  profiledOutvals <- attr(profLik,"solution")

  LRTori <- 2*(MSL$maxlogL - profLik[[1]]) # [[1]] remove sticky attributes
  if (MSL_updating <- (LRTori< -1e-4)) {
    init <- MSL$MSLE
    init[profiledOutPars] <- profiledOutvals 
    init[names(h0)] <- h0
    newmax <- MSL(object, init=init, eval_RMSEs=FALSE, CIs=FALSE) # with a NEW $MSL environment 
    for (st in ls(newmax$MSL)) object$MSL[[st]] <- newmax$MSL[[st]] ## but keeping the environment (~pointer) unchanged.
    LRTori <- 2*(object$MSL$maxlogL - profLik)
  }
  df <- length(h0)
  resu$basicLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=1-pchisq(LRTori,df=df))
  if (nsim) {
    if (is.null(BGP)) {
      if (length(profiledOutPars)) {
        # use nuisance-parameter estimates from the 'profLik' profile for h0 on the original data
        profiledoutpars <- attr(profLik,"solution") # profile.SLik_j(object,value=h0, return.optim=TRUE)$solution
        BGP <- object$MSL$MSLE
        BGP[names(h0)] <- h0
        BGP[names(profiledoutpars)] <- profiledoutpars
      } else {
        BGP <- h0
      }
    }
    h0_pars <- names(h0)
    # bootstrap  simul under the null hypo (where nuisance param have been re-estimated)
    LRboot <- get_LRboot(object, h0_pars=h0_pars, nsim=nsim, BGP=BGP, which=NULL, 
                         bootCI.args= list(type=type,level=level), ...) 
    ecdf_2lr <- LRboot$ecdf_2lr
    meanLR <- mean(ecdf_2lr)
    LRTcorr <- LRTori *df / meanLR
    resu$BartBootLRT <- data.frame(chi2_LR=LRTcorr,df=df,p_value=1-pchisq(LRTcorr,df=df))
    resu$bootCI <- LRboot$bootCI
    if ("gamma" %in% variants) { # private option
      dat <- data.frame(y=ecdf_2lr)
      gfit <- fitme(y ~1, family=Gamma(log), data=dat )
      mu <- exp(fixef(gfit))
      gammap <- 1-pgamma(LRTori, shape= 1 / gfit$phi, scale=mu*gfit$phi)
      resu$gammaCorrLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=gammap)
    }
    rawPvalue <- (1+sum(ecdf_2lr>=LRTori))/(nsim+1L) ## DavisonH, p.141
    resu$rawBootLRT <- data.frame(chi2_LR=LRTori,df=df,p_value=rawPvalue)
    if ("safeBartBootLRT" %in% variants) {
      if (meanLR < df) {
        resu$safeBartBootLRT <- resu$rawBootLRT
      } else resu$safeBartBootLRT <- resu$BartBootLRT
    }
    if ("pstarBoot" %in% variants) {
      pstarBoot <- .fit_2lr(chi2_obs = LRTori,ecdf_2lr = ecdf_2lr, df=df) 
      resu$pstarBoot <- data.frame(chi2_LR=pstarBoot,df=df,p_value=1-pchisq(pstarBoot,df=df))
    }
  }
  attr(resu,"MSL_updated") <- MSL_updating
  resu
}

