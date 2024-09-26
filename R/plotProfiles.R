.MSL_update <- function(object) {
  message(paste0("A new likelihood maximum was found. The input object is being updated,\n",
                 "confidence intervals will have to be recomputed,\n",
                 "and plots will have to be redrawn."))
  # Discard info before calling MSL otherwise MSL will use it !
  if ( ! is.null(object$RMSEs)) {
    object$RMSEs$RMSEs <- NULL 
    object$RMSEs$warn <- paste("RMSE information for logL discarded as a new likelihood maximum was found by a profile-plot function.\n",
                               "Re-run MSL() to regenerate this information.")
  }
  if ( ! is.null(object$CIobject)) {
    for (st in names(object$CIobject)) object$CIobject[[st]] <- NULL
    object$CIobject$warn <-  paste("CI information discarded as a new likelihood maximum was found by a profile-plot function.\n",
                                   "Re-run MSL() to regenerate this information.")
  }
  if ( ! is.null(object$par_RMSEs)) {
    for (st in names(object$par_RMSEs)) object$par_RMSEs[[st]] <- NULL
    object$par_RMSEs$warn <-  paste("Intervals and RMSEs for parameters discarded as a new likelihood maximum was found by a profile-plot function.\n",
                                    "Re-run MSL() to regenerate this information.")
  }
  newmax <- MSL(object, eval_RMSEs=FALSE, CIs=FALSE) # with a NEW $MSL environment 
  for (st in ls(newmax$MSL)) object$MSL[[st]] <- newmax$MSL[[st]] ## but keeping the environment (~pointer) unchanged.
  # so the object$MSL$init_from_prof remains. Better to remove it now:
  object$MSL$init_from_prof <- NULL
} # only environmentswhere updated, so no return value

.calc_prof_coords <- function(st, object, fixedPars, 
                              safe, #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
                              MSLE, template, 
                              lower, upper, oldprofiles, gridSteps, type,
                              CIlevels=NULL) { 
  if (interactive()) {cat(paste("Computing profile for", st, "\n"))}
  prevmsglength <- 0L
  profiledOutPars <- setdiff(object$colTypes$fittedPars,c(st, fixedPars))
  constr_tuning <- FALSE
  objfn <- function(v, template) {
    template[profiledOutPars] <- v 
    - predict(object,template, which="safe", constr_tuning=constr_tuning)[1]
  }
  if ( ! is.null(constr_crits <- object$constr_crits)) {
    neg_ineq_constrfn <- function(v, template) {
      template[profiledOutPars] <- v 
      as.numeric(eval(constr_crits, envir = as.list(template)))
    }
  } else neg_ineq_constrfn <- NULL
  if ( ! safe) {
    gradfn <- .get_gradfn(object, init=MSLE[profiledOutPars])
    gradient <- function(v, template) {
      template[profiledOutPars] <- v 
      gradfn(template)
    }
  }
  init_from_neigh <- object$MSL$MSLE
  profil <- function(z, init=MSLE) {
    template[st] <- z
    ## nested defs so that objfn's template inherits from profil's template 
    
    plower <- lower[profiledOutPars]
    pupper <- upper[profiledOutPars]
    if (inherits(object,"SLik_j")) { # then guess a good init from a high instrumental-distr. density (assuming that logL was used to update this density)
      init <- .safe_init(object=object, given=template[c(st, fixedPars)], plower, pupper,
                         more_inits=init_from_neigh)
    } 
    # optr <- .safe_optim(init[profiledOutPars],objfn=objfn,
    #                   lower=plower,upper=pupper, LowUp=list(), verbose=FALSE, template=template, object=object)
    optr <- .ad_hoc_opt(init[profiledOutPars], objfn=objfn, lower=plower, upper=pupper, gradient=gradient, 
                        safe=safe, LowUp=list(), template=template, object=object, neg_ineq_constrfn=neg_ineq_constrfn)
    profpt <- template
    profpt[profiledOutPars] <- optr$solution
    value <- - optr$objective
    if (value>object$MSL$maxlogL+1e-6) {
      object$MSL$init_from_prof <- profpt
      object$MSL$maxlogL <- value
    }
    list(value=value,profpt=profpt)
  } 
  maxlogL <- object$MSL$maxlogL
  for (level in CIlevels) {
    LRvalue <- qchisq(level, df=1)/2
    if (is.null(ci_info <- .get_ci_info(object$CIobject$CIs,parm=st,level=level))) {
      ci_info <- confint(object,parm=st, level=level, verbose=FALSE)
      ci_info$prof_y <- maxlogL-LRvalue 
      .assign_ci_info(ci_info, object$CIobject, parm=st, level=level)
    }
  } 
  # prof_y handling is a bit clumsy. The constraints are that the block
  # if (length(CIlevels)) { CIs4parm <- object$CIobject$CIs[[st]] ...
  # expects prof_y in each sublist,
  # and we additionally need prof_y in the next ci_info
  ci_info <- .get_ci_info(object$CIobject$CIs,parm=st,level=0.95,
                          prof_y= maxlogL-qchisq(0.95, df=1)/2)

  if (is.null(prof <- oldprofiles[[st]])) {
    dx <- upper-lower
    margin <- dx/10000
    safelow <- lower+margin
    safeup <- upper-margin
    x <- seq(safelow[st],safeup[st],length.out=gridSteps)
    y <- numeric(gridSteps)
    profpts <- matrix(NA, nrow=gridSteps+1L, ncol=length(MSLE))
    colnames(profpts) <- names(MSLE)
    cum_it <- 0L
    ## First profile, initialization from the MSLE
    init_from_neigh <- object$MSL$MSLE
    x_is_below <- (x<init_from_neigh[st])
    # going left from MSLE (excluded)
    for (ptit in rev(which(x_is_below))) {
      prof_it <- profil(x[ptit])   ####### the optimization
      y[ptit] <- prof_it$value
      profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
      if (interactive()) {
        cum_it <- cum_it+1L
        msg <- paste("Already ", cum_it, " profile points computed out of ", gridSteps, "     ", sep="")
        prevmsglength <- .overcat(msg, prevmsglength)
      }
    }
    init_from_neigh <- object$MSL$MSLE # If I comment this I may perform a reverse pass...
    x_is_above <- ! x_is_below
    # going right from MSLE (included -- y has to be computed...)
    for (ptit in which(x_is_above)) {
      prof_it <- profil(x[ptit])   ####### the optimization
      y[ptit] <- prof_it$value
      profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
      if (interactive()) {
        cum_it <- cum_it+1L
        msg <- paste("Already ", cum_it, " profile points computed out of ", gridSteps, "     ", sep="")
        prevmsglength <- .overcat(msg, prevmsglength)
      }
    }
    ## Trying to improve the profile, by initialization from CI bounds
    # We break any of the loop as soon as a point is not improved within it
    # Since it means we don't have a better init for the next point.
    pb_with_ci_lowbnd <- pb_with_ci_upbnd <- FALSE
    if ( ! is.null(ci_info)) {
      if (interactive()) cat("...using CI info...")
      neigh_lo <- neigh_up <- NULL
      init_from_neigh <- ci_info$lowerpar # may be NULL
      if (is.numeric(init_from_neigh)) {
        # going left from CI lower bound
        x_is_below <- (x<init_from_neigh[st])
        for (ptit in rev(which(x_is_below))) {
          prof_it <- profil(x[ptit])   ####### the optimization
          if (prof_it$value> y[ptit]) {
            y[ptit] <- prof_it$value
            profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
          } else { 
            pb_with_ci_lowbnd <- (y[ptit]>ci_info$prof_y) # check whether preexisting profile value should actually be in CI...
            break
          }
        }
      }
      init_from_neigh <- ci_info$upperpar
      if (is.numeric(init_from_neigh)) {
        # going right from CI upper bound
        x_is_above <- (x>init_from_neigh[st])
        for (ptit in which(x_is_above)) {
          prof_it <- profil(x[ptit])   ####### the optimization
          if (prof_it$value> y[ptit]) {
            y[ptit] <- prof_it$value
            profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
          } else { 
            pb_with_ci_upbnd <- (y[ptit]>ci_info$prof_y) # check whether preexisting profile value should actually be in CI...
            break
          }
        }
      }
      which_in_ci1 <- which(! (x_is_below | x_is_above))
      which_in_ci2 <- which(y>ci_info$prof_y)
      which_in_ci <- c(which_in_ci1,which_in_ci2)
      if (length(which_in_ci)) {
        which_in_ci <- range(which_in_ci)
        which_in_ci <- seq(which_in_ci[1L],which_in_ci[2L])
        if ( (! length(which_in_ci1)) || which_in_ci[1L]<which_in_ci1[1L]) {
          init_from_neigh <- profpts[max(1L,which_in_ci[1L]-1L), ]
        } else init_from_neigh <- ci_info$lowerpar # may be NULL
        if ( is.numeric(init_from_neigh)) { 
          # going *right* (within CI) from CI lower bound
          for (it in  seq_along(which_in_ci)) {
            ptit <- which_in_ci[it]
            prof_it <- profil(x[ptit])   ####### the optimization
            if (prof_it$value> y[ptit]) {
              y[ptit] <- prof_it$value
              profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
            } else {
              # look whether eny of the next points has suspiciously low lik
              worth_trying_next_point <- (
                it < length(which_in_ci) &&
                  any(y[which_in_ci[(it+1L):length(which_in_ci)]]< ci_info$prof_y))
              if ( ! worth_trying_next_point) break
            }
          }
        } 
        if ( (! length(which_in_ci1)) || tail(which_in_ci,1L)> tail(which_in_ci1,1L)) {
          init_from_neigh <- profpts[min(length(x),tail(which_in_ci,1L)+1L), ]
        } else init_from_neigh <- ci_info$upperpar # may be NULL
        if ( is.numeric(init_from_neigh)) { 
          rev_which_in_ci <- rev(which_in_ci)
          for (it in  seq_along(rev_which_in_ci)) {
            ptit <- rev_which_in_ci[it]
            prof_it <- profil(x[ptit])   ####### the optimization
            if (prof_it$value> y[ptit]) {
              y[ptit] <- prof_it$value
              profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
            } else {
              # look whether eny of the next points has suspiciously low lik
              worth_trying_next_point <- (
                it < length(rev_which_in_ci) &&
                  any(y[rev_which_in_ci[(it+1L):length(rev_which_in_ci)]]< ci_info$prof_y))
              if ( ! worth_trying_next_point) break
            }
          }
        } 
      } 
    }
    # for (ptit in seq_len(length(x))) {
    #   prof_it <- profil(x[ptit])   ####### the optimization
    #   y[ptit] <- prof_it$value
    #   profpts[ptit, ] <- init_from_neigh <- prof_it$profpt
    #   if (interactive()) {
    #     msg <- paste("Already ", ptit, " profile points computed out of ", gridSteps, "     ", sep="")
    #     prevmsglength <- .overcat(msg, prevmsglength)
    #   }
    # }
    x[gridSteps+1L] <- MSLE[st]
    y[gridSteps+1L] <- maxlogL
    profpts[gridSteps+1L,] <- object$MSL$MSLE
    
    if (length(CIlevels)) {
      CIs4parm <- object$CIobject$CIs[[st]]
      interval_xy <- list()
      if (pb_with_ci_lowbnd) {
        lowers <- lapply(CIs4parm,function(li) {NA})
      } else lowers <- lapply(CIs4parm,function(li) {li$lowerpar})
      if (pb_with_ci_upbnd) {
        uppers <- lapply(CIs4parm,function(li) {NA})
      } else uppers <- lapply(CIs4parm,function(li) {li$upperpar})
      # the list elements are either numeric vectors or asingle NA... 
      ordre <- order(c(seq_len(length(lowers)),seq_len(length(uppers))))
      
      bounds <- c(lowers,uppers)[ordre]
      checkvec <- function(vec) {if (identical(vec,NA)) {return(NULL)} else {return(vec)} }
      whichNAs <- unlist(lapply(bounds,function(vec) {identical(vec,NA)}))
      bounds <- do.call(rbind,bounds[ ! whichNAs, drop=FALSE])
      interval_xy <- lapply(CIs4parm,function(li) {cbind(x=li$interval, y=li$prof_y)})
      interval_xy <- do.call(rbind,interval_xy)[! whichNAs,,drop=FALSE] 
      x <- c(x, interval_xy[,"x"])
      y <- c(y, interval_xy[,"y"])
      profpts <- rbind(profpts, bounds)
    }
    
    order_x <- order(x)
    x <- x[order_x]
    y <- y[order_x]
    profpts <- profpts[order_x,, drop=FALSE]
    
    y <- switch(type,
                "LR" = exp(y-maxlogL),
                "logL" = y,
                y-maxlogL
    )
    list(x=x,y=y, profpts=profpts)
  } else {
    prof # returns the old profile.
  }
}


.plot_itself <- function(profile_obj,
                         type, control,  
                         st, MSL_info, xy_add, 
                         xlab, ylab, .scale,
                         decorations,
                         ...) {
  x <- profile_obj$x
  y <- profile_obj$y
  .control <- list(min=-7.568353, shadow_col="grey70") # using -qchisq(0.9999, df=1)/2
  .control[names(control)] <- control
  if (type %in% c("zoom","dual","ranges")) {
    top_pts <- (y > .control$min)
    zoomable <- (length(unique(top_pts)) == 2L) # do not modify type locally as this would affect plot for next param
  }
  if (type=="dual" && zoomable) {
    shadow_col <- .control$shadow_col
    plot(x,y,xlab="",ylab="",axes=FALSE,frame.plot=FALSE,col=shadow_col)
    xrange <- range(x)
    yrange <- range(y)
    axis(side=3, at = pretty(xrange), col.ticks =shadow_col, col.axis=shadow_col) # col.axis is for the tick labels not the axis! 
    axis(side=4, at = pretty(yrange), col.ticks =shadow_col, col.axis=shadow_col)
    par(new = TRUE)
    plot(x[top_pts],y[top_pts],xlab=xlab,ylab=ylab,axes=FALSE,frame.plot=TRUE)
  } else if (type=="ranges") {
    x <- x[top_pts]
    y <- y[top_pts]
    xrange <- range(x)
    yrange <- range(c(y,.control$min))
    plot(x,y,xlab=xlab,ylab=ylab,axes=FALSE,
         frame.plot=TRUE, xlim=xrange, ylim=yrange)
  } else {
    if (type %in% c("zoom") && zoomable) {
      x <- x[top_pts]
      y <- y[top_pts]
    }
    xrange <- range(x)
    yrange <- range(y)
    plot(x,y,xlab=xlab,ylab=ylab,axes=FALSE,frame.plot=TRUE)
  }
  lines(x,y) ## use plotpar$lty to cancel the effect of this
  x_text <- x[1]
  if (interactive()) cat("\n")
  abline(h=(-qchisq(0.95, df=1)/2), col=2)
  text(x_text, (-qchisq(0.95, df=1)/2), "0.95", pos=3, col=2,offset=0.1) 
  abline(h=(-qchisq(0.99, df=1)/2), col=3)
  text(x_text, (-qchisq(0.99, df=1)/2), "0.99", pos=1, col=3,offset=0.1)
  if (type=="logLR") abline(h=.control$min)
  MSLy <- switch(type,
                 LR = 1,
                 logL = MSL_info$maxlogL,
                 0
  )
  points(MSL_info$MSLE[st],MSLy,pch="+")
  xscale <- .wrap_scale(.scale, values=x, valrange=xrange) 
  decorations(par=st) # if (!is.null(expectation)) abline(v=expectation[st],col="red")
  axis(1, at=xscale$at, labels=xscale$labels)
  axis(2)
  if ( ! is.null(xy_add)) {
    lines(xy_add$x,xy_add$y, ...)
  }
  profile_obj$xscale <- xscale
  profile_obj$xlab <- xlab
  profile_obj$ylab <- ylab
  
  profile_obj
}

.paral_calc_prof_coords <- function(cluster_args,
                                      pars, object, fixedPars, safe, MSLE, 
                                      template, lower, upper, oldprofiles, 
                                      gridSteps, type) {
  
  subobject <- object
  subobject$projectors <- NULL
  attr(subobject$logLs,"Simulate") <- NULL
  attr(subobject$logLs,"control.Simulate") <- NULL
  subobject$reftable_raw <- NULL
  
  cluster_args <- .set_cluster_type(cluster_args=cluster_args) # PSOCK vs FORK
  cores_info <- .init_cores(cluster_args=cluster_args)
  if (cluster_args$type=="FORK") {
    nb_cores <- cores_info$nb_cores
    #if (is.null(mc.silent <- control$mc.silent)) mc.silent <- TRUE 
    mc.silent <- TRUE # no 'control' argument in contrast to spaMM::dopar
    #if (is.null(mc.preschedule <- control$mc.preschedule)) mc.preschedule <- TRUE 
    mc.preschedule <- TRUE # no 'control' argument in contrast to spaMM::dopar    (___F I X M E___?)
    has_progressr <- ("package:progressr" %in% search())
    if (has_progressr) {
      # progressor is the only progress function that 'works' with mclapply
      # although not with load-balancing (mc.preschedule=FALSE)
      # Here we use the default (no balancing), and it is the steps with max value shown below that are reported.  
      prog_fn <- get("progressor", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
      with_fn <- get("with_progress", asNamespace("progressr"), inherits=FALSE) # syntax for using an undeclared package (cf stats:::confint.glm)
      with_fn({
        p <- prog_fn(steps=length(pars))
        p_calc_prof_coords <- function(st, ...) {
          res <- .calc_prof_coords(st, ...) 
          p() # p() call necessary for actual progress report 
          res
        }
        profiles <- try(
          parallel::mclapply(pars, FUN = p_calc_prof_coords, mc.silent=mc.silent, mc.cores=nb_cores,
                             mc.preschedule = mc.preschedule, 
                             object=subobject, fixedPars=fixedPars, safe=safe, MSLE=MSLE, 
                             template=template, lower=lower, upper=upper, oldprofiles=oldprofiles, 
                             gridSteps=gridSteps, type=type)
        )
      })
    } else {
      .warn_once_progressr()
      profiles <- try(
        parallel::mclapply(pars, FUN = .calc_prof_coords, mc.silent=mc.silent,  mc.cores=nb_cores,
                           mc.preschedule = mc.preschedule,
                           object=subobject, fixedPars=fixedPars, safe=safe, MSLE=MSLE, 
                           template=template, lower=lower, upper=upper, oldprofiles=oldprofiles, 
                           gridSteps=gridSteps, type=type)
      )
    }
  } else { # PSOCK
    cl <- cores_info$cl
    packages <- "Infusion"
    parallel::clusterExport(cl, "packages",envir=environment()) ## passes the list of packages to load
    abyss <- parallel::clusterEvalQ(cl, {sapply(packages,library,character.only=TRUE)}) ## snif
    if (cores_info$has_doSNOW) {
      show_pb <- (# verbose$most && 
        ! isTRUE(getOption('knitr.in.progress')))
      if (show_pb) {
        pb <- txtProgressBar(max = length(object), style = 3, char="P")
        progress <- function(n) setTxtProgressBar(pb, n)
        parallel::clusterExport(cl=cl, "progress",envir=environment()) ## slow?
        .options.snow = list(progress = progress)
      } else .options.snow = NULL
      st <- NULL ## otherwise R CMD check complains that no visible binding for global variable 'st'
      foreach_args <- list(
        st = pars, 
        .packages= packages,
        .options.snow = .options.snow,
        .inorder = TRUE, .errorhandling = "remove"
        #                                 "pass"## "pass" to see error messages
      )
      foreach_blob <- do.call(foreach::foreach,foreach_args)
      profiles <- foreach::`%dopar%`(foreach_blob,
                                     .calc_prof_coords(st, object=subobject, fixedPars=fixedPars, safe=safe, MSLE=MSLE, 
                                                       template=template, lower=lower, upper=upper, oldprofiles=oldprofiles, 
                                                       gridSteps=gridSteps, type=type) )
      if (show_pb) close(pb)
    } else { # PSOCK without doSNOW
      pbopt <- pboptions(nout=min(100,2*length(object)),type="timer", char="p")
      profiles <- pblapply(X=pars, FUN = .calc_prof_coords, cl= cl, 
                           object=subobject, fixedPars=fixedPars, safe=safe, MSLE=MSLE, 
                           template=template, lower=lower, upper=upper, oldprofiles=oldprofiles, 
                           gridSteps=gridSteps, type=type)
      pboptions(pbopt)
    }
  }
  .close_cores(cores_info)
  profiles
} 

.fittedPars_plot_pars <- function(np, plotpar, ylab, mai) {
  intsqrt <- floor(sqrt(np))
  if (intsqrt>1) {loccex.axis <- par("cex.axis")*0.6} else {loccex.axis <- par("cex.axis")}
  par_arglist <- list(mfrow=c(ceiling(np/intsqrt), intsqrt), cex.axis=loccex.axis)
  par_arglist$mai <- mai #avoids explicit NULL, not handled
  par_arglist <- c(par_arglist, plotpar)
  if (length(grep("\n", ylab))&& is.null(par_arglist$mgp)) par_arglist$mgp <- c(2,1,0)
  par(par_arglist, no.readonly = TRUE)
}

.reformat_bound <- function(bound, template) {
  if ( ! is.null(bound)) template[names(bound)] <- bound
  template
}

plot1Dprof <- function(object, ## SLik object 
                       pars=object$colTypes$fittedPars, ## for which parameters profiles will be plotted
                       fixed=NULL,
                       type="logLR", ## type of plot: see switch below for possible types  
                       gridSteps=21, ## number of points for plot; 0 produces a 'curve', but may be slow
                       ## 21 is compromise between niceness and computation time
                       xlabs=list(), ## x axis names; non-defaut is a list with names= some of the 'pars' 
                       ylab, ## y-axis name; default deduced from type
                       scales=NULL,
                       plotpar=list(pch=20), ## graphic parameters, a list of valid arguments for par() 
                       control=list(min=-7.568353, shadow_col="grey70"),
                       decorations = function(par) NULL, # function(par) abline(v=expectation[par],col="red")
                       profiles=NULL,
                       add=NULL,
                       safe=TRUE, #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
                       cluster_args=NULL,
                       do_plot=TRUE,
                       CIlevels=NULL,
                       lower=NULL,
                       upper=NULL,
                       verbose=TRUE,
                       ...
) {
  if ( ! is.function(decorations)) decorations <- function(par) NULL
  if (gridSteps <=1L) stop("gridSteps<=1 is obsolete.")
  fixedPars <- names(fixed)
  template <- MSLE <- object$MSL$MSLE
  if (! is.null(fixedPars)) template[fixedPars] <- fixed
  lower <- .reformat_bound(lower, template=object$lower)
  upper <- .reformat_bound(upper, template=object$upper)
  .xlabs <- as.list(pars)
  names(.xlabs) <- pars
  .xlabs[names(xlabs)] <- xlabs
  np <- length(pars)
  .scales <- rep("identity",np)
  names(.scales) <- pars
  .scales[names(scales)] <- scales
  if (missing(ylab)) { ## ylab=NULL has the default ylab=NULL meaning of plot()
    ylab <- switch(type,
                   LR = "profile likelihood ratio",
                   logL = "log profile likelihood",
                   logLR = ,
                   zoom = ,
                   ranges = ,
                   dual = "log profile\n likelihood ratio",
                   stop("Unknown 'type'")
    )
  } 
  if (np>3L) {
    mai <- c(0.45,0.5,0.1,0.1)
  } else mai <- NULL
  opar <- .fittedPars_plot_pars(np=np, plotpar=plotpar, ylab=ylab, mai=mai)
  ##
  oldprofiles <- profiles
  profiles <- vector("list", length(pars))
  names(profiles) <- pars
  
  is_parallelized <- min(length(pars),max(0L, cluster_args$spec))>1L
  if (is_parallelized) {
    if (verbose) {
      if (length(setdiff(names(MSLE), pars))) {
        cat(paste("Computing profiles for", paste(pars, collapse=", "), "\n"))
      } else cat(paste("Computing profiles for all fitted parameters\n"))
    }
    profiles <- .paral_calc_prof_coords(
      cluster_args=cluster_args,
      pars=pars, object=object, fixedPars=fixedPars, 
      safe=safe, #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
      MSLE=MSLE, 
      template=template, lower=lower, upper=upper, oldprofiles=oldprofiles, 
      gridSteps=gridSteps, type=type)
    names(profiles) <- pars
    
    if (do_plot) for (st in pars) {
      # actual plot and returns updated input
      profiles[[st]] <- .plot_itself(profile_obj=profiles[[st]],
                                  type=type, control=control, st=st, 
                                  MSL_info=object$MSL, .scale=.scales[st], xy_add=add[[st]], 
                                  xlab=.xlabs[[st]],
                                  ylab=ylab,
                                  decorations=decorations,
                                  ...)
    }
  } else for (st in pars) {
    profiles[[st]] <- profile_xy <- .calc_prof_coords(
      st=st, object=object, fixedPars=fixedPars, 
      safe=safe,  #  .safe_optim vs nlminb in .ad_hoc_opt; NOT controlling predict(., which)
      MSLE=MSLE, 
      template=template, lower=lower, upper=upper, oldprofiles=oldprofiles, 
      gridSteps=gridSteps, type=type, CIlevels=CIlevels)
    if (do_plot) {
      # actual plot and returns updated input
      profiles[[st]] <- .plot_itself(profile_obj=profile_xy,
                                     type=type, control=control, st=st, 
                                    MSL_info=object$MSL, .scale=.scales[st], xy_add=add[[st]], 
                                     xlab=.xlabs[[st]],
                                     ylab=ylab,
                                     decorations=decorations,
                                     ...)
    }
  }
  par(opar)
  # for (st in pars) {
  #   profiles[[st]]$xscale <- .plot_itself(type, control, y, x, .xlabs, st, ylab, maxlogL, MSLE, .scales, add, ...)
  # }
  MSL_updating <- ( ! is.null(object$MSL$init_from_prof))
  if (MSL_updating) .MSL_update(object)
  ### should be obsolete now as .MSL_update() nullifies init_from_prof:
  # # this updates the MSL environment of the input object...
  # # ... but leaves the MSL$init_from_prof in it, (there may be forgotten reasons for that...? )
  # # so if we redraw the plots the object will be updated again even if no new max has been found, unless 
  # # we remove it now:
  # object$MSL$init_from_prof <- NULL
  ###
  par(opar)
  resu <- list(MSL_updated=MSL_updating, profiles=profiles)
  if (MSL_updating) message("Maximum summary-likelihood was updated within input object\n   as a better value was found during profiling.")
  # class(resu) <- c("prof1D", class(resu))
  invisible(resu)
}

#plot1DProf(slikT)
#plot1DProf(slikT,pars="logPop",gridSteps=0,scales=c(logPop="log10",logSel="log10"),xlabs=c(logPop="Population size"))

.combine_strlist <- function (a, ...) { # function that combines each profile 
  # in a structured list, AND
  ## that handles their "info" attributes, returning the best of them, if any.
  ## We avoid repeated v <- do.call(rbind, list(a, ...)) as follows.
  ## There must be an .init in the foreach arguments, which gives an init 'a' value so that the
  ## combine fn is called once the first item has been computed (rather than the first two).
  dotlist <- list(...)
  v <- c(a,dotlist) # joins preexisting list and new value(-s, potentially. Does this occur?)
  #
  if (is.null(newinfo <- attr(dotlist[[1]], "info"))) { # if NO info on new value...
    newinfo <- attr(a,"info") #keep preexisting info, if any
  } else if ( ! is.null(oldinfo <- attr(a,"info"))) { # else compare to preexisting info, if any...
    if (oldinfo$maxlogL>newinfo$maxlogL) newinfo <- oldinfo # ... and keep the best. 
  }
  attr(v,"info") <- newinfo
  v
}


.calc_2Dprof_xyz <- function(object, parpair, par1, par2, lower, upper, 
                             margefrac, gridSteps, fixedPars, cl, template, MSLE) {
  lob <- lower[par1]
  upb <- upper[par1]
  marge <- margefrac * (upb - lob)
  xGrid <- seq(lob + marge, upb - marge, length.out = gridSteps)
  
  lob <- lower[par2]
  upb <- upper[par2]
  marge <- margefrac * (upb - lob)
  yGrid <- seq(lob + marge, upb - marge, length.out = gridSteps)
  prevmsglength <- 0L
  profiledOutPars <- setdiff(object$colTypes$fittedPars,c(parpair,fixedPars))
  
  if (length(profiledOutPars)) {
    ugly_ad_hoc_fix <- ( is_FORK <- ( ! inherits(cl[[1]],"SOCKnode")) &&
                           inherits(object$jointdens,"dMixmod") ) # ("___F I X M E___": this code not handling Mclust parallelisation?)
    constr_crits <- object$constr_crits
    profil <- function(z) {
      if (ugly_ad_hoc_fix) suppressPackageStartupMessages(requireNamespace("Rmixmod", quietly = TRUE))
      template[parpair] <- z
      ## nested defs so that objfn's template inherits from profil's template ;
      ## but oobnd optimization code check for the template argument so it has to be passed... 
      constr_tuning <- FALSE
      objfn <- function(v, template) {
        template[profiledOutPars] <- v 
        - predict(object,template, which="safe", constr_tuning=constr_tuning)[1]
      }
      if ( ! is.null(constr_crits <- object$constr_crits)) {
        neg_ineq_constrfn <- function(v, template) {
          template[profiledOutPars] <- v 
          as.numeric(eval(constr_crits, envir = as.list(template)))
        }
      } else neg_ineq_constrfn <- NULL
      if (inherits(object,"SLik_j")) {
        init <- .safe_init(object,given = template[c(parpair, fixedPars)], 
                           plower=lower[profiledOutPars],pupper=upper[profiledOutPars], more_inits=object$MSL$MSLE)
      } else init <- MSLE[profiledOutPars]
      optr <- .safe_optim(init, objfn=objfn,
                          lower=lower[profiledOutPars],upper=upper[profiledOutPars], LowUp=list(), 
                          verbose=FALSE, object=object, template=template, neg_ineq_constrfn=neg_ineq_constrfn)
      value <- - optr$objective
      if ( ! is.null(optr$solution)) {
        profpt <- template
        profpt[profiledOutPars] <- optr$solution
        resu <- list(value=value, profpt=profpt)
      } else resu <- list(value=value, profpt=rep(NA_real_, length(template)))
      if (value>object$MSL$maxlogL+1e-6) {
        info <- list(maxlogL=value, init_from_prof=resu$profpt)
        attr(resu,"info") <- info
      }
      return(resu)
    } 
  } else {
    profil <- function(z) {
      template[parpair] <- z
      value <- predict(object,template, which="safe", constr_tuning= Inf)[1]
      resu <- list(value=value)
      if (value>object$MSL$maxlogL+1e-6) {
        next_init <- template
        info <- list(maxlogL=value, init_from_prof=next_init)
        attr(resu,"info") <- info
      }
      return(resu)
    } 
  }
  #
  par_grid <- expand.grid(xGrid,yGrid)
  par_grid <- t(as.matrix(par_grid))

  profile_blob  <- combinepar(par_grid,     profil, fit_env=list(),  
                        cluster=cl, 
                        showpbar = TRUE, 
                        ################################### ..., # no longer ignored. 
                        control=list(.combine=.combine_strlist, .final=identity, .init=NULL))
  profile_blob <- do.call(rbind, profile_blob)
  
  logLvalues <- unlist(profile_blob[,"value"])
  profpts <- array(do.call(cbind,profile_blob[,"profpt"]), 
                   dim=c(length(template),gridSteps,gridSteps))
  dimnames(profpts) <- list(names(template),NULL,NULL)
  if (length(logLvalues) != ncol(par_grid)) {
    stop("Some profile computations appear to have failed.")
    # problem is that if a bug in the depths of profil()
    # combinepar() hides it (_____F I X M E____) (still the case? I wrote better diagnostic code for another issue)
  }
  if ( ! is.null(info <- attr(logLvalues,"info"))) object$MSL$init_from_prof <- info$init_from_prof
  dim(logLvalues) <- c(length(xGrid),length(yGrid))
  
  maxlogL <- object$MSL$maxlogL
  list(x=xGrid, y=yGrid, logLvalues = logLvalues, profpts=profpts)
}

plot2Dprof <- function(object, ## SLik object
                       pars=object$colTypes$fittedPars, ## for which parameters profiles will be plotted
                       fixed=NULL,
                       type="logLR", ## type of plot: see switch below for possible types  
                       gridSteps=17, ## number of grid points in each dimension
                       xylabs=list(), ## x and y axis names; non-defaut is a list with names= some of the 'pars' 
                       main, ## main plot name; default deduced from type
                       scales=NULL,
                       plotpar=list(pch=20), ## graphic parameters, a list of valid arguments for par() 
                       margefrac = 0,
                       decorations = function(par1,par2) NULL,
                       filled.contour.fn = "spaMM.filled.contour",
                       cluster_args=NULL,
                       min_logLR = qchisq(0.95,df=length(object$colTypes$fittedPars))/2 +3,
                       lower=NULL,
                       upper=NULL,
                       color.palette=NULL,
                       profiles=NULL,
                       ...
) {
  if (is.null(color.palette)) color.palette <- .Inf_palette(variant="turbo")
  if ( ! is.function(decorations)) decorations <- function(par1,par2) NULL
  fixedPars <- names(fixed)
  lower <- .reformat_bound(lower, template=object$lower)
  upper <- .reformat_bound(upper, template=object$upper)
  if ((np <- length(lower))<3L) stop(paste0("Only ",np," fitted parameters: plot2Dprof() not meaningful; try plot(<object>, filled=TRUE) instead?"))
  if (margefrac >= 0.5) message("'margefrac' argument too large")
  .xylabs <- as.list(pars)
  names(.xylabs) <- pars
  .xylabs[names(xylabs)] <- xylabs
  if (is.list(pars)) {
    parnames <- unique(unlist(pars))
    parpairs <- expand.grid(pars[[1]],pars[[2]], stringsAsFactors = FALSE)
    parpairs <- as.matrix(parpairs)
    parpairs <- parpairs[ (! parpairs[,1]==parpairs[,2]),, drop=FALSE]
  } else if (is.matrix(pars)) {
    parpairs <- pars[ (! pars[,1]==pars[,2]),, drop=FALSE]
    parnames <- unique(parpairs)
  } else {
    parnames <- setdiff(pars, fixedPars)
    parpairs <- t(combn(parnames,2))
  }
  if (length(wrongnames <- setdiff(parnames,names(lower)))) 
    stop(paste("Invalid parameter names:",paste(wrongnames,collapse=", ")))
  .scales <- rep("identity",length(parnames))
  names(.scales) <- parnames
  .scales[names(scales)] <- scales
  template <- MSLE <- object$MSL$MSLE
  if (! is.null(fixedPars)) template[fixedPars] <- fixed
  maxlogL <- object$MSL$maxlogL
  if (missing(main)) { ## ylab=NULL has the default ylab=NULL meaning of plot()
    main <- switch(type,
                   LR = "profile Likelihood ratio",
                   logLR = "log profile Likelihood ratio",
                   logL = "log profile Likelihood",
                   stop("Unknown 'type'")
    )
  } 
  opar <- par(plotpar, no.readonly = TRUE)

  # before the loop:
  #cl <- .setCluster(cluster_args=cluster_args, iseed=NULL) # spaMM::.setCluster()
  cluster_args <- .set_cluster_type(cluster_args) # PSOCK vs FORK
  cores_info <- .init_cores(cluster_args=cluster_args)
  cl = cores_info$cl
  if (length(cl)>1L) object <- ..shrink(object)
  # 
  oldprofiles <- profiles
  profiles <- vector("list", nrow(parpairs))
  profilenames <- paste0(parpairs[,1],",",parpairs[,2])
  names(profiles) <- profilenames
  for (pairit in seq_len(nrow(parpairs))) {
    parpair <- parpairs[pairit,]
    par1 <- parpair[1]
    par2 <- parpair[2]
    
    if (is.null(profile_xyz <- oldprofiles[[profilenames[pairit]]])) {
      if (interactive()) {cat(paste("Computing profile for (", profilenames[pairit], ")\n",sep=""))}
      profiles[[pairit]] <- profile_xyz <- .calc_2Dprof_xyz(
        object=object, parpair=parpair, par1=par1, par2=par2, lower=lower, upper=upper, 
        margefrac=margefrac, gridSteps=gridSteps, fixedPars=fixedPars, cl=cl, 
        template=template, MSLE=MSLE)
    } 
    
    xGrid <- profile_xyz$x
    yGrid <- profile_xyz$y
    logLvalues <- profile_xyz$logLvalues
    
    Zvalues <- switch(type,
                      LR = exp(logLvalues-maxlogL),
                      logLR = logLvalues-maxlogL,
                      logL = logLvalues
    )      
    if (inherits(filled.contour.fn,"character")) filled.contour.fn <- get(filled.contour.fn)
    if (type=="logLR") {
      min_logLR <- - abs(min_logLR)
      zrange <- range(pmax(Zvalues, min_logLR), finite = TRUE)
      # Before the next pmax() which would otherwise correct the -Inf values:
      Zvalues[is.infinite(Zvalues)] <- NA
      # Using NA instead of -Inf avoids edge effects in the plots.
      # ... and presence of NA's is handled only by ***recent spaMM*** plotting utils 
      Zvalues <- pmax(Zvalues, min_logLR-0.1)
      levels <- pretty(zrange, 20L)
      levels <- c(levels, min_logLR-0.0001)
      levels <- sort(levels)
    } else {
      zrange <- range(Zvalues, finite = TRUE, na.rm=TRUE)
      levels <- pretty(zrange, 20L)
    }
    if (zrange[1]<log(0.05)) {
      levels <- c(levels, log(0.05))
      levels <- sort(levels)
    }
    
    xscale <- .wrap_scale(.scales[par1],xGrid)
    yscale <- .wrap_scale(.scales[par2],yGrid)
    try(filled.contour.fn(x = xGrid, y = yGrid, z = Zvalues,xlab=.xylabs[[par1]],ylab=.xylabs[[par2]],
                      levels=levels, color.palette=color.palette,
                      plot.axes={
                        axis(1, at=xscale$at, labels=xscale$labels)
                        axis(2, at=yscale$at, labels=yscale$labels)
                        eval(decorations(par1=par1,par2=par2))
                      }, 
                      key.axes={
                        axis(4);axis(4,at=c(log(0.1465), log(0.05)),
                             labels=c("1D CI","2D CI"), cex.axis=0.5)}, 
                      main=main
    ))
    cat("\n")
    
  }
  if (length(cl)) parallel::stopCluster(cl)
  MSL_updating <- ( ! is.null(object$MSL$init_from_prof))
  if (MSL_updating) .MSL_update(object) 
  par(opar)
  resu <- list(MSL_updated=MSL_updating, profiles=profiles, 
               fit_info=object[c("lower","upper","MSL","colTypes")])
  if (MSL_updating) message("Maximum summary-likelihood was updated within input object\n   as a better value was found during profiling.")
  invisible(resu)
}

#plot2DProf(slikT,pars=c("logPop","logSel"),scales=c(logPop="log10",logSel="log10"),
#           xylabs=c(logPop="Population size",logSel="Selection coefficient"))
