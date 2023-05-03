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
  # so the object$MSL$init_from_prof remains
} # only environmentswhere updated, so no return value


plot1Dprof <- function(object, ## SLik object
                       pars=object$colTypes$fittedPars, ## for which parameters profiles will be plotted
                       type="logLR", ## type of plot: see switch below for possible types  
                       gridSteps=21, ## number of points for plot; 0 produces a 'curve', but may be slow
                       ## 21 is compromise between niceness and computation time
                       xlabs=list(), ## x axis names; non-defaut is a list with names= some of the 'pars' 
                       ylab, ## y-axis name; default deduced from type
                       scales=NULL,
                       plotpar=list(pch=20), ## graphic parameters, a list of valid arguments for par() 
                       control=list(min=-7.568353, shadow_col="grey70"),
                       decorations = function(par) NULL, # function(par) abline(v=expectation[par],col="red")
                       ...
) {
  lower <- object$lower
  upper <- object$upper
  template <- MSLE <- object$MSL$MSLE
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
                   logLR = "log profile likelihood ratio",
                   zoom = "log profile likelihood ratio",
                   dual = "log profile likelihood ratio",
                   logL = "log profile likelihood",
                   stop("Unknown 'type'")
    )
  } 
  if (FALSE) { # messy and buggy (par(mfrow) not being restored, when - out of Rstudio - it is here modified) 
    # opar <- par(plotpar, no.readonly = TRUE)
    # ##
    # dev <- getOption("device")
    # rstudioMess <-  ( (class(dev)=="character" && dev == "RStudioGD") )  
    # # mfrow marche pas avec rstudio (?) cf OKsmooth::provideDevice
    # if (! rstudioMess) {
    #   intsqrt <- floor(sqrt(np))
    #   if (intsqrt>1) {loccex.axis <- par("cex.axis")*0.6} else {loccex.axis <- par("cex.axis")}
    #   par(mfrow=c(ceiling(np/intsqrt), intsqrt),cex.axis=loccex.axis)
    # }
  } else {
    intsqrt <- floor(sqrt(np))
    if (intsqrt>1) {loccex.axis <- par("cex.axis")*0.6} else {loccex.axis <- par("cex.axis")}
    par_arglist <- c(list(mfrow=c(ceiling(np/intsqrt), intsqrt),cex.axis=loccex.axis), plotpar)
    opar <- par(par_arglist, no.readonly = TRUE)
  }
  ##
  for (st in pars) {
    if (interactive()) {cat(paste("Computing profile for", st, "\n"))}
    prevmsglength <- 0L
    profiledOutPars <- setdiff(object$colTypes$fittedPars,st)
    objfn <- function(v, template) {
      template[profiledOutPars] <- v 
      - predict(object,template, which="safe")[1]
    }
    profil <- function(z, init=MSLE) {
      template[st] <- z
      ## nested defs so that objfn's template inherits from profil's template 
      
      plower <- lower[profiledOutPars]
      pupper <- upper[profiledOutPars]
      if (inherits(object,"SLik_j")) { # then guess a good init from a high instrumental-distr. density (assuming that logL was used to update this density)
        init <- .safe_init(object=object, given=template[st], plower, pupper)
      } 
      optr <- .safe_opt(init[profiledOutPars],objfn=objfn,
                        lower=plower,upper=pupper, LowUp=list(), verbose=FALSE, template=template)
      value <- - optr$objective
      if (value>object$MSL$maxlogL+1e-6) {
        next_init <- template
        next_init[profiledOutPars] <- optr$solution
        object$MSL$init_from_prof <- next_init
        object$MSL$maxlogL <- value
      }
      return(value)
    } 
    if (gridSteps>0) {
      x <- seq(object$lower[st],object$upper[st],length.out=gridSteps)
      y <- numeric(gridSteps)
      for (ptit in seq_len(length(x))) {
        y[ptit] <- profil(x[ptit])
        if (interactive()) {
          msg <- paste("Already ", ptit, " profile points computed out of ", gridSteps, "     ", sep="")
          prevmsglength <- .overcat(msg, prevmsglength)
        }
      }
      maxlogL <- object$MSL$maxlogL
      y <- switch(type,
                     LR = exp(y-maxlogL),
                     logLR = y-maxlogL,
                     zoom = y-maxlogL,
                     dual = y-maxlogL,
                     logL = y
      )
      if (type %in% c("zoom","dual")) {
        .control <- list(min=-7.568353, shadow_col="grey70") # using -qchisq(0.9999, df=1)/2
        .control[names(control)] <- control
        top_pts <- (y > .control$min)
        zoomable <- (length(unique(top_pts)) == 2L) # do not modify type locally as this would affect plot for next param
      }
      if (type=="dual" && zoomable) {
        shadow_col <- .control$shadow_col
        plot(x,y,xlab="",ylab="",axes=FALSE,frame=FALSE,col=shadow_col)
        axis(side=3, at = pretty(range(x)), col.ticks =shadow_col, col.axis=shadow_col) # col.axis is for the tick labels not the axis! 
        axis(side=4, at = pretty(range(y)), col.ticks =shadow_col, col.axis=shadow_col)
        par(new = TRUE)
        plot(x[top_pts],y[top_pts],xlab=.xlabs[[st]],ylab=ylab,axes=FALSE,frame=TRUE)
      } else {
        if (type=="zoom" && zoomable) {
          x <- x[top_pts]
          y <- y[top_pts]
        }
        plot(x,y,xlab=.xlabs[[st]],ylab=ylab,axes=FALSE,frame=TRUE)
      }
      lines(x,y) ## use plotpar$lty to cancel the effect of this
      x_text <- x[1]
      if (interactive()) cat("\n")
    } else {
      stop("gridSteps<0 is obsolete.")
      # profil <- Vectorize(profil)
      # lowx <- object$lower[st]
      # hix <- object$upper[st]  
      # curve(profil(x),from=lowx,to=hix,xlab=.xlabs[st],ylab=ylab,axes=FALSE,frame=TRUE)
      # x_text <- lowx+ (hix-lowx)/40
    }
    abline(h=(-qchisq(0.95, df=1)/2), col=2)
    text(x_text, (-qchisq(0.95, df=1)/2), "0.95", pos=3, col=2,offset=0.1) 
    abline(h=(-qchisq(0.99, df=1)/2), col=3)
    text(x_text, (-qchisq(0.99, df=1)/2), "0.99", pos=1, col=3,offset=0.1)
    MSLy <- switch(type,
                   LR = 1,
                   logLR = 0,
                   zoom = 0,
                   dual = 0,
                   logL = maxlogL
    )
    points(MSLE[st],MSLy,pch="+")
    decorations(st) # if (!is.null(expectation)) abline(v=expectation[st],col="red")
    xscale <- switch(.scales[st],
                     "log" =  makeTicks(exp(x), axis=1,scalefn="log",logticks=TRUE,validRange=c(-Inf,Inf)),
                     "log10"= makeTicks(10^x, axis=1,scalefn="log10",logticks=TRUE,validRange=c(-Inf,Inf)),
                     "identity" = list(),
                     stop(paste("unknown scale '",.scales[st],"' for parameter ",st,sep=""))
    )
    axis(1, at=xscale$at, labels=xscale$labels)
    axis(2)
  }
  MSL_updating <- ( ! is.null(object$MSL$init_from_prof))
  if (MSL_updating) .MSL_update(object) # this updates the MSL environment of the input object...
  # ... but leaves the MSL$init_from_prof in it, (there may be forgotten reasons for that...? )
  # so if we redraw the plots the object will be updated again even if no new max has been found, unless 
  # we remove it now:
  object$MSL$init_from_prof <- NULL
  par(opar)
  return(list(MSL_updated=MSL_updating))
}

#plot1DProf(slikT)
#plot1DProf(slikT,pars="logPop",gridSteps=0,scales=c(logPop="log10",logSel="log10"),xlabs=c(logPop="Population size"))

plot2Dprof <- function(object, ## SLik object
                       pars=object$colTypes$fittedPars, ## for which parameters profiles will be plotted
                       type="logLR", ## type of plot: see switch below for possible types  
                       gridSteps=17, ## number of grid points in each dimension
                       xylabs=list(), ## x and y axis names; non-defaut is a list with names= some of the 'pars' 
                       main, ## main plot name; default deduced from type
                       scales=NULL,
                       plotpar=list(pch=20), ## graphic paratemers, a list of valid arguments for par() 
                       margefrac = 0,
                       decorations = function(par1,par2) NULL,
                       filled.contour.fn = "spaMM.filled.contour",
                       cluster_args=NULL,
                       ...
) {
  lower <- object$lower
  if ((np <- length(lower))<3L) stop(paste0("Only ",np," fitted parameters: plot2Dprof() not meaningful; try plot(<object>, filled=TRUE) instead?"))
  upper <- object$upper
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
    parnames <- pars
    parpairs <- t(combn(pars,2))
  }
  if (length(wrongnames <- setdiff(parnames,names(lower)))) 
    stop(paste("Invalid parameter names:",paste(wrongnames,collapse=", ")))
  .scales <- rep("identity",length(parnames))
  names(.scales) <- parnames
  .scales[names(scales)] <- scales
  template <- MSLE <- object$MSL$MSLE
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
  cl <- .setCluster(cluster_args=cluster_args, iseed=NULL) # spaMM::.setCluster()
  # 
  for (pairit in seq_len(nrow(parpairs))) {
    parpair <- parpairs[pairit,]
    par1 <- parpair[1]
    par2 <- parpair[2]
    
    lob <- lower[par1]
    upb <- upper[par1]
    marge <- margefrac * (upb - lob)
    xGrid <- seq(lob + marge, upb - marge, length.out = gridSteps)
    
    lob <- lower[par2]
    upb <- upper[par2]
    marge <- margefrac * (upb - lob)
    yGrid <- seq(lob + marge, upb - marge, length.out = gridSteps)

    if (interactive()) {cat(paste("Computing profile for (", paste(parpair,collapse=","), ")\n",sep=""))}
    prevmsglength <- 0L
    profiledOutPars <- setdiff(object$colTypes$fittedPars,parpair)
    if (TRUE) {
      ugly_ad_hoc_fix <- ( is_FORK <- ( ! inherits(cl[[1]],"SOCKnode")) &&
                             inherits(object$jointdens,"dMixmod") )
      profil <- function(z) {
        if (ugly_ad_hoc_fix) suppressPackageStartupMessages(requireNamespace("Rmixmod", quietly = TRUE))
        template[parpair] <- z
        ## nested defs so that objfn's template inherits from profil's template 
        objfn <- function(v) {
          template[profiledOutPars] <- v 
          - predict(object,template, which="safe")[1]
        }
        optr <- .safe_opt(MSLE[profiledOutPars],objfn=objfn,
                          lower=lower[profiledOutPars],upper=upper[profiledOutPars], LowUp=list(), verbose=FALSE)
        value <- - optr$objective
        if (value>object$MSL$maxlogL+1e-6) {
          next_init <- template
          next_init[profiledOutPars] <- optr$solution
          info <- list(maxlogL=value, init_from_prof=next_init)
          attr(value,"info") <- info
        }
        return(value)
      } 
      par_grid <- expand.grid(xGrid,yGrid)
      par_grid <- t(as.matrix(par_grid))
      Zvalues <- combinepar(par_grid,     profil, fit_env=list(),  
                            cluster=cl, 
                            showpbar = TRUE, 
                            ################################### ..., # no longer ignored. 
                            control=list(
                              .combine=function (a, ...) { # function that combines each profile in a vector, AND
                                                           # that handles their "info" attributes, returning the best of them, if any.
                                v <- c(a, ...) # oins preexisting vector and new value(-s, potentially. Does this occur?)
                                #
                                dotlist <- list(...)
                                if (is.null(newinfo <- attr(dotlist[[1]], "info"))) { # if NO info on new value...
                                  newinfo <- attr(a,"info") #keep preexisting info, if any
                                } else if ( ! is.null(oldinfo <- attr(a,"info"))) { # else compare to preexisting info, if any...
                                  if (oldinfo$maxlogL>newinfo$maxlogL) newinfo <- oldinfo # ... and keep the best. 
                                }
                                attr(v,"info") <- newinfo
                                v
                              }))
      if ( ! is.null(info <- attr(Zvalues,"info"))) object$MSL$init_from_prof <- info$init_from_prof
      dim(Zvalues) <- c(length(xGrid),length(yGrid))
    }  else {
      if ( ! is.null(cluster_args) ) {
        cluster_args <- .set_cluster_type(cluster_args=cluster_args) # PSOCK vs FORK
        if (cluster_args$type=="FORK") {
          cl <- parallel::makeForkCluster(nnodes = cluster_args$spec) 
        } else cl <- do.call(parallel::makeCluster, cluster_args) # note that _this_ line would make sense for fork clusters too. BUT
      } else cl <- NULL

      profil <- function(z) {
        template[parpair] <- z
        ## nested defs so that objfn's template inherits from profil's template 
        objfn <- function(v) {
          template[profiledOutPars] <- v 
          - predict(object,template, which="safe")[1]
        }
        optr <- .safe_opt(MSLE[profiledOutPars],objfn=objfn,
                          lower=lower[profiledOutPars],upper=upper[profiledOutPars], LowUp=list(), verbose=FALSE)
        value <- - optr$objective
        if (value>object$MSL$maxlogL+1e-6) {
          next_init <- template
          next_init[profiledOutPars] <- optr$solution
          object$MSL$init_from_prof <- next_init
          object$MSL$maxlogL <- value
        }
        return(value)
      } 
      Zvalues <- matrix(ncol=gridSteps,nrow=gridSteps)
      npts <- gridSteps^2
      ij <- 1L
      for (ptit in seq_len(gridSteps)) {
        for (ptjt in seq_len(gridSteps)) {
          Zvalues[ptjt,ptit] <- profil(c(xGrid[ptjt],yGrid[ptit]))
          if (interactive()) {
            msg <- paste("Already ", ij, " profile points computed out of ", npts, "     ", sep="")
            prevmsglength <- .overcat(msg, prevmsglength)
            ij <- ij + 1L
          }
        }
      }
    }
    maxlogL <- object$MSL$maxlogL
    Zvalues <- switch(type,
                      LR = exp(Zvalues-maxlogL),
                      logLR = Zvalues-maxlogL,
                      logL = Zvalues
    )      
    scalefn <- function(name,values) {
      switch(.scales[name],
             "log" =  makeTicks(exp(values), axis=1,scalefn="log",logticks=TRUE,validRange=c(-Inf,Inf)),
             "log10"= makeTicks(10^values, axis=1,scalefn="log10",logticks=TRUE,validRange=c(-Inf,Inf)),
             "identity" = list(),
             stop(paste("unknown scale '",.scales[name],"' for parameter ",name,sep=""))
      )
    }
    xscale <- scalefn(par1,xGrid)
    yscale <- scalefn(par2,yGrid)
    if (inherits(filled.contour.fn,"character")) filled.contour.fn <- get(filled.contour.fn)
    filled.contour.fn(x = xGrid, y = yGrid, z = Zvalues,xlab=.xylabs[[par1]],ylab=.xylabs[[par2]],
                      plot.axes={
                        axis(1, at=xscale$at, labels=xscale$labels)
                        axis(2, at=yscale$at, labels=yscale$labels)
                        eval(decorations(par1=par1,par2=par2))
                      },main=main
    )
    cat("\n")
    
  }
  if (length(cl)) parallel::stopCluster(cl)
  MSL_updating <- ( ! is.null(object$MSL$init_from_prof))
  if (MSL_updating) .MSL_update(object) # this updates the MSL environment of the input object...
  # ... but leaves the MSL$init_from_prof in it, (there may be forgotten reasons for that...? )
  # so if we redraw the plots the object will be updated again even if no new max has been found, unless 
  # we remove it now:
  object$MSL$init_from_prof <- NULL
  par(opar)
  return(list(MSL_updated=MSL_updating))
}

#plot2DProf(slikT,pars=c("logPop","logSel"),scales=c(logPop="log10",logSel="log10"),
#           xylabs=c(logPop="Population size",logSel="Selection coefficient"))
