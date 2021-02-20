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
  newmax <- MSL(object, eval_RMSEs=FALSE, CIs=FALSE)
  for (st in ls(newmax$MSL)) object$MSL[[st]] <- newmax$MSL[[st]] ## but keeping the environment (~pointer) unchanged.
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
                       control=list(min=-7.568353, shadow_col="grey70")
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
    opar <- par(mfrow=c(ceiling(np/intsqrt), intsqrt),cex.axis=loccex.axis, plotpar, no.readonly = TRUE)
  }
  ##
  for (st in pars) {
    if (interactive()) {cat(paste("Computing profile for", st, "\n"))}
    prevmsglength <- 0L
    profiledOutPars <- setdiff(object$colTypes$fittedPars,st)
    profil <- function(z) {
      template[st] <- z
      ## nested defs so that objfn's template inherits from profil's template 
      objfn <- function(v) {
        template[profiledOutPars] <- v 
        - predict(object,template, which="safe")[1]
        ## comparer a OKsmooth:::gridfn
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
  if (MSL_updating) .MSL_update(object)
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
                       margefrac = 0
) {
  lower <- object$lower
  upper <- object$upper
  if (margefrac >= 0.5) message("'margefrac' argument too large")
  .xylabs <- as.list(pars)
  names(.xylabs) <- pars
  .xylabs[names(xylabs)] <- xylabs
  np <- length(pars)
  .scales <- rep("identity",np)
  names(.scales) <- pars
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
  ##
  np <- length(pars)
  for (it in seq_len(np-1)) {
    par1 <- pars[it]
    lob <- lower[par1]
    upb <- upper[par1]
    marge <- margefrac * (upb - lob)
    xGrid <- seq(lob + marge, upb - marge, length.out = gridSteps)
    for (jt in (it+1):np) {
      par2 <- pars[jt]
      lob <- lower[par2]
      upb <- upper[par2]
      marge <- margefrac * (upb - lob)
      yGrid <- seq(lob + marge, upb - marge, length.out = gridSteps)
      parpair <- c(par1,par2)
      if (interactive()) {cat(paste("Computing profile for (", paste(parpair,collapse=","), ")\n",sep=""))}
      prevmsglength <- 0L
      profiledOutPars <- setdiff(object$colTypes$fittedPars,parpair)
      profil <- function(z) {
        template[parpair] <- z
        ## nested defs so that objfn's template inherits from profil's template 
        objfn <- function(v) {
          template[profiledOutPars] <- v 
          - predict(object,template, which="safe")[1]
          ## comparer a OKsmooth:::gridfn
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
      spaMM.filled.contour(x = xGrid, y = yGrid, z = Zvalues,xlab=.xylabs[[par1]],ylab=.xylabs[[par2]],
                                  plot.axes={
                                    axis(1, at=xscale$at, labels=xscale$labels)
                                    axis(2, at=yscale$at, labels=yscale$labels)
                                  },main=main
      )
      cat("\n")
    }
  }
  MSL_updating <- ( ! is.null(object$MSL$init_from_prof))
  if (MSL_updating) .MSL_update(object)
  par(opar)
  return(list(MSL_updated=MSL_updating))
}

#plot2DProf(slikT,pars=c("logPop","logSel"),scales=c(logPop="log10",logSel="log10"),
#           xylabs=c(logPop="Population size",logSel="Selection coefficient"))
