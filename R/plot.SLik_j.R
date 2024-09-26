.log10_1p <- function(x) {log10(1+x)}

.wrap_scale <- function(fnname, values, valrange=range(values), invfnname, logticks, validRange) {
  xscale <- switch(
    fnname,
    "log" =  makeTicks(exp(values), axis=1,scalefn="log",logticks=TRUE,validRange=c(-Inf,Inf)),
    "log10"= makeTicks(10^values, axis=1,scalefn="log10",logticks=TRUE,validRange=c(-Inf,Inf)),
    "log1p"= makeTicks(expm1(values), axis=1,scalefn="log1p",logticks=TRUE,validRange=c(-Inf,Inf)),
    "log10_1p"= makeTicks(10^values-1, axis=1,scalefn=get(".log10_1p"),logticks=TRUE,validRange=c(-Inf,Inf)),
    "identity" = {
      at <- pretty(valrange)
      list(at=at, labels=at)
    },
    { 
      invfn <- get(invfnname)
      scalefn <- get(fnname)
      makeTicks(invfn(values), axis=1,scalefn=scalefn,logticks=logticks,validRange=validRange)
    }
  )
}

.filled.contour.plot <- function(
    object, xyvars, color.palette=.Inf_palette(variant="turbo"), 
    grid_args=list(), plot.title, 
    .scale=c("identity","identity") # I implemented '.scale' arg but default is always used.  
  ) {
  grillelist <- list()
  grid_args$values <- object$logLs[,xyvars[1L]]
  grillelist[[xyvars[1L]]] <- do.call(".gridfn",grid_args) 
  grid_args$values <- object$logLs[,xyvars[2L]]
  grillelist[[xyvars[2L]]] <- do.call(".gridfn",grid_args) 
  grille <- expand.grid(grillelist) 
  z <- predict(object, grille, which="safe", constr_tuning=Inf)
  z <- exp(z-object$MSL$maxlogL)
  xyz <- as.surface(grillelist, z, order.variables = "xy")
  varVals <- object$MSL$MSLE[xyvars]
  has_CI_info <- ( ( ! is.null(object$CIobject)) && is.null(object$CIobject$warn))
  xscale <- .wrap_scale(.scale[1], values=xyz$x) 
  yscale <- .wrap_scale(.scale[2], values=xyz$y) 

  plot.axes <- quote({
    if (!is.null(object$latestPoints)) points(object$logLs[object$latestPoints,xyvars],pch=".",cex=2);
    points(varVals[xyz$xlab],varVals[xyz$ylab],pch="+",cex=1.5) 
    if (has_CI_info) points(object$CIobject$bounds,pch=19,col="red",cex=1.50);
    if (has_CI_info) points(object$CIobject$bounds,pch=20,col="white",cex=0.75);
    points(t(object$MSL$MSLE),pch=19,cex=1.50,col="cyan");
    points(t(object$MSL$MSLE),pch=20,cex=0.75,col="white");
    contour(x=grillelist[[xyvars[1L]]],y=grillelist[[xyvars[2L]]],
            z=matrix(z, ncol = length(grillelist[[xyvars[2L]]])),
            add=TRUE, nlevels=1, levels=c(0.05));
    axis(1, at=xscale$at, labels=xscale$labels)
    axis(2, at=yscale$at, labels=yscale$labels)
  }
  ) 
  dev <- getOption("device")
  rstudioMess <-  ( (class(dev)=="character" && dev == "RStudioGD") )  
  if (is.null(color.palette)) color.palette <- .Inf_palette(variant="turbo")
  knitRmess <- isTRUE(getOption('knitr.in.progress')) # (class(dev)=="function" && environmentName(parent.env(environment(dev)))=="imports:knitr")
  if (interactive() && ! (rstudioMess || knitRmess)) plot.new() 
  spaMM.filled.contour(xyz$x, xyz$y, xyz$z,xlab=xyz$xlab,ylab=xyz$ylab,main="Summary-likelihood-ratio surface",
                       color.palette= color.palette,nlevels=50,
                       plot.axes=eval(plot.axes),
                       plot.title=plot.title
  ) # This now returns a list, invisibly.
}

plot.SLik_j <- function(x, # this is the SLik_j object; we will rename it and use x for other purposes, for clarity
                        y, filled = nrow(x$logLs)>5000L, decorations = NULL, 
                       color.palette = NULL, plot.axes = NULL, 
                       plot.title = NULL, 
                       from_refine=FALSE, # only controls whether MSL may be updated following computation of a grid of predictions.
                       plot.slices=TRUE, show_latest=FALSE,
                       # .scale=c("identity","identity"),
                       ... # may contain map.asp...
                       ) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  maxlogL <- object$MSL$maxlogL
  if (is.null(maxlogL)) stop("plot.SLik plots likelihood ratios, hence it requires the ML to have been computed by MSL(.)\n")
  Ztransf <- function(Z) {exp(Z-maxlogL)}
  np <- length(fittedPars)
  if (np==1L) {
    x <- object$logLs[,fittedPars,drop=FALSE]
    profy <- predict(object,newdata=x, which="safe", constr_tuning = Inf)
    maxprofy <- max(profy)
    if (maxprofy>maxlogL) {
      whichmax <- which.max(profy)
      object$MSL$init_from_prof <- x[whichmax, ]
      maxlogL <- maxprofy
    }
    y <- exp(profy-maxlogL) ## 
    ylim <- max(1,y)
    if (is.infinite(ylim)) {ylim <- NULL} else {ylim <- c(0,ylim)}
    plot(x[,1],y,main="Summary Likelihood Ratio",xlab=fittedPars,ylab="SL ratio",
         xlim=range(x),
         ylim=ylim)
    points(object$MSL$MSLE,1,pch="+")
    RESU <- list(x=x,y=y, MSLE=c(x=object$MSL$MSLE,y=1))
    has_CI_info <- ( ( ! is.null(object$CIobject)) && is.null(object$CIobject$warn))
    if (has_CI_info) {
      yci <- rep(exp(-qchisq(0.95,1)/2),NROW(object$CIobject$bounds))
      points(object$CIobject$bounds,y=yci,pch="+")
      RESU$CIbounds <- cbind(x=object$CIobject$bounds,y=yci)
    }
    invisible(RESU)
  } else if (np==2L) {
    if (is.null(plot.title)) plot.title=quote(title(main="Summary-likelihood-ratio surface",
                                              xlab=fittedPars[1],ylab=fittedPars[2]))
    if (filled) {
      RESU <- .filled.contour.plot(object, xyvars=fittedPars, color.palette=NULL, 
                                   # .scale=.scale,
                                   plot.title = eval(plot.title)) 
    } else {
      decos <- quote({
        if (show_latest && !is.null(object$latestPoints)) {points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2)}; # $latestPoints is added by refine()...
        if (!is.null(object$CIobject)) {points(object$CIobject$bounds,pch="+",cex=1.5,col="red")};
        points(t(object$MSL$MSLE),pch="+",cex=2);
        decorations}) ## language object
      if (is.null(color.palette)) color.palette <- .Inf_palette(variant="spaMM_shift3")
      xy <- object$logLs[,fittedPars]
      profz <- predict(object,newdata=xy, which="safe", constr_tuning = Inf)
      maxprofz <- max(profz)
      if (maxprofz>maxlogL) {
        whichmax <- which.max(profz)
        object$MSL$init_from_prof <- xy[whichmax, ]
        maxlogL <- maxprofz
      }
      z <- exp(profz-maxlogL)
      try(spaMMplot2D(x=xy[,1],y=xy[,2],z=z,
                  color.palette=color.palette ,nlevels=50,
                  plot.title=eval(plot.title),
                  key.axes={axis(4);axis(4,at=c(0.1465, 0.05),labels=c("1D CI","2D CI"), cex.axis=0.5)}, 
                  decorations=decos, 
                  ...) # spaMMplot2D may receive its map.asp arg through these dots ; default determined by spaMM:::.calc_plot_dim() 
      ) # try() to automatically handle Rstudio's 'Plot region too large' errors for user-level plot(<SLik...>) calls.  
      RESU <- list(x=xy[,1],y=xy[,2],z=z, CIbounds=object$CIobject$bounds, MSLE=object$MSL$MSLE)
    }
  } else if (plot.slices) {
    RESU <- .calc_all_slices(object,fittedPars,color.palette,plot.axes=plot.axes)
  } else RESU <- NULL # This one makes sense: no plot -> no plot coordinates to return. [refine(., plot.slices=FALSE) leads here]
  if ( (! is.null(object$MSL$init_from_prof)) && (! from_refine) ) .MSL_update(object)
  invisible(RESU)
}
