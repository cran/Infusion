.filled.contour.plot <- function(object, xyvars, color.palette=NULL, grid_args=list(), plot.title) {
  grillelist <- list()
  grid_args$values <- object$logLs[,xyvars[1L]]
  grillelist[[xyvars[1L]]] <- do.call(".gridfn",grid_args) 
  grid_args$values <- object$logLs[,xyvars[2L]]
  grillelist[[xyvars[2L]]] <- do.call(".gridfn",grid_args) 
  grille <- expand.grid(grillelist) 
  z <- predict(object, grille, which="safe")
  z <- exp(z-object$MSL$maxlogL)
  xyz <- as.surface(grillelist, z, order.variables = "xy")
  varVals <- object$MSL$MSLE[xyvars]
  has_CI_info <- ( ( ! is.null(object$CIobject)) && is.null(object$CIobject$warn))
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
    axis(1); axis(2); 
  }
  ) 
  dev <- getOption("device")
  rstudioMess <-  ( (class(dev)=="character" && dev == "RStudioGD") )  
  knitRmess <- isTRUE(getOption('knitr.in.progress')) # (class(dev)=="function" && environmentName(parent.env(environment(dev)))=="imports:knitr")
  if (interactive() && ! (rstudioMess || knitRmess)) plot.new() 
  if (is.null(color.palette)) {
    #color.palette <- function(n){spaMM.colors(n,redshift=3)}
    color.palette <- function(n){adjustcolor(.viridisOpts(n,bias=2),offset = c(0.5, 0.5, 0.3, 0))}
  }
  spaMM.filled.contour(xyz$x, xyz$y, xyz$z,xlab=xyz$xlab,ylab=xyz$ylab,main="Summary-likelihood-ratio surface",
                       color.palette= color.palette,nlevels=50,
                       plot.axes=eval(plot.axes),
                       plot.title=plot.title
  )
}

plot.SLik_j <- function(x, y, filled = nrow(x$logLs)>5000L, decorations = NULL, 
                       color.palette = NULL, plot.axes = NULL, 
                       plot.title = NULL, from_refine=FALSE, 
                       plot.slices=TRUE, ...) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  maxlogL <- object$MSL$maxlogL
  if (is.null(maxlogL)) stop("plot.SLik plots likelihood ratios, hence it requires the ML to have been computed by MSL(.)\n")
  Ztransf <- function(Z) {exp(Z-maxlogL)}
  np <- length(fittedPars)
  if (np==1L) {
    x <- object$logLs[,fittedPars,drop=FALSE]
    profy <- predict(object,newdata=x, which="safe")
    maxprofy <- max(profy)
    if (maxprofy>maxlogL) object$MSL$init_from_prof <- maxlogL <- maxprofy
    y <- exp(profy-maxlogL) ## 
    ylim <- max(1,y)
    if (is.infinite(ylim)) {ylim <- NULL} else {ylim <- c(0,ylim)}
    plot(x[,1],y,main="Summary Likelihood Ratio",xlab=fittedPars,ylab="SL ratio",
         xlim=range(x),
         ylim=ylim)
    points(object$MSL$MSLE,1,pch="+")
    has_CI_info <- ( ( ! is.null(object$CIobject)) && is.null(object$CIobject$warn))
    if (has_CI_info) {
      yci <- rep(exp(-qchisq(0.95,1)/2),NROW(object$CIobject$bounds))
      points(object$CIobject$bounds,y=yci,pch="+")
    }
  } else if (np==2L) {
    if (is.null(plot.title)) plot.title=quote(title(main="Summary-likelihood-ratio surface",
                                              xlab=fittedPars[1],ylab=fittedPars[2]))
    if (filled) {
      .filled.contour.plot(object, xyvars=fittedPars, color.palette=NULL, plot.title = eval(plot.title))
    } else {
      decos <- quote({
        if (!is.null(object$latestPoints)) {points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2)};
        if (!is.null(object$CIobject)) {points(object$CIobject$bounds,pch="+",cex=1.5,col="red")};
        points(t(object$MSL$MSLE),pch="+",cex=2);
        decorations}) ## language object
      xy <- object$logLs[,fittedPars]
      profz <- predict(object,newdata=xy, which="safe")
      maxprofz <- max(profz)
      if (maxprofz>maxlogL) object$MSL$init_from_prof <- maxlogL <- maxprofz
      spaMMplot2D(x=xy[,1],y=xy[,2],z=exp(profz-maxlogL),
                  color.palette=function(n){spaMM.colors(50,redshift=3)},nlevels=50,
                  plot.title=eval(plot.title),
                  decorations=decos,
                  ...)
    }
  } else if (plot.slices) {.calc_all_slices(object,fittedPars,color.palette,plot.axes=plot.axes)}
  if ( (! is.null(object$MSL$init_from_prof)) && (! from_refine) ) .MSL_update(object)
  invisible(object)
}
