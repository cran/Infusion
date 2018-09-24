.filled.contour.plot <-function(object,xyvars,color.palette=NULL,grid_args=list()) {
  grillelist <- list()
  grid_args$values <- object$logLs[,xyvars[1L]]
  grillelist[[xyvars[1L]]] <- do.call(".gridfn",grid_args) 
  grid_args$values <- object$logLs[,xyvars[2L]]
  grillelist[[xyvars[2L]]] <- do.call(".gridfn",grid_args) 
  grille <- expand.grid(grillelist) 
  z <- predict(object, grille)
  z <- exp(z-object$MSL$maxlogL)
  xyz <- as.surface(grillelist, z, order.variables = "xy")
  varVals <- object$MSL$MSLE[xyvars]
  plot.axes <- quote({
    if (!is.null(object$latestPoints)) points(object$logLs[object$latestPoints,xyvars],pch=".",cex=2);
    points(varVals[xyz$xlab],varVals[xyz$ylab],pch="+",cex=1.5) 
    if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch=19,col="red",cex=1.50);
    if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch=20,col="white",cex=0.75);
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
  if (interactive() && ! rstudioMess) plot.new() 
  if (is.null(color.palette)) {
    #color.palette <- function(n){spaMM.colors(n,redshift=3)}
    color.palette <- function(n){adjustcolor(.viridisOpts(n,bias=2),offset = c(0.5, 0.5, 0.3, 0))}
  }
  spaMM.filled.contour(xyz$x, xyz$y, xyz$z,xlab=xyz$xlab,ylab=xyz$ylab,main="Summary-likelihood-ratio surface",
                       color.palette= color.palette,nlevels=50,
                       plot.axes=eval(plot.axes)
  )
}

plot.SLik_j <-function(x, y, filled = nrow(x$logLs)>5000L, decorations = NULL, 
                       color.palette = NULL, plot.axes = NULL, 
                       plot.title = NULL, ...) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  maxlogL <- object$MSL$maxlogL
  if (is.null(maxlogL)) stop("plot.SLik plots likelihood ratios, hence it requires the ML to have been computed by MSL(.)\n")
  Ztransf <- function(Z) {exp(Z-maxlogL)}
  np <- length(fittedPars)
  if (np==1L) {
    x <- object$logLs[,fittedPars,drop=FALSE]
    y <- Ztransf(predict(object,newdata=x)) ## <>1
    plot(x[,1],y,main="Summary Likelihood Ratio",xlab=fittedPars,ylab="SL ratio",
         xlim=range(x),
         ylim=c(0,max(1,y)))
    points(object$MSL$MSLE,1,pch="+")
    if (!is.null(object$CIobject)) {
      yci <- rep(exp(-qchisq(0.95,1)/2),NROW(object$CIobject$bounds))
      points(object$CIobject$bounds,y=yci,pch="+")
    }
  } else if (np==2L) {
    if (filled) {
      .filled.contour.plot(object,xyvars=fittedPars,color.palette=NULL)
    } else {
      decos <- quote({
        if (!is.null(object$latestPoints)) {points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2)};
        if (!is.null(object$CIobject)) {points(object$CIobject$bounds,pch="+",cex=1.5,col="red")};
        points(t(object$MSL$MSLE),pch="+",cex=2);
        decorations}) ## language object
      xy <- object$logLs[,fittedPars]
      spaMMplot2D(x=xy[,1],y=xy[,2],z=Ztransf(predict(object,newdata=xy)),
                  color.palette=function(n){spaMM.colors(50,redshift=3)},nlevels=50,
                  plot.title=title(main="Summary-likelihood-ratio surface",
                                         xlab=fittedPars[1],ylab=fittedPars[2]),
                  decorations=decos,
                  ...)
    }
  } else if (np>2L) {.calc_all_slices(object,fittedPars,color.palette,plot.axes=plot.axes)}
  invisible(object)
}
