viridisOpts <- function (n, alpha = 1, begin = 0, end = 1, option = "D",...) 
{
  if (begin < 0 | end < 0 | begin > 1 | end > 1) {
    stop("begin and end must be in [0,1]")
  }
  option <- switch(option, A = "A", magma = "A", B = "B", inferno = "B", 
                   C = "C", plasma = "C", D = "D", viridis = "D", {
                     warning(paste0("Option '", option, "' does not exist. Defaulting to 'viridis'."))
                     "D"
                   })
  map <- viridis::viridis.map[viridis::viridis.map$opt == option, ]
  map_cols <- grDevices::rgb(map$R, map$G, map$B)
  fn_cols <- grDevices::colorRamp(map_cols, space = "Lab", 
                                  interpolate = "spline",...)
  cols <- fn_cols(seq(begin, end, length.out = n))/255
  grDevices::rgb(cols[, 1], cols[, 2], cols[, 3], alpha = alpha)
}

plot.SLik <-function(x,y,filled=FALSE,
                     decorations=NULL, ## ADD to the default decorations
                     color.palette=NULL,
                     plot.axes=NULL, plot.title=NULL,
                     ...) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  maxlogL <- object$MSL$maxlogL
  if (is.null(maxlogL)) stop("plot.SLik plots likelihood ratios, hence it requires the ML to have been computed by MSL(.)\n")
  Ztransf <- function(Z) {exp(Z-maxlogL)}
  np <- length(fittedPars)
  if (np==1L) {
    xf <- object$obspred[,fittedPars]
    yf <- Ztransf(object$obspred[,attr(object$obspred,"fittedName")]) ## < 1
    x <- object$logLs[,fittedPars]
    y <- Ztransf(object$logLs[,object$colTypes$logLname]) ## <>1
    plot(xf,yf,main="Summary Likelihood Ratio",xlab=fittedPars,ylab="SL ratio",
         xlim=range(x),
         ylim=c(0,max(1,y)))
    points(x=x,y=y,pch=20,cex=0.5,col="red")
    points(object$MSL$MSLE,1,pch="+")
    if (!is.null(object$CIobject)) {
      yci <- rep(exp(-qchisq(0.95,1)/2),NROW(object$CIobject$bounds))
      points(object$CIobject$bounds,y=yci,pch="+")
    }
    #plot(object$fit,which="predict")
  } else if (np==2L) {
    if (inherits(object$fit,"HLfit")) {
      if (is.null(plot.title)) {
        ## also a language expression 
        plot.title <- substitute(title(main="Summary likelihood ratio surface",
                               xlab=xlabv,ylab=ylabv),list(xlabv=fittedPars[1],ylabv=fittedPars[2]))
      }
      if (filled) {
        decos <- quote({if (!is.null(object$latestPoints)) points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2);
          points(object$logLs[ ! object$logLs[,"isValid"],fittedPars],pch=20,cex=0.8);
          eval(decorations);
          if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch=19,col="red",cex=1.50);
          if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch=20,col="white",cex=0.75);
          points(t(object$MSL$MSLE),pch=19,cex=1.50,col="cyan");
          points(t(object$MSL$MSLE),pch=20,cex=0.75,col="white");
        }
        ) 
        if(is.null(plot.axes)) {
          # a language expression for filled.mapMM, as returned by quote() or substitute()
          plot.axes <- quote({axis(1);axis(2);
            # this is evaluatedat the end of spaMM::filled.mapMM
            contour(x=xGrid,y=yGrid,z=Zvalues,add=TRUE, nlevels=1, levels=c(0.05))})
        }
        if (is.null(color.palette)) {
          color.palette <- function(n){adjustcolor(viridisOpts(n,bias=2),offset = c(0.5, 0.5, 0.3, 0))}
        }
        decosf <- substitute({x;y},list(x=quote(points(object$logLs[,fittedPars],cex=0.8)),
                                        y=decos)) ## joins language objects
        filled.mapMM(object$fit,Ztransf=Ztransf,
                     color.palette=color.palette,nlevels=50,
                     plot.title=plot.title, ## language expr
                     plot.axes=plot.axes, ## language expr
                     decorations=eval(decosf), ## eval before an inner 'eval' is called,cf filled.mapMM code
                     ...)
      } else {
        decos <- quote({if (!is.null(object$latestPoints)) points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2);
          points(object$logLs[ ! object$logLs[,"isValid"],fittedPars],pch=20,cex=0.8);
          eval(decorations);
          if (!is.null(object$CIobject)) points(object$CIobject$bounds,pch="+",col="black",cex=1.50);
          points(t(object$MSL$MSLE),pch=3,cex=1.5,lwd=3); ## 3 rather than "+" for exact positioning
          points(t(object$MSL$MSLE),pch=3,col="white",cex=1.25);
        }
        ) 
        if(is.null(plot.axes)) {
          # a language expression for filled.mapMM, as returned by quote() or substitute()
          plot.axes <- quote({axis(1);axis(2)}) 
          # this is evaluated in spaMMplot2D where no grid of z values is available for contour
        }
        if (is.null(color.palette)) {
          color.palette <- function(n){spaMM.colors(n,redshift=3)}
        }
        mapMM(object$fit,Ztransf=Ztransf,
              color.palette=color.palette,nlevels=50,
              plot.axes=plot.axes, ## appears to accept eval(.) too 
              plot.title=plot.title, ## appears to accept eval(.) toobut best not (cf contours)
              decorations=decos,
              ...)
      }
    } 
  } else if (np>2L) {
    if (is.null(color.palette)) color.palette <- function(n){spaMM.colors(n,redshift=1/2)}
    intsqrt <- floor(sqrt(np))
    if (intsqrt>1) {loccex.axis <- par("cex.axis")*0.6} else {loccex.axis <- par("cex.axis")}
    ## mfrow marche pas avec rstudio (?) cf OKsmooth::provideDevice
    dev <- getOption("device")
    rstudioMess <-  ( (class(dev)=="character" && dev == "RStudioGD") )  
    if (! rstudioMess) opar <- par(cex.axis=loccex.axis, no.readonly = TRUE)
    #  # if (! rstudioMess) opar <- par(mfrow=c(ceiling(np/intsqrt), intsqrt), cex.axis=loccex.axis, no.readonly = TRUE)
    #  ## cf blackbox::gridfn, makeplot, etc
    # mais en fait migraine necase pas plusieurs filled plots sur un device; ici ca ferait planter car plot.new() -> figure margins too large
    margefrac <- 0.001
    gridsteps <- 40
    np <- length(fittedPars)
    grillelist <- list()
    ranges <- apply(object$fit$data[,fittedPars],2,range)
    for (it in seq_len(np-1)) {
      xvar <- fittedPars[it]
      range <- ranges[,xvar]
      lob <- range[1]
      upb <- range[2]
      marge <-  (upb-lob)/(4*gridsteps) ## just enough to see the maximum on the edge
      grillelist[xvar] <- list(seq(lob-marge, upb+marge, length.out=gridsteps)) ## uniform on log scale if relevant
      for (jt in (it+1):np) {
        yvar <- fittedPars[jt]
        fixedPars <- setdiff(fittedPars,c(xvar,yvar))
        grillelist[fixedPars] <- NULL
        fixedVals <- object$MSL$MSLE[fixedPars]
        range <- ranges[,yvar]
        lob <- range[1]
        upb <- range[2]
        marge <-  (upb-lob)/(4*gridsteps) 
        grillelist[yvar] <- list(seq(lob-marge, upb+marge, length.out=gridsteps)) ## uniform on log scale if relevant
        ## oddly, order in grillelist is not always well controlled at this point hence
        grillelist <- grillelist[c(xvar,yvar)] ## simply reorder elements according to this order
        grille <- expand.grid(grillelist) 
        grille <- cbind(grille,t(fixedVals))
        grille <- grille[,fittedPars] ## simply reorder grille elements according to fittedNames order
        z <- predict(object, grille)
        xyz <- as.surface(grillelist, z, order.variables = "xy")
        main <- paste("logL slice for",paste(fixedPars,signif(fixedVals,4),sep="=",collapse=", "))
        varVals <- object$MSL$MSLE[c(xvar,yvar)]
        if (interactive() && ! rstudioMess) plot.new() # provideDevice(bbdefaultPars=TRUE) does not work by itself, cf providePlotFile
        if(is.null(plot.axes)) {
          plot.axes <- quote(points(varVals[xyz$xlab],varVals[xyz$ylab],pch="+",cex=1.5)) 
        }
        spaMM.filled.contour(xyz$x, xyz$y, xyz$z,xlab=xyz$xlab,ylab=xyz$ylab,main=main,
                             color.palette= color.palette,
                             plot.axes=eval(plot.axes)
        )
      } 
    }
    if ( ! rstudioMess) par(opar)
  }
  invisible(object)
}

