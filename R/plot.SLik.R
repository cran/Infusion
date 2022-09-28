.viridisOpts <- function (n, alpha = 1, begin = 0, end = 1, option = "D",...) 
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

.gridfn <- function(values, gridsteps=40, margefrac=NULL) {
  if (is.null(margefrac)) margefrac <- 1/(2*(gridsteps-1))
  if(margefrac>=0.5) stop("'margefrac' argument too large")
  .range <- range(values)
  marge <- margefrac*diff(.range)
  Grid <- seq(.range[1L]+marge, .range[2L]-marge, length.out=gridsteps) ## uniform on log scale if relevant
  return(Grid)
}

plot.SLik <-function(x,y,filled=FALSE,
                     decorations=NULL, ## ADD to the default decorations
                     color.palette=NULL,
                     plot.axes=NULL, plot.title=NULL,
                     plot.slices=TRUE,
                     ...) {
  object <- x
  fittedPars <- object$colTypes$fittedPars
  maxlogL <- object$MSL$maxlogL
  if (is.null(maxlogL)) stop("plot.SLik plots likelihood ratios, hence it requires the ML to have been computed by MSL(.)\n")
  Ztransf <- function(Z) {exp(Z-maxlogL)}
  np <- length(fittedPars)
  has_CI_info <- ( ( ! is.null(object$CIobject)) && is.null(object$CIobject$warn))
  if (np==1L) {
    xf <- object$obspred[,fittedPars]
    yf <- Ztransf(object$obspred[,attr(object$obspred,"fittedName")]) ## < 1
    x <- object$logLs[,fittedPars]
    y <- Ztransf(object$logLs[,object$colTypes$logLname]) ## <>1
    ylim <- max(1,y)
    if (is.infinite(ylim)) {
      ylim <- NULL ## most likely an horrendous plot from horrendous input, but avoids an error()
    } else {ylim <- c(0,ylim)}
    plot(xf,yf,main="Summary Likelihood Ratio",xlab=fittedPars,ylab="SL ratio",
         xlim=range(x),
         ylim=ylim)
    points(x=x,y=y,pch=20,cex=0.5,col="red")
    points(object$MSL$MSLE,1,pch="+")
    if (has_CI_info) {
      yci <- rep(exp(-qchisq(0.95,1)/2),NROW(object$CIobject$bounds))
      points(object$CIobject$bounds,y=yci,pch="+")
    }
    #plot(object$fit,which="predict")
  } else if (np==2L) {
    if (inherits(object$fit,"HLfit")) { # note the alternative
      if (is.null(plot.title)) {
        ## also a language expression  (but, much later... note that in plot.SLik_j()  quote() without substitute() appears to work)
        plot.title <- substitute(title(main="Summary-likelihood-ratio surface",
                                       xlab=xlabv,ylab=ylabv),list(xlabv=fittedPars[1],ylabv=fittedPars[2]))
      }
      if (filled) {
        decos <- quote({if (!is.null(object$latestPoints)) points(object$logLs[object$latestPoints,fittedPars],pch=".",cex=2);
          points(object$logLs[ ! object$logLs[,"isValid"],fittedPars],pch=20,cex=0.8);
          eval(decorations);
          if (has_CI_info) points(object$CIobject$bounds,pch=19,col="red",cex=1.50);
          if (has_CI_info) points(object$CIobject$bounds,pch=20,col="white",cex=0.75);
          points(t(object$MSL$MSLE),pch=19,cex=1.50,col="cyan");
          points(t(object$MSL$MSLE),pch=20,cex=0.75,col="white");
        }
        ) 
        if(is.null(plot.axes)) {
          # a language expression for filled.mapMM, as returned by quote() or substitute()
          plot.axes <- quote({axis(1);axis(2);
            # this is evaluated at the end of spaMM::filled.mapMM
            contour(x=xGrid,y=yGrid,z=Zvalues,add=TRUE, nlevels=1, levels=c(0.05))})
        }
        if (is.null(color.palette)) {
          color.palette <- function(n){adjustcolor(.viridisOpts(n,bias=2),offset = c(0.5, 0.5, 0.3, 0))}
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
          if (has_CI_info) points(object$CIobject$bounds,pch="+",col="black",cex=1.50);
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
              plot.title=plot.title, ## appears to accept eval(.) too but best not (cf contours)
              decorations=decos,
              ...)
      }
    } 
  } else if (plot.slices) {.calc_all_slices(object,fittedPars,color.palette,plot.axes=plot.axes)}
  invisible(object)
}

