.calc_all_slices <- function(object,fittedPars,color.palette,plot.axes=NULL) {
  if (is.null(color.palette)) {
    ##               function(v) {spaMM.colors(v, redshift = 1/2)}
    color.palette <- function(n){adjustcolor(.viridisOpts(n,bias=2),offset = c(0.5, 0.5, 0.3, 0))}
  }
  np <- length(fittedPars)
  intsqrt <- floor(sqrt(np))
  if (intsqrt>1) {loccex.axis <- par("cex.axis")*0.6} else {loccex.axis <- par("cex.axis")}
  ## mfrow marche pas avec rstudio (?) cf OKsmooth::provideDevice
  dev <- getOption("device")
  rstudioMess <-  ( (class(dev)=="character" && dev == "RStudioGD") )  
  knitRmess <- isTRUE(getOption('knitr.in.progress')) # (class(dev)=="function" && environmentName(parent.env(environment(dev)))=="imports:knitr")
  if (! rstudioMess) opar <- par(cex.axis=loccex.axis, no.readonly = TRUE)
  #  # if (! rstudioMess) opar <- par(mfrow=c(ceiling(np/intsqrt), intsqrt), cex.axis=loccex.axis, no.readonly = TRUE)
  #  ## cf blackbox::gridfn, makeplot, etc
  # mais en fait migraine necase pas plusieurs filled plots sur un device; ici ca ferait planter car plot.new() -> figure margins too large
  grillelist <- list()
  grid_args <- list(gridsteps=40)
  grid_args$margefrac <- 1/(4*grid_args$gridsteps) ## just enough to see the maximum on the edge
  for (it in seq_len(np-1)) {
    xvar <- fittedPars[it]
    grid_args$values <- object$logLs[,xvar]
    grillelist[[xvar]] <- do.call(".gridfn",grid_args) 
    for (jt in (it+1):np) {
      yvar <- fittedPars[jt]
      fixedPars <- setdiff(fittedPars,c(xvar,yvar))
      grillelist[fixedPars] <- NULL
      fixedVals <- object$MSL$MSLE[fixedPars]
      grid_args$values <- object$logLs[,yvar]
      grillelist[[yvar]] <- do.call(".gridfn",grid_args) 
      ## Order in grillelist is not always well controlled at this point hence
      grillelist <- grillelist[c(xvar,yvar)] ## simply reorder elements according to this order
      grille <- expand.grid(grillelist) 
      grille <- cbind(grille,t(fixedVals))
      grille <- grille[,fittedPars] ## simply reorder grille elements according to fittedNames order
      z <- predict(object, grille)
      xyz <- as.surface(grillelist, z, order.variables = "xy")
      main <- paste("logL slice for",paste(fixedPars,signif(fixedVals,4),sep="=",collapse=", "))
      varVals <- object$MSL$MSLE[c(xvar,yvar)]
      if (interactive() && ! (rstudioMess || knitRmess)) plot.new() 
      if(is.null(plot.axes)) {
        # axis(1); axis(2); ## ? not in plot.SLik()  
        plot.axes <- quote({
          axis(1);axis(2);
          points(varVals[xyz$xlab],varVals[xyz$ylab],pch="+",cex=1.5) # locate the maximum
        }) 
      }
      spaMM.filled.contour(xyz$x, xyz$y, xyz$z,xlab=xyz$xlab,ylab=xyz$ylab,main=main,
                           color.palette=color.palette,
                           nlevels=50,
                           plot.axes=eval(plot.axes)
      )
    } 
  }
  if ( ! rstudioMess) par(opar)
}