
add_reftable <- function(...) { ## fn for ABC like simulation
  add_simulation(...,nRealizations=1L)
} 

add_simulation <- function(simulations=NULL, Simulate, par.grid=NULL,
                           nRealizations=NULL,
                           newsimuls=NULL,verbose=interactive()) {
  old_nRealizations <- Infusion.getOption("nRealizations")
  if (is.null(nRealizations)) {
    nRealizations <- old_nRealizations
  } else {
    Infusion.options(nRealizations=nRealizations)
  }
  if ( ! is.null(simulations)) {
    lowersList <- list(old=attr(simulations,"LOWER"))  
    uppersList <- list(old=attr(simulations,"UPPER"))  
  } else lowersList <- uppersList <- list()
  if ( ! is.null(par.grid)) {
    if ( ! inherits(Simulate,"character")) stop("'Simulate' must be a character string")
    if ( ! inherits(par.grid,"data.frame")) stop("'par.grid' argument is not a data.frame")
    prevmsglength <- 0L
    nsim <- nrow(par.grid)
    gridsimuls <- list()
    for (ii in seq(nsim)) {
      par <- par.grid[ii,,drop=FALSE]
      simuls <- replicate(n=nRealizations,do.call(Simulate,par))
      if (is.null(dim(simuls))) {
        colName <- names(simuls[1])
        dim(simuls) <- c(length(simuls),1)
        colnames(simuls) <- colName
      } else if (nrow(simuls)>1L) simuls <- t(simuls)
      if(inherits(simuls,"numeric")) simuls <- matrix(simuls) ## if scalar summ stat.
      #colnames(simuls) <- stats
      if (is.null(colnames(simuls))) stop("The 'Simulate' function must provide names for the statistics.")
      attr(simuls,"par") <- par
      gridsimuls <- c(gridsimuls,list(simuls)) 
      if (verbose) {
        if (nRealizations>1 || ! ii %% 1000 || ii==nsim ) {
          msg <- paste(ii,"simulations run out of", nsim," ")
          prevmsglength <- .overcat(msg, prevmsglength) ## does not work well in interactive call to knitr -> redirection of stderr to file; no overwriting
        } 
      }
    }
    if (nRealizations>1) {
      simulations <- c(simulations,gridsimuls) 
    } else {
      gridsimuls <- do.call(rbind,gridsimuls)
      gridsimuls <- cbind(par.grid,gridsimuls)
      simulations <- rbind(simulations,gridsimuls)
    }
    if (verbose && prevmsglength>0L) cat("\n")
    lowersList$pargrid <- attr(par.grid,"LOWER") ## not necessarily present
    if (is.null(lowersList$pargrid)) lowersList$pargrid <- apply(par.grid,2,min)
    uppersList$pargrid <- attr(par.grid,"UPPER")
    if (is.null(uppersList$pargrid)) uppersList$pargrid <- apply(par.grid,2,max)
  } 
  if ( ! is.null(newsimuls)) {
    simulations <- c(simulations,newsimuls)
    ## also for SLik_j: otherwise the parNames info cannot be certain...
    newsimuls_pars <- do.call(rbind,lapply(newsimuls,attr,which="par")) ## allows the following check 
    if (is.null(newsimuls_pars)) {
      stop("'par' attribute appears to be missing from each of the new simulations.")
    }
    lowersList$new <- apply(newsimuls_pars,2,min)
    uppersList$new <- apply(newsimuls_pars,2,max)
  }
  attr(simulations,"LOWER") <- do.call(pmin,lowersList)
  attr(simulations,"UPPER") <- do.call(pmax,uppersList)
  attr(simulations,"Simulate") <- Simulate
  if ( nRealizations>1 && ! inherits(simulations,"EDFlist")) class(simulations) <- c("EDFlist",class(simulations))
  Infusion.options(nRealizations=old_nRealizations)
  return(simulations)
}  





