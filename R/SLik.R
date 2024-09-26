predict.SLik <- function (object, newdata=object$logLs[,object$colTypes$fittedPars,drop=FALSE],...) {
  return(predict(object$fit,newdata=newdata,...))
}

refine.SLik <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "REML"
  refine.default(object,surfaceData=object$logLs,method=method,...)
}

# finally not used for its initial purpose:
.get_nbCluster_from_SLik <- function(object) {
  if (inherits(object$jointdens,"dMixmod")) {
    nbCluster <- object$jointdens@nbCluster # for optional THIRD sampling step
  } else if (inherits(object$jointdens,"dMclust")) {
    nbCluster <- object$jointdens$G
  } else if (inherits(object$gllimobj,"gllim")) {
    nbCluster <- length(object$gllimobj$pi)
  } else nbCluster <- NULL
  nbCluster
}

# options("digits") controls digits in print()
## 3 digits on RMSEs makes a lot of digits overall   
summary.SLik <- function(object, ...) { 
  if ( !is.null(object$MSL) ) {
    nbClu <- .get_nbCluster_from_SLik(object)
    if ( ! is.null(nbClu)) {
      nbclustr <- paste(", joint density modeling:",nbClu,"clusters")
    } else nbclustr <- NULL
    cat(paste("*** Summary ML (",max(object$logLs$cumul_iter)," iterations, ",nrow(object$logLs)," points",
              nbclustr,"): ***\n",sep=""))
    print(c(object$MSL$MSLE,"logL"=object$MSL$maxlogL,"RMSE_logL"=unname(get_from(object,"RMSEs")[1L])))

    # # par_RMSEs
    if ( length(object$CIobject$CIs))  {
      if( rmses_ok <- ! is.null(object$par_RMSEs))  { # more complete object, supersede CIobject which was used to construct it
        if (rmses_ok <- is.null(wrn <- object$par_RMSEs$warn)) {  
          par_RMSEs <- get_from(object,"par_RMSEs")
          if ( rmses_ok <- ! is.null(par_RMSEs)) {
            cat("*** Interval estimates and RMSEs ***\n")
            print(get_from(object,"par_RMSEs"))
          } # NULL if no non-NA CIs... or if no attempt to compute it
        } else message(wrn)
      } 
      if ( ! rmses_ok) {  
        if (is.null(CIwrn <- object$CIobject$warn)) {      
          bounds <- .extract_intervals(object,verbose=FALSE) 
          if (length(bounds)) {
            cat("*** Interval estimates ***\n")
            print(bounds)
          }
        } else message(CIwrn)
      }
    }  
  } else {
    cat("SLik object created. Use MSL(.) to obtain point estimates and CIs.\n")
  }
  invisible(object)
}

print.SLik <- function(x,...) {summary.SLik(x,...)}

calc.lrthreshold.SLik <- function(object,dlr=NULL,verbose=interactive(),...) {
  if (is.null(dlr)) dlr <- Infusion.getOption("LRthreshold") # (<0)
  probErr <- object$MSL$predVar *2.326348 # qnorm(0.99,0,1) ## $predVar differs in conception from the pred MSE used in migraine
  if (verbose) {
    cat("Default d(logL) Chi-square threshold and probable prediction error:\n")
    locstring <- paste(.prettysignif(-dlr)," and ",.prettysignif(probErr),"\n",sep="")
    cat(locstring)
  }
  ## expand beyond *predicted* dlr threshold  ## this should be disconnected from GV$hullExpandFactor
  return( object$MSL$maxlogL + dlr*1.2 -probErr )
}




