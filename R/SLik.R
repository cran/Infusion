predict.SLik <- function (object, newdata=object$logLs[,object$colTypes$fittedPars,drop=FALSE],...) {
  return(predict(object$fit,newdata=newdata,...))
}

refine.SLik <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "REML"
  refine.default(object,surfaceData=object$logLs,method=method,...)
}

# options("digits") controls digits in print()
## 3 digits on RMSEs makes a lot of digits overall   
summary.SLik <- function(object, ...) { 
  if ( !is.null(object$MSL) ) {
    cat(paste("*** Summary ML (",max(object$logLs$cumul_iter)," iterations, ",nrow(object$logLs)," points): ***\n",sep=""))
    print(c(object$MSL$MSLE,"logL"=object$MSL$maxlogL,"RMSE_logL"=unname(get_from(object,"RMSEs")[1L])))
    # # par_RMSEs
    if ( ! is.null(object$CIobject))  {
      if( ! is.null(object$par_RMSEs))  { # more complete object, supersede CIobject which was used to construct it
        if (is.null(wrn <- object$par_RMSEs$warn)) {      
          cat("*** Interval estimates and RMSEs ***\n")
          print(get_from(object,"par_RMSEs"))
        } else message(wrn)
      } else if (is.null(CIwrn <- object$CIobject$warn)) {      
        bounds <- .extract_intervals(object,verbose=FALSE) 
        if (length(bounds)) {
          cat("*** Interval estimates ***\n")
          print(bounds)
        }
      } else message(CIwrn)
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