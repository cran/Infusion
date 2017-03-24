predict.SLik <- function (object, newdata=object$logLs[,object$colTypes$fittedPars,drop=FALSE],...) {
  return(predict(object$fit,newdata=newdata,...))
}

refine.SLik <- function(object,method=NULL,...) {
  if (is.null(method)) method <- "REML"
  refine.default(object,surfaceData=object$logLs,method=method,...)
}

print.SLik <-function(x,...) {summary.SLik(x,...)}

# options("digits") controls digits in print()
## 3 digits on RMSEs makes a lot of digits overall   
`summary.SLik` <- function(object, ...) { 
  if ( !is.null(object$MSL) ) {
    cat(paste("*** Summary ML (cumul. iter. =",max(object$logLs$cumul_iter),"): ***\n",sep=""))
    print(c(object$MSL$MSLE,"logL"=object$MSL$maxlogL,"RMSE_logL"=unname(object$RMSEs[1L])))
    #
    if (inherits(object,"SLik_j")) {
      cat("*** Interval estimates ***\n")
      print(object$pars)
    } else if( ! is.null(CIobject <- object$par_RMSEs))  {
      cat("*** Interval estimates and RMSEs ***\n")
      print(object$par_RMSEs)
    } 
  } else {
    cat("SLik object created. Use MSL(.) to obtain point estimates and CIs.\n")
  }
  invisible(object)
}

`summary.SLik_j` <- function(object,...) summary.SLik(object=object,...)

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