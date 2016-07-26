check_logLs_input <- function(object) {
  Naf <- rep(FALSE,length(object))
  for (lit in seq(length(object),1L,-1L)) { 
    if (any( ! is.finite(object[[lit]]))) Naf[lit] <- TRUE ## is.finite false for NA/NaN/Inf
  }
  if(any(Naf)) {
    object[which(Naf)] <- NULL
    message(paste("NA/NaN/Inf found in distribution(s) # ",paste(which(Naf),collapse=" "),".",sep=""))
    message("Any such distribution is removed from the analysis.")
    message("see help(\"handling_NAs\") for more information.")
  }
  return(object)
}