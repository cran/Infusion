def_projectors <- function(reftable, pars, # latent vars ?
                            stats, npp=1L,
                            projNames=NULL, 
                           npp_opt=Infusion.getOption("npp_opt"), 
                           ...) {
  npar <- length(pars)
  projectors <- vector("list", npar*npp)
  if (length(npp)==1L) npp <- rep(npp, npar)
  if (is.null(projNames)) {
    if (any(npp!=1L)) {
      projNames <- pars[rep(seq(npar),npp)]
      projNames <- paste0(projNames,".", unlist(sapply(npp,seq)))
    } else {
      projNames <- paste0("p",pars)
    }
  }
  if (length(intersect(projNames,pars))) stop("Invalid 'projNames' (at least on is identical to a parameter name).")
  if (length(unique(projNames))!=length(projNames)) stop("Invalid 'projNames' (not unique values)")
  if (length(projNames)!=length(projectors)) stop("Invalid 'projNames' (length(projNames)!=length(projectors))")
  names(projectors) <- projNames
  it <- 0L # iterator over all projections over all parameters
  for (parit in seq_along(pars)) {
    respname <- parname <- pars[parit]
    parent <- NULL
    locdata <- reftable 
    locstatnames <- stats
    quantreg <- npp[parit]>1L && npp_opt=="quantiles"
    # max4par
    for (itpp in seq_len(npp[parit])) {
      if (itpp>1L) { 
        parent <- projNames[it]
        if (FALSE) {# concept: cf email 2024/08/30 12:14
          # response and respname unchanged, predictor added
          newname <- .makenewname("pred",names(locdata))
          parentpreds <- .predictWrap(projectors[[parent]], newdata=locdata, is_trainset=TRUE)
          locdata[,newname] <- parentpreds
          locstatnames <- c(locstatnames,newname)
        } else if (npp_opt=="quantiles") {
          # projector unchanged
        } else {
          # D:/home/francois/travail/stats/Infusionplus/caseStudies/physes/main_v2/essai_doubleRF.R
          # new response, same predictors: easier wince predict does not need new predictors.
          respname <- .makenewname(paste0("res",parname),names(locdata))
          parentpreds <- .predictWrap(projectors[[parent]], newdata=locdata, is_trainset=TRUE)
          locdata[,respname] <- locdata[,parname] - parentpreds
          projector <- project(respname, stats=locstatnames, data=locdata, 
                                         quantreg=quantreg, ...)
        }
      } else projector <- project(respname, stats=locstatnames, data=locdata, 
                                            quantreg=quantreg, ...)
      it <- it+1L
      if (itpp>1L) {
        projector <- list(projector=projector,
                          npp_opt=npp_opt,
                          quantile=0.67,
                          respname=respname,
                          parname=parname,
                          parent=parent) # parent projector's name
        class(projector)  <- c("projector",class(projector))
      }  
      projectors[[it]] <- projector
    }
  }
  projectors
}
