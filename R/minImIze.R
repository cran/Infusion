.IC_wrapper <- function(i, data, models, seed=123, which=Infusion.getOption("criterion")) {
  cluObject <- suppressWarnings(.do_call_wrap("mixmodCluster",list(data=as.data.frame(data), nbCluster=i, models=models, 
                                                                   seed=seed, strategy=eval(Infusion.getOption("strategy"))) ))
  if (cluObject@error) {
    return(list(obj=NA, cluObject=NA))
  } else {
    BIC <- cluObject@results[[1]]@criterionValue
    if (which=="BIC") return(list(obj=BIC, cluObject=cluObject))
    logL <- cluObject@results[[1]]@likelihood
    df <- (2*logL+BIC)/(log(cluObject@nbSample))
    AIC <- -2*logL+2*df
    return(list(obj=AIC, cluObject=cluObject))
  }
}

.minImIze <- function(init=NULL, objfn, range, trace=FALSE, ...) {
  fneval <- 0L
  lower <- min(range)
  upper <- max(range)
  range <- upper - lower
  if (is.null(init)) init <- floor(lower+range/2L)
  if (trace) print(c(lower=lower,upper=upper,fneval=fneval,init=init))
  objlist <- list()
  while (range>0L) {
    vcinit <- as.character(init)
    if (is.null(objlist[[vcinit]])) {
      objlist[[vcinit]] <- objfn(init, ...)
      fneval <- fneval+1L
    }
    objinit <- objlist[[vcinit]]$obj
    #
    vp <- min(upper, init+1L)
    vc <- as.character(vp)
    if (is.null(objlist[[vc]])) {
      objlist[[vc]] <- objfn(vp, ...)
      fneval <- fneval+1L
    }
    objp <- objlist[[vc]]$obj
    #
    if (objp<objinit) {
      lower <- vp
    } else {
      upper <- init
      vm <- max(lower, init-1L)
      vc <- as.character(vm)
      if (is.null(objlist[[vc]])) {
        objlist[[vc]] <- objfn(vm, ...)
        fneval <- fneval+1L
      }
      objm <- objlist[[vc]]$obj
      #
      if (objm<objinit) {
        upper <- vm
      } else lower <- init
    }
    #
    range <- upper - lower
    init <- floor(lower+range/2L)
    if (trace) print(c(nextlower=lower,nextupper=upper,fneval=fneval,nextinit=init))
  }
  resu <- objlist[[as.character(lower)]]
  resu$solution <- init
  resu$fneval <- fneval
  return(resu)
}

.minimize_from_upper <- function(objfn, range, ...) {
  lower <- min(range)
  nb_clu <- max(range)
  currbest <- objfn(nb_clu, ...)
  while (nb_clu>lower) {
    objblob <- objfn(nb_clu-1L, ...)
    if (objblob$obj<currbest$obj) {
      currbest <- objblob
      nb_clu <- nb_clu-1L
    } else break
  }
  currbest$solution <- nb_clu
  return(currbest)
}
