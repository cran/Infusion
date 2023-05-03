setClassUnion("NULLorChar", c("NULL","character"))
setClassUnion("NULLorNum", c("NULL","numeric"))

### A class to describe distributions that involve discrete probability masses for some variables and gaussian mixtures
## for other variables conditional on such discrete events. We aim to reuse "MixmodResults" with minimal overhead.
## Formally extending "MixmodResults"  would requires setting an explicit dependency to the Rmixmod package:
# setClass("dMixmod", contains = c("MixmodResults"), slots=c(<first slots below>).
## Instead, the dMixmod class uses S4 and has all the slots of the \code{MixmodResults} class, but does not formally extends it.
## Note the non-standard declaration below, and the non-standard construction by
# class(<mixmodCluster object>) <- structure("dMixmod", package="Infusion")
## But NOT an attempted extension by
#class(resu) <- c("dMixmod", class(resu)) # => Setting class(x) to multiple strings ("dMixmod", "MixmodResults", ...); result will no longer be an S4 object
setClass("dMixmod", 
         slots=c(freq="NULLorNum",
                 varNames="character",
                 statNames="character",
                 simuls_activeBoundaries="NULLorNum",
                 Sobs_activeBoundaries="NULLorChar",
                 # There are other slots inherited from mixmodResults(mixmodCluster?) objects, but
                 # we need to declare here those that are *set* (LHS, not RHS) by some Infusion fns (.marginalize_dMixmod):
                 proba="matrix",
                 parameters="ANY")
)

# Approach assuming that the option is a list of arguments (which is currently not the case):
.get_mixModstrategy <- function(strategy=Infusion.getOption("strategy"), seed=NULL) {
  if (is.null(strategy$seed)) strategy$seed <- seed
  # => "WARNING: Strategy@seed attribute is obsolete! Consider using mixmodCluster function's 'seed' parameter instead"
  .do_call_wrap("mixmodStrategy",arglist=strategy,pack="Rmixmod")
}

# note the two possible return types:
.process_boundaries <- function(data, stat.obs, boundaries) {
  varNames <- colnames(data) # may be modified below
  statNames <- names(stat.obs)  # colnames(data) was oK for the (2017) version of Infusion without modelling of joint distribs
  parameterNames <- setdiff(colnames(data),statNames)
  Sobs_activeBoundaries <- simuls_activeBoundaries <- freq <- NULL
  if ( ! is.null(boundaries) ) {
    simulsMatches <- t(abs(apply(data[,names(boundaries),drop=FALSE],1L,`-`,y=boundaries))<1e-14)
    if (length(boundaries)==1L) { ## apply() fix
      simulsMatches <- t(simulsMatches)
      colnames(simulsMatches) <- names(boundaries)
    }
    simulsAtSomeBoundaries <- matrixStats::rowAnys(simulsMatches)
    SobsMatches <- (stat.obs[names(boundaries)]==boundaries)
    SobsAtSomeBoundaries <- any(SobsMatches)
    if (any(SobsAtSomeBoundaries)) {
      Sobs_activeBoundaries <- boundaries[which(SobsAtSomeBoundaries)]
      activeNames <- names(Sobs_activeBoundaries)
      ## keep only simuls that match all the boundaries that stat.obs matches:
      goodrows <- which(apply(simulsMatches[,activeNames,drop=FALSE],1L,all))
      # This 'goodrows' may have length zero => special handling in calling function
      varNames <- c(parameterNames,setdiff(statNames,activeNames))  # modify previous varNames
      ## estimate proba that all (SobsMatches) are observed in the simulation 
      freq <- (length(goodrows)+1/2)/(nrow(data)+1) ## log(freq) is estimator of log(p) with bias ~1/n^2 
      # and select data to perform mixture modelling on remaining statistics (for (2017) method)
      #     or     data to perform mixture modelling on parameters and remaining statistics (for SLik_j) 
      data <- data[goodrows,varNames,drop=FALSE] ## conditional on match, lower-dimensional
      # For SLik_j, given stat.obs at boundaries, we thus shoult compute joint as 
      #  freq * joint density(param and other stats| simul at boundary)
    } else if (any(simulsAtSomeBoundaries)) {
      simuls_activeBoundaries <- boundaries[which(apply(simulsMatches,2L,any))]
      ## keep only simuls that match none of the the boundaries:
      goodrows <- which( ! apply(simulsMatches[,names(simuls_activeBoundaries),drop=FALSE],1L,any)) 
      ## estimate proba that no (simulsMatches) is observed in the simulation 
      freq <- (length(goodrows)+1/2)/(nrow(data)+1) ## log(freq) is estimator of log(p) with bias ~1/n^2 
      data <- data[goodrows,,drop=FALSE] ## conditional on no match, full-dimensional
      # For SLik_j, given stat.obs not at boundaries, we thus shoult compute joint as 
      #  freq * joint density(all param and stats| simul not at boundary)
    }
  }
  resu <- list(data=data, freq=freq, simuls_activeBoundaries=simuls_activeBoundaries, statNames=statNames,
               Sobs_activeBoundaries=Sobs_activeBoundaries, varNames=varNames) # alternatively, a mixmodCluster object has been returned
} 

# Alternative to (mixmodCluster+selection by AIC), returns a dMixmod object 
.densityMixmod <- function(data,stat.obs,nbCluster=seq_nbCluster(nr=nrow(data)), 
                           models=.do_call_wrap("mixmodGaussianModel",list(listModels=Infusion.getOption("mixmodGaussianModel"))), 
                           seed=Infusion.getOption("mixmodSeed")) {
  # stat.obs useful for boundaries and to handle degenerate distributions:
  nbCluster <- get_nbCluster_range(projdata=data, nbCluster=nbCluster)
  ## 
  blob <- .process_boundaries(data, stat.obs, boundaries=attr(stat.obs,"boundaries"))
  if ( ! nrow(blob$data)) { # ie if length(blob$goodrows)==0L in .process_boundaries()
    ## code has long been invalid, but not tested, here:
    # resu <- new("mixmodCluster",nbCluster=0L) # just does not exist ('m'ixmodCluster is the fn, not the class); better syntax:
    ## new(structure("MixmodCluster",package="Rmixmod"),nbCluster=0L) # class is OK but data argument missing, nbCluster=0 never works anayway
    ## next line was 
    # class(resu) <- structure("dMixmod", package="Infusion")
    ## so let's try
    resu <- new(structure("dMixmod", package="Infusion")) # line individually OK; no 'data' nor 'nbCluster' slots
  } else {
    # data lines matching the 'boundaries' values have now been removed. There may still be NA's in the data (notably, lines of NA's). The doc explains
    # "the inference method ... must thus ignore all empirical distributions including NA/NaN/inf.": 
    data <- .check_data_post_boundaries(blob$data)
    ## 
    if (inherits(data,"data.frame")) {
      checkranges <- sapply(data,range)
    } else checkranges <- t(matrixStats::colRanges(data))
    rangewidths <- checkranges[2L,]-checkranges[1L,]
    degenerates <- rangewidths<1e-6
    if (any(degenerates)) {
      which_dg <- which(degenerates)
      which_out <- (stat.obs[which_dg]>checkranges[2L,which_dg]) | (stat.obs[which_dg]<checkranges[1L,which_dg])
      if (any(which_out)) {
        for (dg in which_dg) {
          dst <- max(abs(stat.obs[dg]-data[,dg]))
          data[,dg] <- data[,dg] + rnorm(n=nrow(data),sd=dst/5) ## makes it reasonably unlikely that obs will be within noised data
        }
      } else  {
        ## stat.obs in degenerate simuls but mixmodCluster needs non-degenerate distribution 
        for (dg in which_dg) data[,dg] <- data[,dg] + rnorm(n=nrow(data),sd=1e-06)
      }
    }
    if ((nc <- ncol(data))>(nr <- nrow(data))) {
      stop(paste0("Clustering expected to fail on data with more columns (",nc,
                  " retained parameters and/or statistics)) than rows (",nr," retained points)."))
    }
    #  x <- mixmodCluster(as.data.frame(data),seq(2*log(nrow(data))),strategy = mixmodStrategy(seed=123))
    # statdoc of mixmod recommends ceiling(nrow(data)^0.3) and refers to Bozdogan93
    # But...... it seems to underfit, hence use larger nbCluster range
    clu_optimizer <- Infusion.getOption("clu_optimizer")
    if (clu_optimizer==".minImIze") {
      resu <- try(.minImIze(init=NULL, objfn=.IC_wrapper, range=nbCluster, trace=FALSE, data=data, models=models, 
                            seed=123)$cluObject@bestResult,silent=TRUE) # F I X M E ? $cluObject@bestResult loses info about $cluObject@strategy
    } else if (clu_optimizer==".minimize_from_upper") {
      resu <- try(.minimize_from_upper(objfn=.IC_wrapper, range=nbCluster, trace=FALSE, data=data, models=models, 
                                       seed=123)$cluObject@bestResult,silent=TRUE) # F I X M E ? $cluObject@bestResult loses info about $cluObject@strategy
    } else { ## default !
      resu <- try(.do_call_wrap("mixmodCluster",list(data=as.data.frame(data), nbCluster=nbCluster, models=models, seed=seed, 
                                                     strategy=eval(Infusion.getOption("strategy")))),
                  silent = TRUE)
    }
    if (inherits(resu,"try-error")) {
      x <- suppressWarnings(.do_call_wrap("mixmodCluster",list(data=as.data.frame(data), nbCluster=nbCluster, models=models, 
                                                               seed=123) ))
      if (length(nbCluster)==1L) {
        if (x@error) {
          if (FALSE) {
            message(paste("mixmodCluster() failed for nbCluster=",
                          paste(nbCluster,collapse=" ")))
            message("Trying decreasing nbCluster...")
          }
          for (nb in rev(seq_len(min(nbCluster)-1L))) {
            x <- suppressWarnings(.do_call_wrap("mixmodCluster",list(data=as.data.frame(data), nbCluster=nb, models=models, 
                                                                     seed=seed) ))
            if ( ! x@error ) break
          }
        } 
        resu <- x@bestResult # MixmodResults
      } 
    } 
    # if (inherits(resu,"try-error")) {
    #   dump.frames("C:/home/francois/travail/stats/InferentialSimulation/Infusion/devtoolsdump",to.file = TRUE)
    # } else 
    if ( inherits(resu,"MixmodCluster")) resu <- .get_best_mixmod_by_IC(resu) # MixmodResults
    if (resu@model=="Gaussian_pk_Lk_C") { 
      ## then cov matrices are proportional, wecheck they are not too heterogeneous
      varmats <- resu@parameters@variance
      vars11 <- unlist(lapply(varmats,`[`,x=1L,y=1L))
      rangelv <- range(log(vars11))
      if ( rangelv[2L]-rangelv[1L]> 30) {
        message("Inferred Gaussian components have very heterogeneous variances.")
        message("This can occur when some summary statistic takes one particular value with high probability,")
        message("which will be poorly fitted by mixmodCluster().")
        warning("Result of densityMixmod with model \"Gaussian_pk_Lk_C\" is suspect.\n See screen messages for further information")
      } 
    } else if (resu@model=="Gaussian_pk_Lk_Ck"){
      if (resu@error=="determinant of matrix too small" ) {
        message("Inferred Gaussian components have too small determinant.")
        message("This can occur when some summary statistic takes one particular value with high probability,")
        message("which will be poorly fitted by mixmodCluster().")
        warning("Result of densityMixmod with model \"Gaussian_pk_Lk_Ck\" is suspect.\n See screen messages for further information.")
      }
    }
    class(resu) <- structure("dMixmod", package="Infusion") 
  } 
  resu@varNames <- blob$varNames
  resu@statNames <- blob$statNames
  resu@simuls_activeBoundaries <- blob$simuls_activeBoundaries
  resu@Sobs_activeBoundaries <- blob$Sobs_activeBoundaries
  resu@freq <- blob$freq
  resu
}

predict.dMixmod <- function(object,
                            newdata, ## should be (reformatted) stat.obs (in canned procedures)
                            tcstat.obs=NULL, ## to avoid checks of arguments
                            solve_t_chol_sigma_list,
                            logproportions,
                            clu_means,
                            log=FALSE, ...) {
  varNames <- object@varNames
  if (is.null(tcstat.obs)) {
    if (is.null(dim(newdata))) { ## less well controlled case, but useful for maximization (which is not performed in canned procedures)
      if ((ns <- length(varNames)) != length(newdata)) {
        stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",
                   paste(varNames,collapse=" ")))
      } else {
        datanames <- names(newdata) # before they are erased by the next dim() assignment (!)
        dim(newdata) <- c(1L, ns)
        if (is.null(datanames)) {
          colnames(newdata) <- varNames
        } else colnames(newdata) <- datanames
      }
    } else if (inherits(newdata,"data.frame")) {
      newdata <- matrix(unlist(newdata,recursive = FALSE, use.names = FALSE), nrow=nrow(newdata),
                        dimnames=list(NULL,colnames(newdata)) ) # newdata <- as.matrix(newdata)
    }
    if ( !is.null(Sobs_activeBoundaries <- object@Sobs_activeBoundaries)) { # if the Sobs matches some boundaries, identified in Sobs_activeBoundaries
      # then looks whether the newdata match all of the boundaries met by Sobs
      boundsdata <-  newdata[,names(Sobs_activeBoundaries),drop=FALSE]
      atb <- apply(boundsdata,1L,`==`,y=Sobs_activeBoundaries)
      # If not (i.e. only partial matching) the object does not predict correctly the density => warning + heuristic patch:
      if (!all(atb)) {warning("'predict.dMixmod' cannot compute joint out-of-boundary density from conditional at-boundary density. ")}
      freq <- object@freq
      freqs <- atb*freq+(1-atb)*(1-freq) ## uses 1-freq instead of the unknown density of the variable(s) in boundaries 
      densitydata <- newdata[,varNames,drop=FALSE]
    } else if ( !is.null(simuls_activeBoundaries <- object@simuls_activeBoundaries)) {
      ## do not reduce data in this case
      densitydata <- newdata 
      ## only for the warning:
      boundsdata <-  newdata[,names(simuls_activeBoundaries),drop=FALSE]
      atb <- apply(boundsdata,1L,`==`,y=simuls_activeBoundaries)
      if (any(atb)) {
        warning("'predict.dMixmod' cannot compute conditional at-boundary density from joint out-of-boundary density.")
        # return value is the same as for predict(,tcstat.obs=<newdata>) 
      }
    } else densitydata <- newdata[,varNames,drop=FALSE]
  } else { 
    if ( !is.null(Sobs_activeBoundaries <- object@Sobs_activeBoundaries)) {
      atb <- TRUE
      freqs <- object@freq
      densitydata <- tcstat.obs[,varNames,drop=FALSE]
    } else densitydata <- tcstat.obs  ## [,statNames,drop=FALSE]
  } 
  nbCluster <- object@nbCluster
  if (log) { 
    if (nbCluster>0L) {
      nr <- nrow(densitydata)
      if (nr==1L) {
        logdensity <- numeric(nbCluster)
        for (k in 1:nbCluster) {
          logdensity[k] <- logproportions[k] + # log(object@parameters["proportions", k]) + 
            .fast_dmvnorm(densitydata, clu_means[, k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=TRUE)
        }
        maxlogs <- max(logdensity)
        normalizedlogs <- logdensity - maxlogs 
        ## exp(normalizedlogs) in (0,1); sum(exp(logLi-maxlog))= exp(-maxlog) sum(exp(logLi))= exp(-maxlog) sum(Li)
        mixture <- log(sum(exp(normalizedlogs))) + maxlogs ## log(sum(Li))
      } else {
        logdensity <- matrix(nrow=nr,ncol=nbCluster)
        for (k in 1:nbCluster) {
          logdensity[,k] <- logproportions[k] + # log(object@parameters["proportions", k]) + 
            .fast_dmvnorm(densitydata, clu_means[, k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=TRUE)
        }
        mixture <- matrixStats::rowLogSumExps(logdensity)
      }
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture + atb*log(freqs)
    } else mixture <- atb*log(freqs)
  } else {
    if (nbCluster>0L) {
      density <- matrix(nrow=nrow(densitydata),ncol=nbCluster)
      for (k in 1:nbCluster) {
        density[,k] <- logproportions[k] + # object@parameters["proportions", k] * 
          .fast_dmvnorm(densitydata, clu_means[, k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=FALSE)
      }
      mixture <- .rowSums(density,m=nrow(density), n=ncol(density)) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture)
}


plot.dMixmod <- function(x, data=NULL, ...) {
  if (requireNamespace("Rmixmod")) {
    class(x) <- "MixmodResults"
    if (is.null(data)) stop("'data' needed: something like <SLik_j object>$logLs")
    Rmixmod::plotCluster(x=x, data=data, ...)
  } # else we are unlikely to need this function (which is used only for debugging)...
}

.simulate.MixmodResults <- function (object, nsim=1, seed=NULL, 
                                     size=1, # number of points for each simulation 
                                     drop=TRUE,
                                     ...) {
  ## RNG stuff copied from simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)
  if ( ! is.null(seed)) {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  #
  prob <- object@parameters["proportions"]
  locrmvnorm <- function(v,rclutable) {
    numv <- as.numeric(v)
    rmvnorm(rclutable[v], 
            mean=object@parameters["mean", numv], 
            sigma= object@parameters["variance",numv])
  }
  onesimfn <- function() {
    rclu <- sample(seq( object@nbCluster),size,replace=TRUE, prob=prob) ## vector of sampled clusters
    rclutable <- table(rclu)
    onesim <- do.call(rbind,lapply(names(rclutable), locrmvnorm,rclutable=rclutable))
    return(onesim)
  }
  simuls <- replicate(nsim,onesimfn(),simplify = FALSE)
  if (drop && nsim==1L) simuls <- simuls[[1L]]
  return(simuls)
}

.do_call_wrap <- function(chr_fnname,arglist, pack="Rmixmod") { # could be sought as .wrap_do_call
  #eval(as.call(c(quote(require),list(package="Rmixmod", quietly = TRUE))))
  if (length(grep(pack,packageDescription("Infusion")$Imports))) {
    ## then the necessary functions must be imported-from in the NAMESPACE  
    do.call(chr_fnname,arglist) ## "stuff"
  } else if (length(grep(pack,packageDescription("Infusion")$Suggests))) {
    ## then the necessary functions cannot be imported-from in the NAMESPACE  (and the package must be written in an appropriate way)
    if ( requireNamespace(pack, quietly = TRUE)) {
      #eval(as.call(c(chr_fnname,arglist))) # quote(stuff)
      myfun <- get(chr_fnname, asNamespace(pack)) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
      do.call(myfun,arglist) ## "stuff"           ## la version "jeroen" (ibid) est int aussi.
    } else {stop(paste("'",pack,"' required but not available.",sep=""))}
  } else { ## package not declared in DESCRIPTION; to circumvent possible archiving of Rmixmod
    if (suppressWarnings(do.call("require",list(package=pack, quietly = TRUE)))) {
      #eval(as.call(c(chr_fnname,arglist))) # quote(stuff)
      do.call(chr_fnname,arglist) ## "stuff"
    } else {stop(paste("'",pack,"' required but not available.",sep=""))}
  }
}

