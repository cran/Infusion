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
                 parameters="ANY",
                 strategy="ANY") # missing slots => validObject() warning may be hard to track ./.
         # This was solved by setting options(warn=1) and applying str() to slik objects in the .GlobalEnv:
         # the warning appeared when str() the $jointdens and $pardens in objects created by an old version. 
)

.get_mixModstrategy <- function(nc, initParam=NULL, # For programming purposes: 
                                          # optional object inheriting from (Rmixmod) class Parameter
                                          # (typically of class GaussianParameter)
                                ...) { # nc may differ from 2 npar (e.g., by latent variables)
  arglist <- list()
  #### (1) apply Infusion defaults (as they may differ from Rmixmod defaults)
  # ___F I X M E____  rethink, but these defaults may be reasonable:
  ## nc>20 typically meaning that >10 parameters are fitted    
  if (nc>20L) arglist$algo <- c("CEM","EM") # generalize ? CEM before EM => more consistent maxim, but probably slower. try further control
  ## ...otherwise default $algo in "EM".
  #
  global_strategy_args <- Infusion.getOption("global_strategy_args")
  initMethod <- global_strategy_args$initMethod # default is NULL, in which case initParam determines initMethod:
  if (is.null(initMethod) && inherits(initParam,"Parameter")) initMethod <- "parameter"
  # Default method in the absence of an initParam:
  if ( ! inherits(initParam,"Parameter") && is.null(initMethod)) initMethod <- "CEM" # CEM strongly reduces NaN's. "SEMMax" generates manyy errors
  if (initMethod!="parameter") initParam <- NULL # :handling possible conflict between explicit initMethod and explicit initParam
  if ( initMethod =="parameter") { # This must now imply that inherits(initParam,"Parameter")
    arglist$initMethod <- "parameter"
    arglist$parameter <- initParam
  } else arglist$initMethod <- initMethod 
  #
  #### (2) apply User non-defaults (default: empty list ; but user may have modified them):
  arglist[names(global_strategy_args)] <- global_strategy_args
  #
  #### (3) apply hypothetical per-fit controls through \dots :
  dotlist <- list(...)
  arglist[names(dotlist)] <- dotlist
  #
  do.call(Rmixmod::mixmodStrategy, arglist) # returns a Strategy object
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

.any_mixmodResult_issue <- function(mmRe_obj) { # may be a dMixmod object
  chk <- inherits(mmRe_obj,"try-error")
  inherits(mmRe_obj,"try-error") || 
    anyNA(mmRe_obj@parameters@mean) ||
    ! length(mmRe_obj@nbCluster)
}


.any_mixmodCluster_w_retries_issue <- function(mmCl_resu) {
  inherits(mmCl_resu,"try-error") || 
    anyNA((bestResult <- mmCl_resu@bestResult)@parameters@mean) ||
    (! length(bestResult@nbCluster)) ||
    min(table(bestResult@partition))<2L 
  # bestResult once had a cluster with a single element. cov matrix was ==0
  # _____F I X M E_____ => perform tests on the cov matrices [compare determinants?] ? or increase threshold ? 
}

.any_mixmodCluster_issue <- function(mmCl_resu, arglist) {
  checks <- logical(3L) # this is rep(FALSE, 3L)
  chk <- ( 
    (checks[1] <- inherits(mmCl_resu,"try-error")) || 
    (checks[2] <- anyNA(mmCl_resu@bestResult@parameters@mean)) ||
    (checks[3] <- ! length(mmCl_resu@bestResult@nbCluster))) # there are case where the results are thus empty (but in cwhich context?)
  if (chk && .Infusion.data$options$mixturing_errorfn()) {
    reason <- c("try-error", "NaNs", "zero-length parameters")[which(checks)]
    # compact display of nbCluster argument:
    nbClu_info <- arglist$nbCluster
    if (length(nbClu_info)>1L) {
      rge <- range(nbClu_info)
      if (identical(nbClu_info,seq(rge[1],rge[2]))) {
        nbClu_info <- paste0(rge[1]," to ",rge[2])
      } else nbClu_info <- paste0(arglist$nbCluster,collapse=",")
    }
    cat(
      cli::col_green(
        paste0("mixmodCluster(nbCluster=", nbClu_info, ") failed (symptom:",
               reason,") for initMethod=", arglist$strategy@initMethod,
               " and algo=",paste0(arglist$strategy@algo,collapse="+"),"\n")))
  }
  chk
}

.mixmodCluster_w_retries <- function(data, locarglist) {
  resu <- try(.do_call_wrap("mixmodCluster", locarglist), silent = TRUE)
  
  if (.any_mixmodCluster_issue(resu, locarglist)) {
    locdefault <- locarglist$strategy@algo
    if (identical(locdefault,"EM")) {
      alt_algo <- c("CEM","EM")
      alt_algo2 <- c("SEM","EM")
    } else if (identical(locdefault,c("CEM","EM"))) {
      alt_algo <- c("EM")
      alt_algo2 <- c("SEM","EM")
    } else if (identical(locdefault,c("SEM","EM"))) {
      alt_algo <- c("EM")
      alt_algo2 <- c("CEM","EM")
    } 
    locarglist$strategy <- .Infusion.data$options$get_mixModstrategy(
      nc=ncol(data), algo=alt_algo)
    resu <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
  }
  if (.any_mixmodCluster_issue(resu, locarglist)) {
    locarglist$strategy <- .Infusion.data$options$get_mixModstrategy(
      nc=ncol(data), algo=alt_algo2)
    resu <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
  }
  
  if (.any_mixmodCluster_issue(resu, locarglist)) {
    # bestResult full of NaN -> it seems irrespective of mixModstrategy(algo) but dependent on 
    # mixModstrategy(initMethod)
    alt_initMethod <- setdiff(c("CEM","SEMMax"), locarglist$strategy@initMethod)[1L]
    locarglist$strategy <- .Infusion.data$options$get_mixModstrategy(
      nc=ncol(data),initMethod=alt_initMethod)
    resu <- try(.do_call_wrap("mixmodCluster",locarglist),silent = TRUE)
    # (___F I X M E___?) I could nest another call to .mixmodCluster_w_retries() retrying all algo
    # but I must avoir an infinite loop on initMethod...
  }
  
  .any_mixmodCluster_issue(resu,locarglist)
  resu
}

# Alternative to (mixmodCluster+selection by AIC), returns a dMixmod object 
.densityMixmod <- function(
    data, stat.obs, nbCluster=NULL, 
    models=.do_call_wrap("mixmodGaussianModel",
                         list(listModels=Infusion.getOption("mixmodGaussianModel"))), 
    seed=Infusion.getOption("mixmodSeed"),
    initParam=NULL) {
  # stat.obs useful for boundaries and to handle degenerate distributions:
  if (is.null(nbCluster)) nbCluster <- get_nbCluster_range(projdata=data)
  ## 
  blob <- .process_boundaries(data, stat.obs, boundaries=attr(stat.obs,"boundaries")) # original role probably irrelevant
  #     but blob is used 
  if ( ! nrow(blob$data)) { # ie if length(blob$goodrows)==0L in .process_boundaries()
    ## code has long been invalid, but not tested, here:
    # resu <- new("mixmodCluster",nbCluster=0L) # just does not exist ('m'ixmodCluster is the fn, not the class); better syntax:
    ## new(structure("MixmodCluster",package="Rmixmod"),nbCluster=0L) # class is OK but data argument missing, nbCluster=0 never works anayway
    ## next line was 
    # class(resu) <- structure("dMixmod", package="Infusion")
    ## so let's try
    resu <- new(structure("dMixmod", package="Infusion")) # line individually OK; no 'data' nor 'nbCluster' slots
  } else {
    data <- .check_data_post_boundaries(blob$data) # removes lines with NA's
    ## 
    if ((nc <- ncol(data))>(nr <- nrow(data))) {
      stop(paste0("Clustering expected to fail on data with more columns (",nc,
                  " retained parameters and/or statistics)) than rows (",nr," retained points)."))
    }
    if (length(nbCluster)==1L && inherits(initParam,"Parameter") && nbCluster==length(initParam@proportions)) {
      locstrategy <- .Infusion.data$options$get_mixModstrategy(initParam=initParam, nc=ncol(data))
    } else locstrategy <- .Infusion.data$options$get_mixModstrategy(nc=ncol(data))
    # Contra the doc, seed must be a mixmodCluster() argument, not a mixmodStrategy() argument
    # Using it as mixmodStrategy() arg will generate a warning, although it may actually be taken into account...
    # But but... the seed is not effective if the initMethod is "parameter"!
    # To enforce use of a specific seed value, one must therefore use something like
    # Infusion.options(mixmodSeed=789, global_strategy_args=list(initMethod="CEM"))
    locarglist <- list(data=as.data.frame(data), nbCluster=nbCluster, models=models, seed=seed, 
                       strategy=locstrategy)
    #  x <- mixmodCluster(as.data.frame(data),seq(2*log(nrow(data))),strategy = mixmodStrategy(seed=123))
    # statdoc of mixmod recommends ceiling(nrow(data)^0.3) and refers to Bozdogan93
    # But...... it seems to underfit, hence use larger nbCluster range
    resu <- .mixmodCluster_w_retries(data=data, locarglist=locarglist)
    # resu should be a MixmodCluster object, with an element@bestResult of class mixmodResult
    if (length(nbCluster)==1L) {
      if (.any_mixmodCluster_w_retries_issue(resu)) {
        # changing the seed seems more efficient than changing the algorithm (!)
        # here a simple RNG with short period should be enough. =>
        # 'random0' - (short period in https://en.wikipedia.org/wiki/Linear_congruential_generator#Parameters_in_common_use
        locarglist$seed <- (8121L * locarglist$seed + 28411L) %% 134456L 
        if (locarglist$strategy@initMethod=="parameter") locarglist$strategy@initMethod <- "CEM" # otherwise seed is ignored.
        resu <- .mixmodCluster_w_retries(data=data, locarglist=locarglist)
      }
      if (.any_mixmodCluster_w_retries_issue(resu)) { # Then we try any smaller size.
        for (nb in rev(seq_len(min(nbCluster)-1L))) {
          locarglist$nbCluster <- nb
          x <- .mixmodCluster_w_retries(data=data, locarglist=locarglist)
          if ( ! .any_mixmodCluster_w_retries_issue(x) ) break
        }
        resu <- x@bestResult # MixmodResults, not mixmodCluster
        class(resu) <- structure("dMixmod", package="Infusion") 
        resu@strategy <- x@strategy
      } 
    } else {
      ## in that case, calling .any_mixmodCluster_w_retries_issue() would still checks @bestResult, 
      ## which may be invalid (possibly irrelevant), and  
      ## .get_best_mixmod_by_IC(resu) is able to extract other, possibly valid or invalid models
      ## =>  This means .get_best_mixmod_by_IC() should check and reject any invalid model,
      ## so that the .get_best_mixmod_by_IC() result below is correct.
      ## I finally got an example from test-reparam -> ini_rpdensv <- infer_SLik_joint(projrpSimuls,stat.obs=projrpSobs)
      ## ... in the formal checks only!
    }
    ## Rmixmod::mixmodCluster and Infusion::dMixmod objects have strategy info, but Rmixmod::mixmodResults don't.
    ## And the finally succesful strategy may differ from the locarglst$ one.
    ## => Awkward copying with conditional initialization: 
    ## 'strategy <- ...' needed whenever we get a mixmodCluster; '... <- strategy' whenever we create a dMixmod.
    if ( inherits(resu,"MixmodCluster")) {
      strategy <- resu@strategy
      resu <- .get_best_mixmod_by_IC(resu) # MixmodResults
    } 
    if ( inherits(resu,"MixmodResults")) {
      class(resu) <- structure("dMixmod", package="Infusion") 
      resu@strategy <- strategy
    } 
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
  } 
  resu@varNames <- blob$varNames # this works on dMixmod objects, not on Rmixmod:: classes
  if (is.null(blob$statNames)) {
    resu@statNames <- NA_character_
  } else resu@statNames <- blob$statNames
  resu@simuls_activeBoundaries <- blob$simuls_activeBoundaries
  resu@Sobs_activeBoundaries <- blob$Sobs_activeBoundaries
  resu@freq <- blob$freq
  resu
}

.fast_as.matrix.df <- function(newdata, varNames) {
  
  ## This is faster than the often recommended [e.g. https://csgillespie.github.io/efficientR/performance.html] data.matrix():
  # newdata <- matrix(unlist(newdata,recursive = FALSE, use.names = FALSE), nrow=nrow(newdata),
  #                  dimnames=list(NULL,colnames(newdata)) ) # newdata <- as.matrix(newdata)
  ## But even better:
  newdata <- newdata[,varNames,drop=FALSE]
  nrnc <- dim(newdata)
  newdata <- unlist(newdata,recursive = FALSE, use.names = FALSE)
  dim(newdata) <- nrnc
  colnames(newdata) <- varNames
  newdata
}


predict.dMixmod <- function(object,
                            newdata, 
                            tcstat.obs=NULL, ## to avoid checks of arguments
                            solve_t_chol_sigma_list,
                            logproportions, # ie log(object@parameters["proportions"])
                            clu_means,
                            log=FALSE, ...) {
  varNames <- object@varNames
  if (is.null(tcstat.obs)) {
    if (nodim <- is.null(dim(newdata))) { ## less well controlled case, but useful for maximization (which is not performed in canned procedures)
      if ((ns <- length(varNames)) != length(newdata)) {
        stop(paste("(!) newdata has incorrect length. It should match the following variables:\n",
                   paste(varNames,collapse=" ")))
      } else {
        datanames <- names(newdata) # before they are erased by the next dim() assignment (!)
        dim(newdata) <- c(1L, ns)
        if (nonames <- is.null(datanames)) {
          colnames(newdata) <- varNames
        } else colnames(newdata) <- datanames
      }
    } else if (isDF <- inherits(newdata,"data.frame")) {
      newdata <- .fast_as.matrix.df(newdata, varNames)
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
    } else if ( ( ! (nodim && nonames)) ||
                  ! (nodim || isDF)) {  # data fram already sus-col-ed
      densitydata <- newdata[,varNames,drop=FALSE]
    } else densitydata <- newdata
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
          logdensity[k] <- logproportions[k] + 
            .fast_dmvnorm(densitydata, clu_means[, k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=TRUE)
        }
        mixture <- matrixStats::logSumExp(logdensity)
      } else {
        logdensity <- matrix(nrow=nr,ncol=nbCluster)
        for (k in 1:nbCluster) {
          logdensity[,k] <- logproportions[k] + 
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
        density[,k] <- exp(logproportions[k])*
          .fast_dmvnorm(densitydata, clu_means[, k], solve_t_chol_sigma= solve_t_chol_sigma_list[[k]],log=FALSE)
      }
      mixture <- .rowSums(density,m=nrow(density), n=ncol(density)) ## sum(Li) 
      if ( !is.null(Sobs_activeBoundaries)) mixture <- mixture*freqs
    } else mixture <- freqs
  }
  return(mixture)
}

.grad.dMixmod <- function(object,
                            newdata, 
                            solve_t_chol_sigma_list,
                            logproportions,
                            clu_means,
                         nbCluster,
                         col_ids,
                            log=TRUE, ...) {

  
  if (log) {
    logdensity <- matrix(nrow=1L,ncol=nbCluster)
    grads <- matrix(nrow=length(col_ids),ncol=nbCluster)
    for (k in 1:nbCluster) {
      sig_k <- solve_t_chol_sigma_list[[k]]
      logdensity[,k] <- logproportions[k] +  
        .fast_dmvnorm(newdata, clu_means[, k], solve_t_chol_sigma= sig_k,log=TRUE)
      gr <- tcrossprod((clu_means[, k] - newdata), sig_k)
      grads[,k] <-  (gr %*% sig_k)[col_ids]
    }
    logdensity <- logdensity - .logsumexp(logdensity)
    
    grad <- rowSums(sweep(grads, MARGIN=2L, exp(logdensity), `*`))
    
  } else {
    stop("code missing here")

  }
  return(grad)
}

# used as if (inherits(currMSL$jointdens, "dMixmod")) plot(currMSL$jointdens, data=currMSL$logLs) 
# in the 1D narrow-peak example.
plot.dMixmod <- function(x, data=NULL, 
                         vars = colnames(data)[1:2],
                         pch= x@partition %% 20, # Rmixmod code did not anticipate more than 20 clusters...
                         ...) {
  if (requireNamespace("Rmixmod")) {
    if ( ! inherits(x, "dMixmod")) warning("class of 'x' is suspect; something like <SLik_j object>$jointdens is expected here")
    class(x) <- "MixmodResults"
    if (is.null(data)) stop("'data' needed: something like <SLik_j object>$logLs")
    varnames <- intersect(x@varNames, colnames(data)) # notably to remove constant parameters...
    Rmixmod::plotCluster(x=x, data=data[, varnames, drop=FALSE], 
                         variable1=vars[1], variable2=vars[2],
                         pch=pch, ...)
  } # else we are unlikely to need this function (which is used only for debugging)...
}

.wrap_rmvnorm <- function(nsim, parameters, clu_i) {
  rmvnorm(n=nsim, 
          mean=parameters["mean", clu_i], 
          sigma= parameters["variance",clu_i])
}

.wrap_rmvt <- function(nsim, parameters,clu_i) {
  rmvt(n=nsim, 
       delta= parameters["mean", clu_i], 
       sigma= parameters["variance",clu_i])
}


.simulate.MixmodResults <- function (object, seed=NULL, 
                                     size=1, # number of points for each simulation ~user-level nsim
                                     drop=TRUE,
                                     norm_or_t, # .wrap_rmvnorm or .wrap_rmvt
                                     n_tables=1L, 
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
  onesimfn <- function() { # 'returns 'size' vectors
    rclu <- sample(seq( object@nbCluster),size,replace=TRUE, prob=prob) ## vector of sampled clusters
    rclutable <- table(rclu) # samples clusters -> rclutable is # ofsampled points from eachcluster
    onesim <- vector("list", length(rclutable))
    names(onesim) <- clu_ids <- names(rclutable)
    for (clu_it in clu_ids) {
      onesim[[clu_it]] <- norm_or_t(nsim=rclutable[clu_it], 
                                    parameters=object@parameters, clu_i=as.numeric(clu_it))
    }
    onesim <- do.call(rbind,onesim)
    return(onesim)
  }
  simuls <- replicate(n_tables,onesimfn(),simplify = FALSE) # typically n_tables=1 and size is large
  if (drop && n_tables==1L) simuls <- simuls[[1L]]
  colnames(simuls) <- colnames(object@parameters["mean",])
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

