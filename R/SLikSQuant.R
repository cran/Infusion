
.calc_EI <- function(summInferObject,points,Qmax=NULL) {  ## for both SLik and SLikp
  trypred <- predict(summInferObject,newdata=points,variances=list(linPred=TRUE,dispVar=TRUE))
  tryVar <- attr(trypred,"predVar")
  if (any(tryVar<0))  { ## anticipating numerical problems (ignore also spuriously large tryVar values)
    return( summInferObject$fit$family$linkfun(as.numeric(trypred[,1])) )
  } else {
    trySE <- sqrt(tryVar)
    #dmudeta <- summInferObject$fit$family$mu.eta(as.numeric(trypred)) ## as.numeric also removes attributes, quite useful
    #tryVar <- tryVar * dmudeta^2  
    ## EI on linear predictor scale vs prediction on response scale
    tryQ <- summInferObject$fit$family$linkfun(as.numeric(trypred[,1])) + 1.96*trySE ## improvement function for candidate points
    if (is.null(Qmax)) Qmax <- max(tryQ)
    return( trySE*dnorm((Qmax-tryQ)/trySE)+(tryQ-Qmax)*pnorm((tryQ-Qmax)/trySE) )## 7.5 p. 121              
  }
}

"calc.lrthreshold" <- function(object, ...) UseMethod("calc.lrthreshold") ## makes it easy to develop new inference methods

# has an SLik method

calc.lrthreshold.default <- function(object,dlr=NULL,verbose=interactive(),...) {
  stop("No default method for calc.lrthreshold. Provide a method.")
}

.rparam_from_SLik <- function(object, logLname, fittedPars, level, n, useEI, EIsampleFactor, useCI){
  if (inherits(object,"SLikp")) {
    surfaceData <- object$tailp
  } else {
    surfaceData <- object$logLs
  } 
  ## replace pairs with low predicted lik by pairs with high predicted lik
  pred <- predict(object$fit,binding=logLname) 
  ## : corrected 11/07/2016: pred was predict(object,.) using object$logLs, not object$fit, thereby not matching uli
  uli <- .ULI(object$fit$data[,fittedPars])
  table_uli <- table(uli)
  CIthreshold <- object$MSL$maxlogL - (qchisq(level,df=1)/2)
  singletsLvls <- names(table_uli)[table_uli==1L]
  singletsBool <- uli %in% singletsLvls
  single_pred <- pred[singletsBool,] ## (non-redundant) table of predicted values
  sort_single_pred <- single_pred[order(single_pred[,logLname],decreasing=TRUE),]
  doubletsBool <- ! uli %in% singletsLvls
  double_pred <- pred[doubletsBool,]
  unique_pred <- unique(double_pred)
  sort_unique_pred <- unique_pred[order(unique_pred[,logLname],decreasing=FALSE),]
  ncomp <- min(nrow(sort_unique_pred),nrow(sort_single_pred))
  if (ncomp>0L) {
    n_sub <- max(0,which(sort_single_pred[1:ncomp,logLname] > sort_unique_pred[1:ncomp,logLname]))
  } else n_sub <- 0L ## number than can be substituted
  if (n_sub>0L) {
    n_top <- length(which(sort_single_pred[,logLname]> CIthreshold))
    if (n_top< n_sub) {
      goodrows <- seq(n_sub) ## takes all highest
    } else goodrows <- sample(n_top,n_sub) ## random sampling of top  
    goodrows <- rownames(sort_single_pred[goodrows,,drop=FALSE])
    newpairs <- surfaceData[goodrows,object$colTypes$allPars,drop=FALSE]
    ## REMOVE ONE REPLICATE OF EACH "POOR" PAIR IN surfaceData
    surfaceData <- surfaceData[ ! rownames(surfaceData) %in% rownames(sort_unique_pred)[seq(n_sub)],]
  } else newpairs <- NULL
  # and the really slow part: 
  trypoints <- do.call(Infusion.getOption("rparamfn"),
                       list(object=object,n=n,useEI=useEI,tryn=EIsampleFactor*n,useCI=useCI,level=level,verbose=FALSE))
  trypoints <- rbind(trypoints,newpairs)
  return(list(trypoints=trypoints, surfaceData=surfaceData))
}

.update_seed <- function(object) {
  cl_seed <- attr(object$logLs,"workflow_env")$cl_seed
  if ( ! is.null(cl_seed)) {
    if (FALSE) {
      # if I want the cl_seed to control entirely the parallel operations, I must use it (rather that a global RNG)
      # to determine the next value of the parallel RNG: 
      R.seed <- get(".Random.seed", envir = .GlobalEnv) # save parent RNG state
      ori <- RNGkind("L'Ecuyer-CMRG")
      set.seed(cl_seed)  ## full control of parallel RNG by cl_seed
      cl_seed <- sample(.Machine$integer.max, 1) # using (and altering) global RNG to change parallel RNG
      RNGkind(ori[1])
      assign(".Random.seed", R.seed, envir = .GlobalEnv) # restore parent RNG state
      # BUT then: either I update the stored cl_seed in the SLik object => non repeatable code à a different cl_seed is used in input
      # or I don't => the same cl_seed will always be produced in output of .update_seed(), which looks suspect. Hence...
    } else { #...  I should used the global RNG to control reproducibly the .update_seed() result 
      cl_seed <- sample(.Machine$integer.max, 1) ## using (and altering) the global RNG
    }
  }
  cl_seed
}

"refine" <- function(object, ...) UseMethod("refine") ## makes it easy to develop new inference methods

# has SLikp and SLik method

refine_nbCluster <- function(nr, onlymax=7L) {
  res <- seq_nbCluster(nr=nr)
  maxres <- max(res)
  if (maxres>onlymax) res <- maxres
  res
}

`refine.default` <- local({
  update_warned <- FALSE
  ## si Simulate est exterieure, il faut que l'utilisateur puisse decomposer la fn et sample_volume doit être public...
  function(object, ## SLikp or SLik
           surfaceData, ## object$logLs or object$tailp, with stat.obs attribute, etc
           Simulate=attr(surfaceData,"Simulate"),
           maxit=1,n=NULL,useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE),newsimuls=NULL, trypoints=NULL,
           CIs=useCI, useCI=TRUE, level=0.95,
           verbose=list(most=interactive(),final=NULL,movie=FALSE,proj=FALSE),
           precision = Infusion.getOption("precision"),
           nb_cores=NULL, packages=attr(object$logLs,"packages"), env=attr(object$logLs,"env"),
           method, ## "GCV" and HLfit methods for Slik objects; mixmodCluster or... for SLik_j objects but the latter better controlled by "using"
           using=object$using, ## 
           eval_RMSEs=TRUE,
           update_projectors=FALSE,
           cluster_args=list(),
           cl_seed=.update_seed(object),
           nbCluster=quote(refine_nbCluster(nr=nrow(data))),
           ...) {
    fittedPars <- object$colTypes$fittedPars
    #
    if (!is.list(verbose)) verbose <- as.list(verbose)
    if (is.null(names(verbose))) names(verbose) <- c("most","movie","final","proj")[seq_len(length(verbose))]
    if (is.null(verbose$most)) verbose$most <- interactive()
    if (is.null(verbose$final)) verbose$final <- (interactive() && length(fittedPars)<3L)
    if (is.null(verbose$movie)) verbose$movie <- FALSE
    if (is.null(verbose$proj)) verbose$proj <- FALSE
    #
    if (is.null(using)) using <- Infusion.getOption("mixturing")
    #
    if (is.null(packages)) packages <- packages$add_simulation
    reftable_cluster_args <- .lookup(cluster_args,try_in="reftable")
    reftable_cluster_args[["RMSE"]] <- NULL
    RMSE_cluster_args <- .lookup(cluster_args,try_in="RMSE")
    RMSE_cluster_args[["reftable"]] <- NULL
    
    if (is.null(reftable_cluster_args$spec)) reftable_cluster_args$spec <- nb_cores # which means that cluster_args$spec overrides nb_cores
    if (is.null(RMSE_cluster_args$spec)) RMSE_cluster_args$spec <- 1L # which means that cluster_args$spec overrides nb_cores
    it <- 0L
    previous_cumul_iter <- max(object$logLs$cumul_iter)
    RMSEs <- get_from(object,"RMSEs") # LR_RMSE+ logLik() MSE ## 
    if (is.null(RMSEs)) RMSEs <- 1e10
    stat.obs <- attr(surfaceData,"stat.obs")
    if ( ! is.null(newsimuls) && maxit>1) stop("'maxit'>1 is incompatible with user-provided 'newsimuls'") 
    if  ( target_reached <- ( length(RMSEs) && 
                              ( ! (anyNA_RMSE <- anyNA(RMSEs))) && # but currently (with the reftable method at least) there is no NA in the RMSE table
                              all(RMSEs<precision))) {
      cat("Target precision appears to be already reached in input object.\n") ## nevertheless continue for one iteration 
    }
    EIsampleFactor <- 30
    while( it==0L ##always perform one iteration on request  
           || (it <maxit && ( ! target_reached)) ) {
      if (anyNA_RMSE) { # does not occur, at least in routine reftable usage
        message(paste0("CI bounds for some parameters not available: iterations will continue."))
      }
      if (maxit>1L && verbose$most) cat(crayon::yellow(paste("iter = ",it+1L," (total=",previous_cumul_iter+it+1L,"):\n",sep="")))
      logLname <- object$colTypes$logLname
      ## (1) Provide newsimuls
      if ( is.null(newsimuls)) {
        ## (1.1) generate parameter points
        if (inherits(object,"SLik_j")) {
          if (is.null(trypoints)) {
            freq_good <- list()
            samplingType <- Infusion.getOption("samplingType")
            locblob <- .rparam_from_SLik_j(object,frac=samplingType["default"],fittedPars=fittedPars,level=level, target_size=n) 
            trypoints <- locblob$trypoints
            attr(trypoints,"freq_good") <- locblob$freq_good
            freq_good$default <- locblob$freq_good
            if ((frac <- samplingType["posterior"])>0) {
              if (inherits(object$jointdens,"Mclust")) {
                object$postdens <- .conditional_mclust(object$jointdens,fittedPars=fittedPars,given=stat.obs)
              } else object$postdens <- .conditional_Rmixmod(object$jointdens,#fittedPars=fittedPars,
                                                             given=stat.obs)
              locblob <- .rparam_from_SLik_j_post(object,frac=samplingType["posterior"]) 
              freq_good$posterior <- locblob$freq_good
              trypoints <- rbind(trypoints,locblob$trypoints)
            } 
          } else freq_good <- attr(object$logLs,"freq_good") # keep older value
          surfaceData <- object$logLs
        } else {
          if (is.null(trypoints)) {
            if (verbose$most) cat("\nDesigning new parameter points (may be slow)...\n")
            locblob <- .rparam_from_SLik(object, logLname, fittedPars, level, n, useEI, EIsampleFactor, useCI=(CIs && useCI))
            trypoints <- locblob$trypoints
          }
          surfaceData <- locblob$surfaceData
        }
        if (is.null(Simulate)) {
          if (verbose) message("No 'Simulate' function specified: only parameter points are returned.")
          return(trypoints) #      RETURN                 but it is difficult for the user to reproduce "manually" the call to add_...() below; 
        }
        ## (1.2) Simulate for given parameter points
        if (inherits(object,"SLik_j")) {
          newsimuls <- add_reftable(Simulate=Simulate,par.grid=trypoints,verbose=verbose$most,
                                    control.Simulate=attr(surfaceData,"control.Simulate"),
                                    cluster_args=reftable_cluster_args, packages=packages$add_simulation, env=env,
                                    # Simulate_input=attr(surfaceData,"Simulate_input"), 
                                    cl_seed=cl_seed)     
        } else {
          newsimuls <- add_simulation(Simulate=Simulate,par.grid=trypoints,verbose=verbose$most,
                                      control.Simulate=attr(surfaceData,"control.Simulate"),
                                      cluster_args=reftable_cluster_args, packages=packages$add_simulation, env=env,
                                      # Simulate_input=attr(surfaceData,"Simulate_input"), 
                                      cl_seed=cl_seed)   
        }
        trypoints <- NULL # otherwise they would be used in next iteration (cf is.null(trypoints))
      } else freq_good <- list(default=attr(newsimuls,"freq_good")) ## _FIXME_ quick patch to avoid bug
      ## (2) Regenerate surface object 
      if (inherits(object,"SLik_j")) {
        checkNA <- apply(newsimuls,1L,anyNA) # newsimuls of class "reftable" and "data.frame"
        newsimuls <- newsimuls[!checkNA,,drop=FALSE] ## FR->FR quick patch
        if ( ! (is.matrix(newsimuls) || is.data.frame(newsimuls)) ) {
          stop("'newsimuls' must be a matrix or data.frame for refine.Slik_j() method.")
        }
        if ( ! is.null(projectors <- object$projectors)) { # the following assumes it's an environment
          # F I X M E which attributes "should" newsimuls have ? LOWER is sufficient for .update_raw_data()) code, but
          # more generally add_reftable (maybe, though perhaps not always) expects also UPPER and Simulate 
          if (is.null(attr(newsimuls,"LOWER"))) attr(newsimuls,"LOWER") <- attr(object$raw_data,"LOWER") # project(newsimuls,...) needs it
          raw_data <- .update_raw_data( object$raw_data, newsimuls)
          project_methodArgs <- .lookup(cluster_args, try_in="project")
          if (update_projectors) {
            proj_names <- ls(object$projectors) # list of objects in environment
            if ( length(proj_names)) {
              for (st in proj_names) {
                projector <- get(st,envir = projectors)
                projcall <- attr(projector,"project_call")
                projcall$verbose <- verbose$proj
                projcall$data <- raw_data[,names(projcall$data)]
                methodArgs <- projcall$methodArgs
                methodArgs[names(project_methodArgs)] <- project_methodArgs
                projcall$methodArgs <- methodArgs
                if ( (! update_warned) && ! inherits(projector, c("HLfit", "randomForest","ranger"))) {
                  warning(paste0(c("You may be updating a projection using a method\n",
                                   "for which Infusion has no built-in method to avoid overfitting (see help('project'))."))) 
                  update_warned <<- TRUE
                }
                projector <- eval(projcall) #,envir=environment(projcall))
                assign(st, value=projector, envir=projectors)
              }
            } else warning("non-NULL $projectors but ls() has zero length: check code or object")
            jointEDF <- project(raw_data,projectors=projectors, is_trainset=TRUE) ## AFAICS is_trainset is true here => passed to .predictWrap() to use oob predictions without check.
                                                                                  ## No benefit of project_methodArgs here. 
            jointEDF$cumul_iter <- c(object$logLs$cumul_iter, rep(previous_cumul_iter + it +1L, nrow(newsimuls)))
            stat.obs <- project(attr(stat.obs,"raw_data"),projectors=projectors)
            if (FALSE) { # only for diagnostic purposes...
              zzut <- (raw_data$mu-4)^2<0.01 & (raw_data$s2-1.25)^2<0.1
              plot(raw_data$var, jointEDF$VAR)
              abline(0,1)
              points(raw_data[zzut,]$var, jointEDF[zzut,]$VAR, col="red")
            }
          } else {
            newsimuls <- project(newsimuls,projectors=projectors, 
                                 methodArgs=project_methodArgs, # need to control parall independently for simuls and for project. 
                                 use_oob=update_projectors) # i.e. use_oob=FALSE: since projectors were not updated, 
                                                                                         # the new data are out of the training set => no oob predictions.
            newsimuls$cumul_iter <- previous_cumul_iter + it +1L
            jointEDF <- rbind(object$logLs,newsimuls)
          }
        } else {
          newsimuls$cumul_iter <- previous_cumul_iter + it +1L
          jointEDF <- rbind(object$logLs,newsimuls)
        }
        jointEDF <- structure(jointEDF, allPars=object$colTypes$allPars)
        object <- infer_SLik_joint(data = jointEDF,
                                   stat.obs = stat.obs,
                                   Simulate = attr(object$logLs,"Simulate"),
                                   using=using,
                                   nbCluster=nbCluster,
                                   verbose=verbose$most)
        object$latestPoints <- nrow(jointEDF)+1-seq_len(nrow(newsimuls)) ## for plots
        if ( ! is.null(projectors)) object$raw_data <- raw_data
      } else {
        if ( ! is.null(projectors <- object$projectors)) newsimuls <- project(newsimuls,projectors=projectors)
        if ( ! inherits(newsimuls,"list") ) {
          stop("'newsimuls' must be a _list_ of empirical distributions.")
        }
        if (inherits(object,"SLikp")) {
          mc <- attr(object$tailp,"call") 
        } else {
          mc <- attr(object$logLs,"call") # call to infer_logLs
        } 
        mc$object <- newsimuls
        mc$verbose <- verbose$most
        mc$`stat.obs` <- attr(object$logLs,"stat.obs") ## bc otherwise mc$`stat.obs` stores a promise such as 'Sobs'
        # if (inherits(object,"SLikp")) arglist$refDensity <- object$refDensity 
        if (length(unlist(cluster_args))) mc$cluster_args <- reftable_cluster_args ## do not erase the value in the initial call stored in the object
        if ( ! is.null(packages)) mc$packages <- packages$logL_method ## idem
        newlogLs <- eval(mc)
        if (verbose$most) cat ("\n")
        successrate <- length(which(newlogLs$isValid>0))/nrow(newlogLs)
        EIsampleFactor <- EIsampleFactor * 0.98/successrate
        newlogLs$cumul_iter <- previous_cumul_iter + it +1L
        surfaceData <- rbind(surfaceData,newlogLs)
        itmethod <- method[min(length(method),it+1L)] ## may be overriden below when hat(nu) is low.
        ## tests whether resmoothing can yield substantial improvements:
        if (object$fit$spaMM.version<"2.4.26") {
          corrPars1 <- object$fit$corrPars[["1"]]
        } else corrPars1 <- get_ranPars(object$fit,which="corrPars")[["1"]]
        allFix <- c(corrPars1[c("nu","rho")],list(lambda=object$fit$lambda[1],phi=object$fit$phi[1],beta=fixef(object$fit)))
        previousRho <- allFix$rho ## full length vector of scale params
        if (is.null(previousRho)) {
          stop("is.null(previousRho) in 'refine.default'. Check code (trRho?).")
        }
        smoothingOK <- (allFix$nu>3.95) ## provisional assessment
        if (smoothingOK) { ## further assessments
          # Test whether the old fit predicts well the new points
          prednewfromold <- predict(object,newdata=newlogLs[,fittedPars,drop=FALSE],
                                    variances=list(linPred=TRUE,disp=TRUE))  
          msepred <- mean((prednewfromold[,1]-newlogLs[,logLname])^2)
          respvar <- mean(attr(prednewfromold,"predVar")) + object$fit$phi ## FR->FR does not account forprior weights
          relerr <- respvar/(1e-6+msepred)
          #print(paste("relerr:",relerr))
          smoothingOK <- (relerr>0.8 && relerr<1.25) ## practically always, given distrib of estimator 'relerr'  
          #
          if (FALSE) {
            # assessment by perturbing rho
            currsurf <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=allFix)
            currp_bv <- currsurf$fit$APHLs$p_bv
            testCorr <- allFix
            testCorr$rho <- previousRho*1.1
            newp_bv <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=testCorr)$fit$APHLs$p_bv
            if ( smoothingOK <- (newp_bv < currp_bv) ) { ## no progress by increasing rho
              testCorr$rho <- previousRho/1.1
              newp_bv <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=testCorr)$fit$APHLs$p_bv
              smoothingOK <- (newp_bv < currp_bv) ## no progress by decreasing rho
            }
          }
        } else if (itmethod=="GCV") itmethod <- "REML" ## (overrides GCV user option) enforces REML when GCV does not estimate nu correctly
        if (smoothingOK) {
          #print(paste("smoothingOK:",smoothingOK))
          # object <- infer_surface(surfaceData,method="newdata",verbose=FALSE,allFix=allFix) 
          ## newdata => source des biais évidents: pas assez de pointsenhaut de la surface ?
          object <- infer_surface(surfaceData,method=itmethod,verbose=FALSE,allFix=allFix)
          ## : ie infer surface from purgedlogLs from new data,method=itmethod, with allFix=<previous pars>
        } else {
          object <- infer_surface(surfaceData,method=itmethod,verbose=verbose$most) ## new smoothing;could reuse corrpars 
          if (verbose$most) {
            if (object$fit$spaMM.version<"2.4.26") {
              corrPars1 <- object$fit$corrPars[["1"]]
            } else corrPars1 <- get_ranPars(object$fit,which="corrPars")[["1"]]
            vranfix <- unlist(corrPars1[c("nu","rho")])
            cat(paste(paste(paste(names(vranfix),"=",signif(vranfix)),collapse=", ")," (estimated by ",itmethod,")\n",sep=""))
          }
        }
        object$latestPoints <- nrow(surfaceData)+1-seq_len(nrow(newlogLs)) ## for plots
      }    
      # (3) maximization, new CIs, new MSEs
      # ! ! ! ! ! above calls to infer_surface() potentially change rownames of object$logLs relative to those of surfaceDate
      # bc infer_surface() -> make.names(). One solution is to get surfaceData from the new MSL object: 
      object <- MSL(object, CIs=CIs, level=level, verbose=verbose$most, eval_RMSEs=eval_RMSEs, cluster_args=RMSE_cluster_args)
      if (verbose$movie) {
        plot(object, from_refine=TRUE, ...) # from which RMSEs may be NULLified 
        if (! is.null(object$MSL$init_from_prof)) {
          cat("Better likelihood found. Maximizing again...\n")
          object <- MSL(object, CIs=CIs, level=level, verbose=verbose$most, eval_RMSEs=eval_RMSEs, cluster_args=RMSE_cluster_args)
        } 
      }
      if (inherits(object,"SLik_j")) {attr(object$logLs,"freq_good") <- freq_good} 
      RMSEs <- get_from(object,"RMSEs") # LR_RMSE+ logLik() MSE ## 
      target_reached <- ( length(RMSEs) && 
                            ( ! (anyNA_RMSE <- anyNA(RMSEs))) && # but currently (with the reftable method at least) there is no NA in the RMSE table
                            all(RMSEs<precision))
      it <- it+1
      newsimuls <- NULL ## essential for the test
    } ## end while() loop
    
    if (verbose$most && it <= maxit && target_reached) cat("\nIterations terminated because all(RMSEs<precision) where precision =",precision,"\n")

    if (( ! verbose$movie) && verbose$final) {
      plot(object, from_refine=TRUE, ...) # from which RMSEs may be NULLified 
      if (! is.null(object$MSL$init_from_prof)) {
        cat("Better likelihood found. Maximizing again...\n")
        object <- MSL(object, CIs=CIs, level=level, verbose=verbose$most, eval_RMSEs=eval_RMSEs, cluster_args=cluster_args)
        if (inherits(object,"SLik_j")) {attr(object$logLs,"freq_good") <- freq_good} # long missing 
      } 
    }
    if (verbose$most) { # as meant by the doc.) 
      if  ( any(is.na(RMSEs))) {
        cat("Target precision could not be ascertained.\n")
      } else if (is.null(get_from(object,"RMSEs"))) {
        cat("Precision has not been evaluated. Use argument 'eval_RMSEs=TRUE' to force its computation.\n")
      } else if  ( any(RMSEs>precision) ) {cat("\nTarget precision does not appear to be reached.\n")}
    }
    return(object)
  }
})
