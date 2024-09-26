
.calc_EI <- function(summInferObject,points,Qmax=NULL) {  ## for both SLik and SLikp
  trypred <- predict(summInferObject,newdata=points,variances=list(predVar=TRUE))
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
  trypoints <- do.call(Infusion.getOption("rparamFn_SLik"),
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
      RNGkind(ori[1]) # or do.call("RNGkind", as.list(ori))
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



reproject <- function(object, eval_RMSEs=NULL, CIs=NULL, ...) {
  if (is.null(eval_RMSEs)) eval_RMSEs <- length(object$par_RMSEs$par_RMSEs)>0L 
  if (is.null(CIs)) CIs <- length(object$CIobject$CIs)>0L
  # : make sure they are logical: numbers have distinct interpretations.
  refine(object = object,..., ntot=0L, maxit=1L, eval_RMSEs=eval_RMSEs, CIs=CIs, update_projectors=TRUE)
}

recluster <- function(object, eval_RMSEs=NULL, CIs=NULL, 
                      update_projectors=FALSE, ...) {
  if (is.null(eval_RMSEs)) eval_RMSEs <- length(object$par_RMSEs$par_RMSEs)>0L 
  if (is.null(CIs)) CIs <- length(object$CIobject$CIs)>0L
  # : make sure they are logical: numbers have distinct interpretations.
  refine(object = object,..., ntot=0L, maxit=1L, eval_RMSEs=eval_RMSEs, CIs=CIs, 
         update_projectors=update_projectors)
}

## There is a documented reproject()
# .reforest <- function(object, eval_RMSEs=NULL, CIs=NULL, ...) {
#   refine(object = object,..., ntot=0L, maxit=1L, eval_RMSEs=FALSE, CIs=FALSE, update_projectors=TRUE)
# }

# Next line does not make sense
# .update_projector <- function(object, ...) UseMethod("update_projector") # internal generic... not used...

.update_projector.SLik_j <- function(object, parm=NULL, proj,  # arguments work as in plot_proj 
                                     verbose=FALSE, 
                                     reftable_raw,
                                     methodArgs=list(), 
                                     ...) {
  if (is.null(parm)) {
    projector <- object$projectors[[proj]] # there must be proj
    # parm <- paste(projector$call$formula[[2]])
    parm <- attr(projector,"project_call")$x
  } else {
    projnames <- ls(object$projectors)
    # pars <- sapply(projnames, function(st) paste(object$projectors[[st]]$call$formula[[2]])) # may be removed...
    pars <- sapply(projnames, function(st) attr(object$projectors[[st]],"project_call")$x) # ___F I X M E___ formula not always present
    proj <- names(which(pars==parm))
    projector <- object$projectors[[proj]]
  }
  .update_projector.default(projector, verbose=verbose, reftable_raw=reftable_raw, 
                           methodArgs=methodArgs, ...)
}

.update_projector.default <- function(projector, 
           verbose=FALSE, reftable_raw, 
           methodArgs=list(), # optional modif of the methodArgs args stored in the project_call:
             # they are for project.character(), not directly ranger() 
             # => project's default splitrule = "extratrees" can be reversed by setting it = "variance"
           ...) {
    projcall <- attr(projector,"project_call") # will call project.character
    projcall$verbose <- verbose
    projcall$data <- reftable_raw[,names(projcall$data)]
    projcall$methodArgs[names(methodArgs)] <- methodArgs
    if ( (! .one_time_warnings$update_warned) && ! inherits(projector, c("HLfit", "randomForest","ranger"))) {
      warning(paste0(c("You may be updating a projection using a method\n",
                       "for which Infusion has no built-in method to avoid overfitting (see help('project'))."))) 
      .one_time_warnings$update_warned <- TRUE
    }
    projector <- eval(projcall) #,envir=environment(projcall))    
    projector
  }

deforest_projectors <- function(object) {
  projs <- object$projectors 
  for (proj_it in names(projs)) {
    projs[[proj_it]]$forest <- NULL
    projs[[proj_it]]$call <- NULL
  }
  "*Input* object has been internally modified." # 'object$projectors' environment modified
}

.missing_forests <- function(object) {
  projectors <- object$projectors
  if (is.null(projectors)) return(FALSE)
  projnames <- ls(projectors)
  missing_forests <- logical(length(projnames))
  names(missing_forests) <- projnames
  for (st in projnames) missing_forests[st] <- (inherits(projectors[[st]], "ranger") && 
                                                 is.null(projectors[[st]]$forest$split.varIDs))
  missing_forests
}


# Currently this is documentation but it makes easy to reimplement this old approach.
.process_target_LR_v2.1.85 <- function(target_LR, level, maxit, target_cumns, ...) {
  ref <- qchisq(level,df=2)/2
  template_target_LR <- c(rep(ref, maxit-1L),2*ref)
  names(template_target_LR) <- target_cumns

  if (is.null(target_LR)) {
    target_LR <- template_target_LR
  } else {
    if (is.null(names_target_LR <- names(target_LR))) { # interpret nameless 'target_LR' argument
      len_tLR <- length(target_LR)
      if (len_tLR==1L) {
        target_LR <- rep(target_LR, maxit-1L)
        names(target_LR) <- target_cumns[-maxit]
      } else if (len_tLR==2L) {
        target_LR <- c(rep(target_LR[1], maxit-1L),target_LR[2])
      } else if (len_tLR>maxit) target_LR <- target_LR[seq_len(maxit)]
      names(target_LR) <- target_cumns[seq_along(target_LR)]
    } else if (length(setdiff(names_target_LR,as.character(target_cumns)))) {
      stop("'target_LR' names do not match cumulative reftable sizes.")
    }
    template_target_LR[names(target_LR)] <- target_LR
  }
  target_LR <- c(template_target_LR, c("Inf"=2*ref)) # maxit +1L terms
}

.process_target_LR <- function(target_LR, # modifier of the default  template_target_LR
                               ref, 
                               level,
                               maxit, target_cumns, workflow_design, prev_size, ...) {
  npar <- workflow_design$npar
  subblock_size <-  workflow_design$subblock_size
  refsizes <- max(prev_size, subblock_size)+target_cumns # >= subblock_size
  template_target_LR <- ref/log(1+refsizes/subblock_size) # <= ref/log(2)
  min_value <- eval(.Infusion.data$options$target_LR_min_value) # uses 'level' (and is ad hoc)
  template_target_LR <- pmax(template_target_LR, min_value) # >= min_value
  names(template_target_LR) <- target_cumns
  
  if (is.null(target_LR)) {
    target_LR <- template_target_LR
  } else { # rethink later ___F I X M E___ well this seems OK
    if (is.null(names_target_LR <- names(target_LR))) { # interpret nameless 'target_LR' argument
      len_tLR <- length(target_LR)
      if (len_tLR==1L) {
        target_LR <- rep(target_LR, maxit)
        names(target_LR) <- target_cumns[-maxit]
      } else if (len_tLR==2L) {
        target_LR <- c(rep(target_LR[1], maxit-1L),target_LR[2])
      } else if (len_tLR>maxit) target_LR <- target_LR[seq_len(maxit)]
      names(target_LR) <- target_cumns[seq_along(target_LR)]
    } else if (length(setdiff(names_target_LR,as.character(target_cumns)))) {
      stop("'target_LR' names do not match cumulative reftable sizes.")
    }
    template_target_LR[names(target_LR)] <- target_LR
  }
  target_LR <- c(template_target_LR, c("Inf"=min_value)) # maxit +1L term, all >= min_value
  target_LR # *named*; controls .rparam_SLik_j_in_out() and .get_subrows_for_proj()
}

.get_subrows_for_proj <- function(object, init_reft_size, reftable_raw, target_LR) {
  target_minsize <- .Infusion.data$options$target_size_factor*length(object$colTypes$fittedPars)
  nr <- nrow(reftable_raw)
  if (nr<= target_minsize) {
    return(seq_len(nr))
  } else {
    pred_logls <- predict(object, newdata=reftable_raw, which="safe")
    max_logl <- max(pred_logls)
    # initial reftable + top:
    subrows <- unique(c(seq(init_reft_size),
                        which(pred_logls>max_logl-target_LR)))
    if (length(subrows)<target_minsize) { # if target_minsize not reached by above criterion,
      # we bypass the target_LR criterion to sample a large enough "top"
      actual_topsize <- min(target_minsize, length(pred_logls))-init_reft_size 
      order_logls <- order(pred_logls, decreasing = TRUE)
      subrows <- unique(c(seq(init_reft_size),
                          order_logls[seq_len(actual_topsize)]))
    }
    subrows
  }
}


.process_update_projectors <- function(update_projectors, maxit, target_cumns, 
                                       workflow_design, prev_size, ...) {
  if (is.logical(update_projectors)) {
    arg_len <- length(update_projectors)
    if ((miss_its <- maxit-arg_len)>0L) 
      update_projectors <- update_projectors[c(seq_len(arg_len), rep(arg_len, miss_its))] # logical vector
    update_projectors <- prev_size + target_cumns[update_projectors] # reftable sizes
  } else if (is.numeric(update_projectors)) {
    # assuming that user gave reftable size(s) at which projectors should be evaluated.
    if (maxit > 1L) {
      update_projectors <- update_projectors[update_projectors>prev_size]
    } else update_projectors <- update_projectors[update_projectors>=prev_size]
  } else { # NULL update_projectors (the default): 
    ## sets all refsizes as updating thresholds (as TRUE implies) 
    ##  BUT the kept 'user_update_projectors' will allow further control in refine.default().
    ## Also a fix for failures to reach a prescribed prev_size:
    minimal_prev_size_for_reprojection <- max(prev_size, workflow_design$subblock_size) 
    update_projectors <- minimal_prev_size_for_reprojection+target_cumns 
  }
  update_projectors <- c(update_projectors, Inf) # reftable sizes
  update_projectors
}

.reformat_refine_verbose <- function(verbose, fittedPars) {
  if (!is.list(verbose)) verbose <- as.list(verbose)
  if (is.null(names(verbose))) names(verbose) <- (names(formals(refine.default)$verbose)[-1])[seq_len(length(verbose))]
  if (is.null(verbose$most)) verbose$most <- (interactive() || .is_devel_session())
  if (is.null(verbose$notable)) verbose$notable <- TRUE
  if (is.null(verbose$final)) verbose$final <- (interactive() && length(fittedPars)<3L)
  if (is.null(verbose$movie)) verbose$movie <- FALSE
  if (is.null(verbose$proj)) verbose$proj <- FALSE
  if (is.null(verbose$MAF)) verbose$MAF <- verbose$most # MAF=2L => most verbose
  if (is.null(verbose$rparam)) verbose$rparam <- .is_devel_session()
  if (is.null(verbose$progress_bars)) verbose$progress_bars <- verbose$most
  if (is.null(verbose$debug)) verbose$debug <- FALSE # to debug try(plot...)
  verbose$from_refine <- TRUE
  verbose
}

# 
.update_timeStamps <- function(timeStamps=list(), 
                               previous, # environment, if the previous object is updated (projectors...)
                                         # or NULL if the previous object is not to be updated (proj_stats are not) 
                               ...) {
  if (is.null(timeStamps)) timeStamps <- as.list(previous) # all NULLs become list()
  dotlist <- list(...)
  if (is.environment(previous)) 
    for (st in names(dotlist)) assign(st,dotlist[[st]], envir=previous) # modified 'previous' envir
  timeStamps[names(dotlist)] <- dotlist
  timeStamps # returned _list_
}

.wrap_update_projectors <- function(projectors, sub_reftable, 
                                    verbose, project_methodArgs) {
  proj_names <- ls(projectors) # list of objects in environment
  np <- length(proj_names)
  if (np) {
    if (verbose$most) {
      if (.is_devel_session()) {
        devel_mess <- cli::col_green(paste(nrow(sub_reftable), "samples used.\n"))
      } else devel_mess <- "\n"
      cat("Updating projectors... ", devel_mess)
    } 
    if (verbose$proj) {
      plot.new()
      if (np>3L) {
        mai <- c(0.6,0.6,0.1,0.1) # b l u r
      } else mai <- NULL
      opar <- .fittedPars_plot_pars(np=np, 
                                    plotpar=list(pch=20L), 
                                    ylab=NULL,
                                    mai=mai) 
    }
    for (st in proj_names) {
      projector <- projectors[[st]]
      if ( ! is.null(projector)) {
        methodArgs <- attr(projector,"project_call")$methodArgs
        methodArgs[names(project_methodArgs)] <- project_methodArgs # the two user-level args override projcall$methodArgs
        projectors[[st]] <- .update_projector.default( # INPUT OBJECT'S 'projectors' ENVIRONMENT IS MODIFIED, but 
          # ... only the return object's logLs are updated. So, while plot_proj(upsliks[["14K"]], parm="log.Na.") looks like 
          # ... it plots the updated projections, it actually uses the not-updated $logLs... Really confusing. =>
          # ... a system of timeStamp's has been implemented to track and warn abotu such issues. 
          projector, verbose=verbose$proj, reftable_raw=sub_reftable, 
          methodArgs=methodArgs) 
      }
    }
    if (verbose$proj) par(opar)
  } else warning("non-NULL $projectors but ls() has zero length: check code or object")
}


.plot_cloud <- function(object, latestPoints, rparam_blob, verbose, fittedPars) {
  cloud_parm <- verbose$cloud_parm
  if ( ! is.null(rparam_blob)) {
    newpts <- object$logLs[latestPoints,,drop=FALSE]
    fill_info <- rparam_blob$fill_info
    locpts <- newpts[1:fill_info$NROW_selected,]
    plot(locpts[,cloud_parm],
         predict(object, newdata=locpts[,fittedPars, drop=FALSE]),
         col="red")
    if (fill_info$fill) {
      locpts <- newpts[fill_info$NROW_selected +(1:fill_info$fill),]
      points(locpts[,cloud_parm],
             predict(object, newdata=locpts[,fittedPars, drop=FALSE]),
             col="blue")
    }
    graphics::legend(x='bottom', legend=c("pardens","I.postdens"), fill =c("red","blue"))
  }
  plot(object$logLs[,cloud_parm],
       predict(object, newdata=object$logLs[,fittedPars]),
       #ylim=slik$MSL$maxlogL+c(-10,1), 
       #xlim=c(0.15,0.2), 
       pch=20, cex=0.5, xlab=cloud_parm, ylab="logL",
       col=object$logLs$cumul_iter)
  axis(3)
}

.get_initParam <- function(object) {
  if ( inherits(object$jointdens,"MAF")) {
    NULL 
  } else object$jointdens@parameters
}
  
## si Simulate est exterieure, il faut que l'utilisateur puisse decomposer la fn et sample_volume doit être public...
refine.default <- function(
  object,  ## SLik_j, SLikp or SLik
  ##       reference table simulations  
  Simulate=attr(surfaceData,"Simulate"),
  control.Simulate=attr(surfaceData,"control.Simulate"),
  newsimuls=NULL,
  ##       CIs
  CIs=workflow_design$reftable_sizes[useCI], 
  useCI=prod(dim(object$logLs))<12000L, level=0.95,
  ##       workflow design
  workflow_design=get_workflow_design(npar=length(fittedPars),
                                      n_proj_stats = length(statNames),
                                      n_latent=length(latentVars)), 
  maxit, ntot= maxit*.get_size_first_iter(object), n=NULL,
  ##       termination conditions 
  precision = Infusion.getOption("precision"),
  eval_RMSEs=workflow_design$reftable_sizes, 
  ##       verbosity
  verbose=list(notable=TRUE, most=interactive(),final=NULL,
               movie=FALSE,proj=FALSE,rparam=NULL, progress_bars=interactive()),
  ##       projection controls
  update_projectors=NULL,
  methodArgs=list(),
  ##       Likelihood surface modeling (up-to-date workflow)
  using=object$using, ## mclust...
  nbCluster=quote(refine_nbCluster(nr=nrow(data))),
  ##       parallelisation
  cluster_args=list(), nb_cores=NULL, env=get_from(object,"env"), 
  packages=get_from(object,"packages"), cl_seed=.update_seed(object),
  ##       obscure stuff
  target_LR=NULL,  
  ##       not explicitly needed in up-to-date workflow
  trypoints=NULL,
  surfaceData, ## object$logLs or object$tailp, with stat.obs attribute, etc
  method, ## "GCV" and HLfit methods for Slik objects; mixmodCluster or... for SLik_j objects but the latter better controlled by "using"
  useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE),
  rparamFn = Infusion.getOption("rparamFn"),
  ## 
  ... # map.asp may be passed through these dots
) {
  if (is.null(object)) stop("Input 'object' is NULL: check refine() call.")
  timeStamps <- as.list(object$timeStamps)
  colTypes <- object$colTypes
  # The ree next ones need to be in the closure envir when workflow_design is evaluated.
  fittedPars <- colTypes$fittedPars
  statNames <- colTypes$statNames 
  latentVars <- colTypes$latentVars 
  #
    verbose <- .reformat_refine_verbose(verbose, fittedPars)
    #
    if (is.null(using)) using <- Infusion.getOption("mixturing") # but in up to date objects, using=object$using is not NULL.
    #
    if (is.null(packages)) packages <- packages$add_simulation
    reftable_cluster_args <- .lookup(cluster_args,try_in="reftable")
    #reftable_cluster_args[["RMSE"]] <- NULL
    RMSE_cluster_args <- .lookup(cluster_args,try_in="RMSE")
    #RMSE_cluster_args[["reftable"]] <- NULL
    rparam_cluster_args <- .lookup(cluster_args,try_in="rparam")
    
    if (is.null(reftable_cluster_args$spec)) reftable_cluster_args$spec <- nb_cores # which means that cluster_args$spec overrides nb_cores
    if (is.null(RMSE_cluster_args$spec)) RMSE_cluster_args$spec <- 1L # which means that cluster_args$spec overrides nb_cores
    previous_cumul_iter <- max(object$logLs$cumul_iter)
    RMSEs <- get_from(object,"RMSEs") # LR_RMSE+ logLik() MSE ## 
    if (is.null(RMSEs)) RMSEs <- 1e10
    stat.obs <- attr(surfaceData,"stat.obs")
    if  ( target_reached <- ( length(RMSEs) && 
                              ( ! (anyNA_RMSE <- anyNA(RMSEs))) && # but currently (with the reftable method at least) there is no NA in the RMSE table
                              all(RMSEs<precision))) {
      cat("Target precision appears to be already reached in input object.\n") ## nevertheless continue for one iteration 
    }
    previous_reftable_size <- nrow(object$logLs)
    if ( ! is.null(newsimuls)) {
      ntot <- 0L # ntot can also be set to 0L by user to recompute the projections... (as doc'ed)
      if ( (! missing(maxit)) && maxit!=1L) message("'maxit'!=1L is ignored with user-provided 'newsimuls'.") 
      maxit <- 1L
      cumn_over_maxit <- FALSE
      newrawstats <- newsimuls # not NULL here
    } else { 
      cumn_over_maxit <- workflow_design$cumn_over_maxit 
      ## Get ntot:
      if (previous_cumul_iter==1L) {
        if (missing(ntot) || is.null(ntot)) {
          if (is.null(n)) {
            is_default_wf <- missing(workflow_design)
            ntot <- workflow_design$first_refine_ntot
            if (is.null(ntot)) {
              ntot <- 2L*.get_size_first_iter(object)
            } else if (verbose$most && is_default_wf) {
              message(paste0("Using default workflow_design: 1st-refine 'ntot'= ", ntot,";"))
              info <- workflow_design[c("refine_blocksize","final_reft_size")]
              names(info) <- c("default for next-refines 'ntot'", "suggested minimal final size")
              message(paste0(names(info),"=",info, collapse="; "),".")
            }
          } else ntot <- 2L*n 
        }
      } else {
        if (missing(ntot) || is.null(ntot)) {
          if (is.null(n)) {
            ntots <- workflow_design$reftable_sizes - previous_reftable_size
            ntots <- ntots[ntots>0L]
            if (length(ntots)) {
              ntot <- min(ntots)
            } else ntot <- workflow_design$refine_blocksize
            # ntot <- workflow_design$refine_blocksize
            if (verbose$most) message(paste("refine() with missing 'ntot': 'ntot' set to",ntot,
                                            "(as controlled by 'workflow_design')."))
          } else ntot <- n
        }
      }
      ## Get maxit:
      if (missing(maxit) || is.null(maxit)) {
        if (previous_cumul_iter==1L) {
          maxit <- max(1L, as.integer(log(ntot, 2)-log(.get_size_first_iter(object)/20, 2)))
          if (verbose$most) message(paste0("First refine() with missing 'maxit': 'maxit' set to ",maxit,"."))
        } else {
          maxit <- ceiling(ntot*workflow_design$subblock_nbr/workflow_design$refine_blocksize)
          if (verbose$most) message(paste0("refine() with missing 'maxit': 'maxit' set to ",maxit,"."))
        }
      } 
    } 
    
    if (previous_cumul_iter==1L) {
      if (maxit>2L) {
        seq_pow <- maxit-(1:(maxit))
        seq_pow <- seq_pow[seq_len(max(0,maxit-2L))] 
        target_nsims <- as.integer(ntot/2^(seq_pow))
        target_cumns <- c(cumsum(target_nsims[c(1,seq_along(target_nsims))]), ntot)
      } else if (maxit==2L) {
        target_cumns <- c(as.integer(ntot/2L), ntot)
      } else target_cumns <- ntot
    } else {
      target_nsims <- rep(as.integer(ntot/maxit), maxit-1L)
      target_cumns <- c(cumsum(target_nsims), ntot)
    }
    # cumn and target_cumns compare to ntot and all start from 0 within a refine
    
    npar <- workflow_design$npar
    target_LR_rparam <- .process_target_LR(
      target_LR=target_LR, 
      ref=qchisq(level,df=npar)/2, #  /2 seems useful...
      level=level,
      maxit=maxit, target_cumns=target_cumns,
      workflow_design=workflow_design,
      prev_size=previous_reftable_size)
    target_LR_proj <- .process_target_LR(
      target_LR=target_LR, 
      ref=qchisq(level,df=npar)/2, # without /2 from v2.1.109.5 to v2.1.112: does not seem good (B_13from17)
                                   # some replicates err (in the final iteration: estimates at bound, high lik...)
      level=level,
      maxit=maxit, target_cumns=target_cumns,
      workflow_design=workflow_design,
      prev_size=previous_reftable_size)
    target_LR_it <- 1L
    
    user_update_projectors <- update_projectors
    update_projectors <- .process_update_projectors(update_projectors=update_projectors, 
                                                    maxit=maxit, target_cumns=target_cumns,
                                                    workflow_design=workflow_design,
                                                    prev_size=previous_reftable_size)
    update_projectors_it <- 1L
    prev_sub_trainset <- object$proj_sub_trainset
    proj_trainset <- object$proj_trainset
    if (is.logical(CIs)) {
      arg_len <- length(CIs)
      if ((miss_its <- maxit-arg_len)>0L) CIs <- CIs[c(seq_len(arg_len), rep(arg_len, miss_its))]
      CIs <- c(target_cumns[CIs], Inf)
    } else {
      # assuming that user gave reftable size(s) at which CIs should be evaluated.
      CIs <- CIs-previous_reftable_size
      if (ntot > 0L) {
        CIs <- CIs[CIs>0L]
      } else CIs <- CIs[CIs>=0L]
      CIs <- c(CIs,Inf)
    } 
    eval_CIs_it <- 1L
    
    if (is.logical(eval_RMSEs)) {
      arg_len <- length(eval_RMSEs)
      if ((miss_its <- maxit-arg_len)>0L) eval_RMSEs <- eval_RMSEs[c(seq_len(arg_len), rep(arg_len, miss_its))]
      eval_RMSEs <- c(target_cumns[eval_RMSEs], Inf)
    } else {
      # assuming that user gave reftable size(s) at which RMSEs should be evaluated.
      eval_RMSEs <- eval_RMSEs-previous_reftable_size
      if (ntot > 0L) {
        eval_RMSEs <- eval_RMSEs[eval_RMSEs>0L]
      } else eval_RMSEs <- eval_RMSEs[eval_RMSEs>=0L]
      eval_RMSEs <- c(eval_RMSEs,Inf)
    } 
    eval_RMSEs_it <- 1L
    
    EIsampleFactor <- 30
    it <- 0L
    it_index <- 1L
    cumn <- 0L
    # d_maxlogL <- -Inf
    while( it==0L || ##always perform one iteration on request  
           ( # necessary conditions for continuing if it>1L:
             ( # "not enough points, eventually even if maxit has been reached"
               ( (! cumn_over_maxit) && it < maxit ) || # If maxit most important crit (old design), check maxit
               cumn<ntot # If cumn most important crit (new design), check cumn
             ) && 
             ! target_reached #ie target reached is a sufficient condition for terminating after it=1
           )
         ) {
      # previous_maxlogL <- object$MSL$maxlogL
      if (anyNA_RMSE) { # does not occur, at least in routine reftable usage
        message(paste0("CI bounds for some parameters not available: iterations will continue."))
      }
      if (maxit>1L) {
        if (verbose$most) cat(cli::col_yellow(paste("iter = ",it+1L," (total=",previous_cumul_iter+it+1L,"):\n",sep="")))
        target_nsim <- target_cumns[it_index]-cumn # do not set a new target when the previous one has not been reached.
      } else target_nsim <- ntot-cumn
      logLname <- object$colTypes$logLname
      ## (1) Provide newrawstats
      if (target_nsim>0L) { 
        ## (1.1) generate parameter points
        if (inherits(object,"SLik_j")) {
          if (is.null(trypoints)) {
            rparam_blob <- rparamFn( # default is .rparam_SLik_j_in_out(), alternative is .rparam_SLik_j_B_postdens()
              object,fittedPars=fittedPars,level=level, target_size=target_nsim, 
              target_LR=target_LR_rparam[target_LR_it], 
              verbose=verbose$rparam,
              # cluster_args=rparam_cluster_args,
              safe=.Infusion.data$options$safe_prof4rparam #  default=> .safe_optim (vs nlminb) in .ad_hoc_opt; 
                                                           # safe=FALSE is faster but was not helpful enough 
                                                           # nlminb call aims to use gradfn but does not handle(., which=safe) (?)
            ) 
            ## 
            trypoints <- rparam_blob$trypoints
            if (verbose$notable &&
                NROW(trypoints) < target_nsim) {
              locmess <- paste0("'target_size' was ",target_nsim," but only ", NROW(trypoints)," candidate points\n",
                                "retained after rejection step (and satisfying parameter constraints).")
            } 
            cumn <- cumn + nrow(trypoints)
            rparam_info <- rparam_blob$rparam_info
          } else rparam_info <- object$rparam_info # keep older value
          surfaceData <- object$logLs
        } else {
          if (is.null(trypoints)) {
            if (verbose$most) cat("\nDesigning new parameter points (may be slow)...\n")
            rparam_blob <- .rparam_from_SLik(object, logLname, fittedPars, level, target_nsim, useEI, EIsampleFactor, 
                                         useCI=( ! is.null(object$CIobject) && useCI))
            trypoints <- rparam_blob$trypoints
          }
          surfaceData <- rparam_blob$surfaceData
        }
        if (is.null(Simulate)) {
          if (verbose$most) message("No 'Simulate' function specified: only parameter points are returned.")
          return(trypoints) #      RETURN                 but it is difficult for the user to reproduce "manually" the call to add_...() below; 
        }
        ## (1.2) Simulate for given parameter points
        if (inherits(object,"SLik_j")) {
          newrawstats <- add_reftable(Simulate=Simulate, parsTable=trypoints,
                                      verbose=verbose$progress_bars, 
                                    control.Simulate=control.Simulate,
                                    cluster_args=reftable_cluster_args, packages=packages$add_simulation, env=env,
                                    # Simulate_input=attr(surfaceData,"Simulate_input"), 
                                    cl_seed=cl_seed)     
        } else {
          newrawstats <- add_simulation(Simulate=Simulate, parsTable=trypoints,verbose=verbose$most,
                                      control.Simulate=control.Simulate,
                                      cluster_args=reftable_cluster_args, packages=packages$add_simulation, env=env,
                                      # Simulate_input=attr(surfaceData,"Simulate_input"), 
                                      cl_seed=cl_seed)   
        }
        trypoints <- NULL # otherwise they would be used in next iteration (cf is.null(trypoints))
      } else {
        newrawstats <- newsimuls
        rparam_info <- object$rparam_info # keep older value
      }
      ## (2) Regenerate surface object 
      eval_CIs_here <- cumn>=CIs[eval_CIs_it]
      eval_RMSEs_here <- cumn>=eval_RMSEs[eval_RMSEs_it]
      if (inherits(object,"SLik_j")) {
        reftable_raw <- .get_reft_raw(object)
        NROW_newrawstats <- NROW(newrawstats) # newrawstats may be NULL in the case I regenerate projectors without new simuls
        if (NROW_newrawstats) {
          if ( ! (is.matrix(newrawstats) || is.data.frame(newrawstats)) ) {
            stop("'newsimuls' must be a matrix or data.frame for refine.Slik_j() method.")
          }
          newrawstats <- na.omit(newrawstats)
          NROW_newrawstats <- nrow(newrawstats)
          reftable_raw <- .update_reft_raw( reftable_raw, newrawstats) # without $cumul_iter
          attr(newrawstats,"LOWER") <- attr(reftable_raw,"LOWER") # (only) LOWER with latentVars needed for project 
        } 
        if ( ! is.null(projectors <- object$projectors)) { # the following assumes it's an environment
          # Which attributes "should" newrawstats have ? LOWER is sufficient for .update_reft_raw()) code, but
          # more generally add_reftable (maybe, though perhaps not always) expects also UPPER and Simulate 
          
          project_methodArgs <- .lookup(cluster_args, try_in="project")
          project_methodArgs[names(methodArgs)] <- methodArgs # user-level methodArgs override cluster_args, 
          
          ## eval update_projectors_here
          update_projectors_here <- nrow(reftable_raw)>=update_projectors[update_projectors_it] 
          if (update_projectors_here) { 
            # We'll need this in many cases:
            newest_sub_trainset <- .get_subrows_for_proj(object, init_reft_size=workflow_design$init_reft_size, 
                                             reftable_raw=reftable_raw[,fittedPars, drop=FALSE], 
                                             target_LR=target_LR_proj[target_LR_it])
            if (is.null(user_update_projectors) && 
                (subrows_overlap_thr <- .Infusion.data$options$subrows_overlap_thr)<1) {
              # The compared trainsets include the 1st reftable. 
              # Points are not necess ordered as in the reftable bc -> .get_subrows_for_proj -> (under some condition) ...order()
              subrows_overlap <- length(intersect(newest_sub_trainset, prev_sub_trainset))/length(newest_sub_trainset)
              update_projectors_here <- (subrows_overlap < subrows_overlap_thr) 
            } # else the explicit user control determines when projectors are updated.
            if (.is_devel_session()) {
              if ( ! update_projectors_here) {
                cat(cli::col_green("... Not updating projectors bc subrows not distinct enough ('subrows_overlap_thr' option).\n"))
              } 
            }
          }
          if (update_projectors_here) { # add safety check to avoid by default reprojection on huge sub_reftable 
            if (is.null(user_update_projectors)) 
              update_projectors_here <- length(newest_sub_trainset)<.Infusion.data$options$upd_proj_subrows_thr 
          }
          ## do not change update_projectors_here any further within current iteration!
          
          if ( ! update_projectors_here) { # second check for possible reforest
            missing_forests <- (it==0L & .missing_forests(object))
            any_missing_forest <- any(missing_forests)
            if (any_missing_forest) {
              if (ntot==0L) { # ntot=0 <=> reforest/recluster but within this block this is recluster
                newest_sub_trainset <- setdiff(object$proj_trainset, prev_sub_trainset)
                # using the alternative might make sense, except if both projectors and the densities are missing.
                # By using the old proj_trainset we avoid complicated code depending on 'using'...
              } else {
                newest_sub_trainset <- .get_subrows_for_proj(object, init_reft_size=workflow_design$init_reft_size, 
                                               reftable_raw=reftable_raw[,fittedPars, drop=FALSE], 
                                               target_LR=target_LR_proj[target_LR_it])
              }
            } 
          } 

          if (update_projectors_here || 
              (ntot>0L && any_missing_forest) # Not always for 'any_missing_forest' 
              # allows recluster  from saved 'deforested' object without reprojection.
              # But if ntot>0L this is not simply recluster
              ) {
            #### provide projectors
            # nth fix: avoid periodic behavior by using:
            proj_trainset <- unique(c(prev_sub_trainset, newest_sub_trainset))
            # So messages about subrows are not about all rows finally used.
            # Without it the whole iteration sometimes seems to alternate between two inferences of the 
            # liksurf (+associated sampling, projections...with newest_sub_trainset surprizingly variable over iterations). 
            sub_reftable <- reftable_raw[proj_trainset, ] # keeps attributes... 
            prev_sub_trainset <- newest_sub_trainset
            proj_trainset <- proj_trainset
            if (update_projectors_here) {
              .wrap_update_projectors(projectors, sub_reftable=sub_reftable, 
                                      verbose, project_methodArgs)
              timeStamps <- .update_timeStamps(timeStamps, previous=object$timeStamps, projectors=Sys.time())
              update_projectors_it <- update_projectors_it+1L
            } else if (any_missing_forest) {
              # either forest is missing (as per ranger(., write.forest=FALSE))
              # or I removed one of its bulky elements.
              if (verbose$notable) message("'forest' information not available from some projector(s). Regenerating it...")
              for (st in names(which(missing_forests))) {
                projectors[[st]] <- .update_projector.SLik_j(
                  object, proj=st, reftable_raw=sub_reftable, methodArgs = list(write.forest=TRUE))
              }
            } 
            
            #### update jointEDF object: (with $cumul_iter)
            ## fill will oob and non-oob predictions:
            jointEDF <- .project_reftable_raw(reftable_raw, projectors=projectors, 
                                is_trainset=FALSE, 
                                use_oob=FALSE)
            jointEDF[proj_trainset,] <- .project_reftable_raw(sub_reftable,projectors=projectors, 
                                                is_trainset=TRUE, 
                                                use_oob=TRUE)[,]
            timeStamps <- .update_timeStamps(timeStamps, previous=NULL, proj_stats=Sys.time())
            ##
            jointEDF$cumul_iter <- c(object$logLs$cumul_iter, 
                                     rep(previous_cumul_iter + it +1L, NROW_newrawstats))
            stat.obs <- .project_reftable_raw(attr(stat.obs,"raw_data"),projectors=projectors, use_oob=FALSE)
          } else if (NROW_newrawstats) { # we did not regenerate projectors, but we have new simuls 
            #  to add to the reftable (no need to reproject the old simuls)
            newlogLstats <- .project_reftable_raw(newrawstats, projectors=projectors, 
                                    methodArgs=project_methodArgs, # need to control parall independently for simuls and for project. 
                                    use_oob=FALSE, # bc only newdata are projected and they are not in training set of re-used predictors
                                    ext_projdata=object$projdata
                                    )
            newlogLstats$cumul_iter <- previous_cumul_iter + it +1L
            jointEDF <- try(rbind(object$logLs, newlogLstats), silent=TRUE)
            if (inherits(jointEDF,"try-error")) {
              preexcols <- colnames(object$logLs)
              newcols <- colnames(newlogLstats)
              if ( length(preexcols) != length(newcols) || ! all(preexcols==newcols) ) {
                message("colnames(object$logLs):")
                print(colnames(object$logLs))
                message("colnames(newlogLstats):")
                print(colnames(newlogLstats))
                stop("columns of pre-existing reference table do not match those of of new simulations.")
              } else stop(attr(jointEDF,"condition")$message)
            }
            attr(jointEDF,"projectors") <- attr(newlogLstats,"projectors")
            timeStamps <- .update_timeStamps(timeStamps, previous=NULL, proj_stats=Sys.time()) # most relevant time stamp is when new proj stats were computed,
            # i.e. when project.default is called
          } else { # recluster
            jointEDF <- object$logLs 
            attr(jointEDF,"projectors") <- projectors
          }
          ####
        } else { # NO projectors, possibly recluster
          newlogLstats <- newrawstats
          if (NROW_newrawstats) newlogLstats$cumul_iter <- previous_cumul_iter + it +1L
          jointEDF <- rbind(object$logLs,newlogLstats)
        }
        attr(jointEDF,"allPars") <- object$colTypes$allPars
        if (NROW_newrawstats) {
          latestPoints <- nrow(jointEDF)-NROW_newrawstats + seq_len(NROW_newrawstats) 
        } else latestPoints <- object$latestPoints
        object <- .infer_SLik_joint(data = jointEDF,
                                   stat.obs = stat.obs,
                                   Simulate = get_from(object,"Simulate"),
                                   using=using,
                                   nbCluster=nbCluster,
                                   verbose=verbose$most,
                                   constr_crits=object$constr_crits,
                                   initParam=.get_initParam(object),
                                   latentVars=colTypes$latentVars)
        # (3) maximization, new CIs, new MSEs
        object <- MSL(object, CIs=eval_CIs_here, level=level, verbose=verbose, 
                      eval_RMSEs=eval_RMSEs_here, cluster_args=RMSE_cluster_args)
        if (.is_clustering_suspect(object)) { 
          if (.is_devel_session()) cat(cli::col_cyan("... reclustering as clustering is suspect.\n"))
          ## Condition the use of initParam on the previous clustering not being 'suspect'
          ## The MLE 'belonging" to an unlikely and peaked cluster would be suspect.
          object <- .infer_SLik_joint(data = jointEDF,
                                     stat.obs = stat.obs,
                                     Simulate = get_from(object,"Simulate"),
                                     using=using,
                                     nbCluster=nbCluster,
                                     verbose=verbose$most,
                                     constr_crits=object$constr_crits,
                                     initParam=NULL,
                                     latentVars=colTypes$latentVars)
          object <- MSL(object, CIs=eval_CIs_here, level=level, verbose=verbose, 
                        eval_RMSEs=eval_RMSEs_here, cluster_args=RMSE_cluster_args)
        } 
        ## 
        object$latestPoints <- latestPoints
        if ( ! is.null(projectors)) {
          object$reftable_raw <- reftable_raw
          object$proj_sub_trainset <- prev_sub_trainset # subset of:
          object$proj_trainset <- proj_trainset
        }
        
        if (! is.null(verbose$cloud_parm)) {  # private: the name of the x variable induces the plot
          if ( ! exists("rparam_blob")) rparam_blob <- NULL
          .plot_cloud(
            object=object, latestPoints=latestPoints, rparam_blob=rparam_blob, 
            verbose=verbose, fittedPars=fittedPars)
        }
      } else { # primitive workflow
        if ( ! is.null(projectors <- object$projectors)) {
          newlogLstats <- .project_reftable_raw(newrawstats,projectors=projectors)
        } else newlogLstats <- newrawstats
        if ( ! inherits(newlogLstats,"list") ) {
          stop("'newsimuls' must be a _list_ of empirical distributions.")
        }
        if (inherits(object,"SLikp")) {
          mc <- attr(object$tailp,"call") 
        } else {
          mc <- get_from(object,"call") # call to infer_logLs (primitive workflow !)
        } 
        mc$object <- newlogLstats
        mc$verbose <- verbose$most
        mc$`stat.obs` <- get_from(object,"stat.obs") ## bc otherwise mc$`stat.obs` stores a promise such as 'Sobs'
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
        corrPars1 <- get_ranPars(object$fit,which="corrPars")[["1"]]
        allFix <- c(corrPars1[c("nu","rho")],list(lambda=object$fit$lambda[1],
                                                  phi=object$fit$phi[1],
                                                  beta=fixef(object$fit)))
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
            corrPars1 <- get_ranPars(object$fit,which="corrPars")[["1"]]
            vranfix <- unlist(corrPars1[c("nu","rho")])
            cat(paste(paste(paste(names(vranfix),"=",signif(vranfix)),collapse=", ")," (estimated by ",itmethod,")\n",sep=""))
          }
        }
        object$latestPoints <- nrow(surfaceData)+1-seq_len(nrow(newlogLs)) ## for plots
        # (3) maximization, new CIs, new MSEs
        # ! ! ! ! ! above calls to infer_surface() potentially change rownames of object$logLs relative to those of surfaceDate
        # bc infer_surface() -> make.names(). One solution is to get surfaceData from the new MSL object: 
        object <- MSL(object, CIs=eval_CIs_here, level=level, verbose=verbose, 
                      eval_RMSEs=eval_RMSEs_here, cluster_args=RMSE_cluster_args)
      }    
      # d_maxlogL <- object$MSL$maxlogL - previous_maxlogL
      cumn_here_as_char <- names(target_LR_rparam[target_LR_it])
      # The eval(parse)) allows using "Inf" as name:
      if (cumn>=as.integer( eval(parse(text=cumn_here_as_char)))) target_LR_it <- target_LR_it+1L
      if (verbose$movie) {
        if (verbose$debug) {
          plot(object, from_refine=TRUE, ...) # from which RMSEs may be NULLified   # map.asp may be passed through these dots
        } else {
          chk <- try(plot(object, from_refine=TRUE, ...), # from which RMSEs may be NULLified   # map.asp may be passed through these dots
                     silent=TRUE)
          if (inherits(chk,"try-error")) {
            warning("Figure not produced because error\n'",
                    attr(chk,"condition")$message,
                    "'\n occurred, but computation continues.")
          }         
        }
        if (! is.null(object$MSL$init_from_prof)) {
          cat("Better likelihood found. Maximizing again...\n")
          object <- MSL(object, CIs=eval_CIs_here, level=level, verbose=verbose$most, 
                        eval_RMSEs=eval_RMSEs_here, cluster_args=RMSE_cluster_args)
        } 
      }
      if (inherits(object,"SLik_j")) { object$rparam_info <- rparam_info} 
      RMSEs <- get_from(object,"RMSEs") # LR_RMSE+ logLik() MSE ## 
      target_reached <- ( length(RMSEs) && 
                            ( ! (anyNA_RMSE <- anyNA(RMSEs))) && # but currently (with the reftable method at least) there is no NA in the RMSE table
                            all(RMSEs<precision)) # ___F I X M E___?  better condition if not RMSE for CI bounds ?
      if (eval_CIs_here) eval_CIs_it <- eval_CIs_it+1L
      if (eval_RMSEs_here) eval_RMSEs_it <- eval_RMSEs_it+1L
      it <- it+1L
      it_index <- min(it+1L,maxit)
      newsimuls <- NULL ## essential for the test  ! is.null(newsimuls) && maxit>1L
    } ## end while() loop
    
    if (verbose$most && it < maxit && target_reached) cat("\nIterations terminated because all(RMSEs<precision) where precision =",precision,"\n")

    if (( ! verbose$movie) && verbose$final) {
      if (verbose$debug) {
        plot(object, from_refine=TRUE, ...) # from which RMSEs may be NULLified   # map.asp may be passed through these dots
      } else {
        chk <- try(plot(object, from_refine=TRUE, ...), # from which RMSEs may be NULLified   # map.asp may be passed through these dots
                   silent=TRUE)
        if (inherits(chk,"try-error")) {
          warning("Figure not produced because error\n'",
                  attr(chk,"condition")$message,
                  "'\n occurred, but computation continues.")
        }         
      }
      if (! is.null(object$MSL$init_from_prof)) {
        cat("Better likelihood found. Maximizing again...\n")
        object <- MSL(object, CIs=eval_CIs_here, level=level, verbose=verbose$most, 
                      eval_RMSEs=eval_RMSEs_here, cluster_args=cluster_args)
        if (inherits(object,"SLik_j")) { object$rparam_info <- rparam_info}
      } 
    }
    if (verbose$most) { # as meant by the doc.) 
      if  ( any(is.na(RMSEs))) {
        cat("Target precision could not be ascertained.\n")
      } else if (is.null(get_from(object,"RMSEs"))) {
        cat("Precision has not been evaluated. Use argument 'eval_RMSEs=TRUE' to force its computation.\n")
      } else if  ( any(RMSEs>precision) ) {
        cat("\nTarget precision does not appear to be reached.\n")
      } else cat("\nTarget precision appears to be reached.\n")
    }
    object$timeStamps <- list2env(timeStamps, parent=emptyenv())
    return(object)
  }
