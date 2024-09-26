project <- function(x, ...) UseMethod("project")

.train_cP_size_fn <- function(method_string,stats) {
  if (method_string %in% c("REML")) {
    400
  } else if (method_string == "ranger") { # but fn may not be called in that case ?
    Inf # floor(10000*log2(length(stats)+1)) 
  } else if (method_string %in% c("keras","fastai")) {
    Inf
  } else {floor(1000*log2(length(stats)+1))}
}

.trainingsize_fn <- function(method_string,stats) {
  if (method_string =="REML") {
    Infusion.getOption("projKnotNbr")
  } else if (method_string == "ranger") { # but fn may not be called in that case ?
    Inf # floor(10000*log2(length(stats)+1)) 
  } else if (method_string %in% c("keras","fastai")) {
    Inf
  } else {floor(1000*log2(length(stats)+1))}
}



.fastai <- function(data, stats, x) {
  nr <- nrow(data)
  tot <- 1:nr
  tr_idx <- sample(nr, 0.8 * nr)
  ts_idx <- tot[!tot %in% tr_idx]
  ###        TabularDataTable(data, procs = list(), cat_names=c(), cont_names=stats, y_names = x, splits = list(tr_idx, ts_idx) )
  locfn <- get("TabularDataTable", asNamespace("fastai"))
  args <- list(df=data, procs=list(), cat_names=c(), cont_name =stats,
               y_names = x, splits = list(tr_idx, ts_idx) )
  tdt <- do.call(locfn,args) # opaque data formatting
  ###         dls = tdt %>% dataloaders(bs = 50)
  locfn <- get("dataloaders", asNamespace("fastai"))
  args <- list(object=tdt, bs=50L)
  dls <- do.call(locfn,args) 
  ###         projector = dls %>% tabular_learner(layers=c(200,100,100,200), config = tabular_config(embed_p = 0.3, use_bn = FALSE), metrics = list(accuracy))
  #
  locfn <- get("tabular_config", asNamespace("fastai"))
  args <- list(embed_p = 0.3, use_bn = FALSE)
  config <- do.call(locfn,args) 
  #
  locfn <- get("tabular_learner", asNamespace("fastai"))
  accuracy. <- get("accuracy", asNamespace("fastai"))
  args <- list(dls=dls, config = config, metrics = list(accuracy.)) # metrics=list() results in an error in project.default-> predict() -> ... 
  projector <- do.call(locfn,args) 
  ###          lrs <- capture.output(projector %>% lr_find()) 
  locfn <- get("lr_find", asNamespace("fastai"))
  args <- list(object=projector)
  lrs <- capture.output(do.call(locfn,args)) 
  ###          projector %>% fit(5, lr = 0.005)
  lr <- as.numeric(strsplit(strsplit(lrs,"=")[[1]][2], ",")[[1]][1])
  locfn <- get("fit", asNamespace("generics")) # calling "fit.fastai.tabular.learner.TabularLearner" is not sufficient 
  args <- list(object=projector, 5, lr=lr) 
  abyss <- do.call(locfn,args) 
  ###
  attr(projector,"stats") <- stats 
  # base::plot(cbind(data[[x]],projector %>% predict(data)))
  projector
}

.keras <- function(data, stats, x, methodArgs) {
  method <- get("fit.keras.engine.training.Model", asNamespace("keras"))
  build_model <- methodArgs$build_model
  if (is.null(build_model)) build_model <- .default_build_model
  formalNames <- names(formals(method))
  methodArgs[setdiff(names(methodArgs),formalNames)] <- NULL ## so that only arguments of 'method' remain in methodArgs
  if ("x" %in% formalNames && is.null(methodArgs$x)) methodArgs$x <- as.matrix(data[,stats,drop=FALSE]) 
  if ("y" %in% formalNames && is.null(methodArgs$y)) methodArgs$y <- data[,x] 
  methodArgs$object <- do.call(build_model, args=list(input_shape=dim(methodArgs$x)[[2]]))
  if (is.null(methodArgs$epochs)) methodArgs$epochs <- 20L 
  abyss <- do.call("fit",methodArgs) # the LHS is not the projector
  projector <- methodArgs$object
  attr(projector,"stats") <- stats
  list(projector=projector, methodArgs=methodArgs)
}


project.character <- function(x,
           stats,
           data,
           trainingsize= eval(Infusion.getOption("trainingsize")),
           train_cP_size= eval(Infusion.getOption("train_cP_size")),
           method, methodArgs=eval(Infusion.getOption("proj_methodArgs")), 
           #nb_cores= Infusion.getOption("nb_cores"), # but see methodArgs$num.threads for ranger
           verbose=TRUE,
           keep_data=TRUE,
           ...) {
    if (x %in% stats) stop(paste0("Given parameter name '",x,"' is one of the statistics' names."))
    if (missing(method)) {
      method <- "ranger"
      #if (requireNamespace(method,quietly=TRUE)) {
        if ( ! .one_time_warnings$RF_warned ) {
          message("Selecting 'ranger' as default method ")
          .one_time_warnings$RF_warned <- TRUE
        }
      # } else {
      #   if ( ! .one_time_warnings$RF_warned ) {
      #     message(paste("  If the 'ranger' package were installed, fast projection would be possible.\n",
      #                   "  Instead, a slow method (REML) will be used."))
      #     .one_time_warnings$RF_warned <- TRUE
      #   }
      #   method <- "REML"
      # }
    }
    if (!is.character(method)) {
      stop("'method' should be a character string, \n in particular a function name rather than a function")
    } else method_string <- method
    if (is.data.frame(data)) {
      totsim <- data
    } else if (inherits(data,"list")) {
      pardata <- lapply(data,function(ll) {
        cbind(ll,attr(ll,"par")[x],row.names=NULL) ## attr(ll,"par") is a (1-row) data.frame and cbind(matrix,data.frame) is data.frame!
      }) ## puts the parameter values into the data
      totsim <- do.call(rbind,pardata)
    } else {
      stop("'data' are neither a data frame nor a list of data frames.")
    }
    # ranger does not generally handles missing data, hence we need this: (hopefully without sideeffects on the use of row names)
    raw_form <- as.formula(paste(x," ~ ",paste(stats,collapse="+")))
    totsim <- model.frame(raw_form, data=totsim)
    #
    rownames(totsim) <- make.names(rownames(totsim),unique = TRUE)
    nr <- nrow(totsim)
    if (nr>trainingsize) {
      message("projection data reduced according to 'trainingsize' argument.")
      trainsample <- sample(nr,trainingsize)
    } else trainsample <- seq(nr)
    if (method_string %in% c("GCV","REML")) {
      form <- as.formula(paste(x," ~ 1 + Matern(1|",paste(stats,collapse="+"),")"))
      # utiliser totsim dans ce qui suit => plantage severe
      if (method_string=="REML") {
        ## hummmm prend ~ 3 minutes pour 300 points
        ## pas de replicat de X available for projection -> estim phi/lambda pb et pas de GCV possible...
        if ((nrs <- length(trainsample))>train_cP_size) {
          train_cP_sample <- trainsample[sample(nrs,train_cP_size)] ## row indices of totsim...
        } else train_cP_sample <- trainsample
        essfit <- fitme(form,data=totsim[train_cP_sample,], fixed=list(nu=4), 
                          method=method, init=list(rho=rep(NA,length(stats)))) 
        corrPars1 <- get_ranPars(essfit,which="corrPars")[["1"]]
        ranfix <- c(corrPars1,list(lambda=essfit$lambda,phi=essfit$phi))
      } else if (method_string=="GCV") {## no difference between trainingSize and knotNbr
        forGCV <- prepareData(data=totsim[trainsample,],ParameterNames=stats,
                              respName=x)
        oldsmoo <- list(minSmoothness=blackbox.getOption("minSmoothness"),maxSmoothness=blackbox.getOption("maxSmoothness"))
        do.call(blackbox.options,list(minSmoothness=4,maxSmoothness=4)) ## all the more important if there are no replicates 
        gcvres <- calcGCV(forGCV)
        ## reestimate lambda and phi (otherwise one should take into account that lambda has different meaning in GCV and hglm notations)
        ranfix <- list(rho=1/gcvres$CovFnParam[stats],
                       nu=gcvres$CovFnParam[["smoothness"]]) 
        do.call(blackbox.options,oldsmoo) ## restaure previous values 
      }    
      #
      ## FR->FR heteroskedas ?
      projector <- fitme(form,data=totsim[trainsample,], fixed=ranfix, method="REML") 
    } else if (method_string=="fastai") {
      nr <- nrow(totsim)
      tot <- 1:nr
      tr_idx <- sample(nr, 0.8 * nr)
      ts_idx <- tot[!tot %in% tr_idx]
      ###        TabularDataTable(totsim, procs = list(), cat_names=c(), cont_names=stats, y_names = x, splits = list(tr_idx, ts_idx) )
      locfn <- get("TabularDataTable", asNamespace("fastai"))
      args <- list(df=totsim, procs=list(), cat_names=c(), cont_name =stats,
                   y_names = x, splits = list(tr_idx, ts_idx) )
      tdt <- do.call(locfn,args) # opaque data formatting
      ###         dls = tdt %>% dataloaders(bs = 50)
      locfn <- get("dataloaders", asNamespace("fastai"))
      args <- list(object=tdt, bs=50L)
      dls <- do.call(locfn,args) 
      ###         projector = dls %>% tabular_learner(layers=c(200,100,100,200), config = tabular_config(embed_p = 0.3, use_bn = FALSE), metrics = list(accuracy))
      #
      locfn <- get("tabular_config", asNamespace("fastai"))
      args <- list(embed_p = 0.3, use_bn = FALSE)
      config <- do.call(locfn,args) 
      #
      locfn <- get("tabular_learner", asNamespace("fastai"))
      accuracy. <- get("accuracy", asNamespace("fastai"))
      args <- list(dls=dls, config = config, metrics = list(accuracy.)) # metrics=list() results in an error in project.default-> predict() -> ... 
      projector <- do.call(locfn,args) 
      ###          lrs <- capture.output(projector %>% lr_find()) 
      locfn <- get("lr_find", asNamespace("fastai"))
      args <- list(object=projector)
      lrs <- capture.output(do.call(locfn,args)) 
      ###          projector %>% fit(5, lr = 0.005)
      lr <- as.numeric(strsplit(strsplit(lrs,"=")[[1]][2], ",")[[1]][1])
      locfn <- get("fit", asNamespace("generics")) # calling "fit.fastai.tabular.learner.TabularLearner" is not sufficient 
      args <- list(object=projector, 5, lr=lr) 
      abyss <- do.call(locfn,args) 
      ###
      attr(projector,"stats") <- stats 
      # base::plot(cbind(data[[x]],projector %>% predict(data)))
    } else if (method_string=="keras") {
      blob <- .keras(data=totsim[trainsample,], stats, x, methodArgs)
      projector <- blob$projector
      methodArgs <- blob$methodArgs # modified by .keras()
    } else {
      if (method_string=="ranger") { # availability already checked
        method <- ranger 
      } else method <- eval(parse(text=method)) 
      formalNames <- names(formals(method))  
      methodArgs[setdiff(names(methodArgs),formalNames)] <- NULL ## so that only arguments of 'method' remain in methodArgs
      dotargs <- match.call(expand.dots = FALSE)$... ## produce a pairlist of (essentially) promises. No quote() needed
      for (st in names(dotargs)) methodArgs[[st]] <- dotargs[[st]]
      #
      methodArgs_for_reproject <- methodArgs # save original copy before modifs
      #
      if ("formula" %in% formalNames && is.null(methodArgs$formula)) methodArgs$formula <- raw_form 
      if ("data" %in% formalNames && is.null(methodArgs$data)) methodArgs$data <- totsim[trainsample,] # for ranger, notably
      if (method_string=="ranger") {
        if (is.null(methodArgs$importance)) {
          methodArgs$importance <- "permutation"
          methodArgs_for_reproject$importance <- "none"        
        }
        if (is.null(methodArgs$num.threads)) methodArgs$num.threads <- 
            max(1L,Infusion.getOption("nb_cores")) # ranger's default is to use all cores!!!  => NULL is distinct from 1L for ranger
        if ( ! .one_time_warnings$nThr_warned && methodArgs$num.threads==1L && nrow(methodArgs$data)>2000L ) {
          message(paste("Parallelisation might be useful for ranger() calls. See e.g. 'cluster_args' argument of refine().\n"))
          .one_time_warnings$nThr_warned <- TRUE
        }
        if (is.null(methodArgs$splitrule)) {
          if ("splitrule" %in% names(methodArgs)) { 
            warning("Explicit NULL 'splitrule' in 'methodArgs' is interpreted as 'splitrule=\"extratrees\"', contrary to the ranger() default.'",
                    immediate.=TRUE)
            methodArgs$splitrule <- "extratrees" 
          } else methodArgs$splitrule <- "extratrees" # Implicit NULL also has a different meaning for ranger and for Infusion 
        }
        # ---- 'second best' in alex's simulation study, yet providing oob predictions: (replace <- FALSE with sample.fraction <- 1 would prevent them)
        if (is.null(methodArgs$mtry)) {
          if (methodArgs$splitrule=="extratrees") {
            mtry <- function(n) n ## often better in Alex's simulations
          } else mtry <- function(n) floor(n/3) ## Breiman's advice for *regression*
        }
        if (is.null(methodArgs$replace)) {
          if (methodArgs$splitrule=="extratrees") {
            methodArgs$replace <- FALSE
          } else methodArgs$replace <- TRUE
        }
        if (is.null(methodArgs$sample.fraction)) {
          if (methodArgs$splitrule=="extratrees") {
            methodArgs$sample.fraction <- 0.632
          } 
        }
        if (is.null(methodArgs$num.trees)) methodArgs$num.trees <- 1000L
      } else {
        # ranger allows 'data'+'formula' or 'x'+'y' as alternative sets of arguments. Here we ignore 'x'+'y'.
        if ("x" %in% formalNames && is.null(methodArgs$x)) methodArgs$x <- totsim[trainsample,stats,drop=FALSE] 
        if ("y" %in% formalNames && is.null(methodArgs$y)) methodArgs$y <- totsim[trainsample,x] ## drop is important...
      }
      projector <- do.call(method,methodArgs)
      if ("formula" %in% formalNames) attr(projector,"stats") <- stats ## otherwise NULL, see use in project.default
    }
    if (verbose) {
      OKplot <- FALSE
      if (inherits(projector,"HLfit")) {
        yy <- projector$data[,x]
        OKplot <- TRUE
      } else if (inherits(projector,"train")) {
        yy <- projector$trainingData$.outcome
        OKplot <- TRUE
      } else if (inherits(projector,"ranger")) {
        yy <- totsim[trainsample,x] # but the calls could be used ?    
        OKplot <- TRUE
      } else if (inherits(projector,"randomForest")) {
        yy <- totsim[trainsample,x] ## (which is not returned in the randomForest object) # but the calls could be used
        OKplot <- TRUE
      } else if (inherits(projector,"fastai.tabular.learner.TabularLearner")) {
        yy <- totsim[trainsample,x]
        OKplot <- TRUE
      } else if (inherits(projector,"keras.engine.training.Model")) {
        yy <- methodArgs$y
        OKplot <- TRUE
      } 
      if (OKplot) { ## should try to plot oob.
        if (inherits(projector,"fastai.tabular.learner.TabularLearner")) {
          chk <- try(plot(predict(projector, totsim[,stats,drop=FALSE])[,1],yy,xlab=paste("predicted",x),ylab=paste("true",x)))
        } else if (inherits(projector,"keras.engine.training.Model")) {
          chk <- try(plot(predict(projector, x=methodArgs$x),yy,xlab=paste("predicted",x),ylab=paste("true",x)))
        } else if (inherits(projector,"ranger")) {
          chk <- try(plot(projector$predictions,yy,xlab=paste("oob-predicted",x),ylab=paste("true",x))) 
        } else if (inherits(projector,"randomForest")) {
          chk <- try(plot(predict(projector),yy,xlab=paste("oob-predicted",x),ylab=paste("true",x)))
        } else chk <- try(plot(predict(projector,newdata=totsim[trainsample,stats,drop=FALSE]),yy,
                               xlab=paste("predicted",x),ylab=paste("true",x)))
        if ( ! inherits(chk, "try-error")) abline(0,1)
      }
    }
    if (method_string=="ranger") {
      projector$forest$child.nodeIDs <- lapply(projector$forest$child.nodeIDs,function(li) lapply(li,as.integer))
      projector$forest$split.varIDs <- lapply(projector$forest$split.varIDs,as.integer)
      if ( ! keep_data) projector$call$data <- "removed"
    }
    
    ###### 
    ## The project_call attribute is the info used by reproject() or .update_projector()
    ## It has zero-line $data, and called function 'project.character' rather than the method (this actually saves memory).
    project_call <- match.call()
    if (is.data.frame(data)) project_call$data <- data[FALSE,] ## <0 rows> data frame that keeps the colnames info that will be useful for reprojection
    ## : where data is not always a data frame (found while testing obsolete examples, specifically nnet).
    ## We also need to secure x = parName, stats = statNames too:
    ## (Next lines do not work: # does not work: this may store x=local value of "somevar" but the promise remains x= somevar where somevar is not a variable in the list  
    #    call_env <- list2env(lapply(as.list(match.call(expand.dots=TRUE))[-1],eval)) 
    #    environment(mc) <- call_env
    ## )
    project_call$x <- x # the focal parameter names
    project_call$stats <- stats
    project_call$method <- method_string # not its interpretation as a function!
    project_call$trainingsize <- trainingsize
    project_call$train_cP_size <- train_cP_size
    if (method_string=="ranger") project_call$methodArgs <- methodArgs_for_reproject # no data, in particular, otherwise future data would not be taken into account. 
    attr(projector,"project_call") <- project_call 
    ######
    
    return(projector) 
  }

# wrapper initially for handling names but now also with oob/trainset control
.predictWrap <- function(oneprojector,newdata, use_oob=Infusion.getOption("use_oob"), 
                         is_trainset=FALSE, methodArgs=NULL,
                         ext_projdata,
                         ...) {
  if (inherits(newdata,c("numeric","integer"))) { ## not data.frame...
    ## utils:::str.default(oneprojector) -> function with attributes py_object and project_call
    stats <- attr(oneprojector,"stats") ## should be non-null if projector.character used a formula
    ## or for keras (whose predict method does not handle extra columns!)
    ## else names should (?) not be required and then stats is ideally NULL
    if( ! is.null(stats)) {
      datanames <- names(newdata)
      dim(newdata) <- c(1,length(newdata)) ## converts to 1-row matrix faster than t(x)
      if (is.null(datanames)) {
        colnames(newdata) <- stats
      } else {
        colnames(newdata) <- datanames
      }
    }   
    # now still numeric, but matrix and sure to have colnames
  } ## else colnames should already be set
  if (inherits(oneprojector,"keras.engine.training.Model")) {
    if (is.data.frame(newdata)) {
      stats <- attr(oneprojector,"stats") 
      newdata <- newdata[,stats,drop=FALSE]
      newdata <- as.matrix(newdata)
    }
    predict(oneprojector,x=newdata,...)
  } else if (inherits(oneprojector,"fastai.tabular.learner.TabularLearner")) {
    # vector has been converted to matrix
    newdata <- as.data.frame(newdata)
    predict(oneprojector,newdata,...)
  } else if (inherits(oneprojector,"ranger")) { 
    # use_oob is TRUE by default
    if (is_trainset && use_oob) {
      pred <- oneprojector$predictions
    } else {
      num.threads <- methodArgs$num.threads
      if (is.null(num.threads)) num.threads <- max(1L,Infusion.getOption("nb_cores")) # ranger's default is to use all cores!!! 
                                              # => NULL is distinct from 1L for ranger, and must be avoided.
      pred <- predict(oneprojector,data=newdata,
                      num.threads=num.threads, 
                      ...)$predictions 
      
      if (use_oob &&  ! is.null(dim(newdata))) { # costly check so (is_trainset && use_oob) should be used whenever TRUE
        # not-oob predictions are retained for:
        #  * possibly part of the reftable if it is not fully used for training. 
        #   (=> e.g. if projections are not updated in a iteration: common case when the ranger data information must be kept  
        #     as it is distinct from the SLik data) => next line is quite speculative
        if ( is.null(dim(projdata <- oneprojector$call$data))) projdata <- ext_projdata
        if ( is.null(projdata)) stop("Neither 'oneprojector$call$data' nor 'ext_projdata' are available.")
        #  * the stat.obs; and 
        #  * new simulations from goftest
        x <- projdata[,oneprojector$forest$independent.variable.names, drop=FALSE] # data.frame `==`
        #(without the drop arg, a 1-col df is reduced to a vector => pb for toy example projecting a 1-dim stat)
        # removed awful old dist()-based code that memory-failed for large matrices...
        # Next line is slow hence 'use_oob' and 'is_trainset' controls were implemented.
        posinold <- match(data.frame(t(newdata[,colnames(x),drop=FALSE])),data.frame(t(x)) ) # https://stackoverflow.com/questions/12697122/in-r-match-function-for-rows-or-columns-of-matrix
        newinold <- na.omit(posinold)
        if (length(newinold)) {
          oob_positions <-  ! is.na(posinold)
          pred[oob_positions] <- oneprojector$predictions[na.omit(posinold)] # using out-of-bag predictions.
          attr(pred,"oob_positions") <- oob_positions
        }
      }
    }
    return(pred)
  } else if (inherits(oneprojector,"randomForest")) { 
    pred <- predict(oneprojector,newdata=newdata,...)
    
    if (use_oob &&  ! is.null(dim(newdata))) {
      x <- oneprojector$call$x
      posinold <- match(data.frame(t(newdata[,colnames(x),drop=FALSE])),data.frame(t(x)) ) # https://stackoverflow.com/questions/12697122/in-r-match-function-for-rows-or-columns-of-matrix
      newinold <- na.omit(posinold)
      if (length(newinold)) {
        oobpred <- predict(oneprojector,...) # ie oneprojector$predicted
        pred[ ! is.na(posinold)] <- oobpred[na.omit(posinold)] 
      } 
    }
    return(pred)
  } else {
    predfn <- getS3method("predict",class(oneprojector)[1L])
    if ("newdata" %in% names(formals(predfn))) {
      resu <- predict(oneprojector,newdata=newdata,...)
    } else resu <- predict(oneprojector,data=newdata,...) 
    resu <- unname(resu) # caret::predict.train found to return dummy names 2023/07...
  }
    
}

# x attributes used late in code, hence x should not be modified
.project_reftable_raw <- function (x, projectors, use_oob=Infusion.getOption("use_oob"), is_trainset=FALSE, 
                             methodArgs=list(), ext_projdata, 
                             ...) {
  #
  if (inherits(projectors,"list")) {
    projectors <- list2env(projectors,
                           parent=emptyenv())
  } else if (! inherits(projectors,"environment")) stop("'projectors' argument must be a environment or a list.")
  #
  if (inherits(x,"list")) { ##for list of EDFs where pars are not in the table (ie old workflow, not a reftable)
    ly <- lapply(x, function(lt) {
      class(lt) <- c(class(lt),"EDF")
      tmp <- .project_reftable_raw(lt,projectors=projectors)
      attr(tmp,"par") <- attr(lt,"par")
      tmp
    })
    ## copy 'infusion-specific' attributes from x to ly
    attrx <- attributes(x) # list
    attrx["names"] <- NULL
    attrx["class"] <- NULL
    for (li in names(attrx)) attr(ly,which=li) <- attrx[[li]]
    attr(ly,"projectors") <- projectors 
    return(ly)
  } #### ELSE 
  #
  if (is.vector(x)) {
    checknames <- names(x)
  } else {
    checknames <- colnames(x)
    # check moved from another position, all conditions may may no longer be most adjusted
    if ( ! inherits(x,"EDF")) {
      if ( ( ! is.null(dim(x))) && is.null(attr(x,"na.action"))) { # is.null(...) possible in direct call of project(<projector>, x=<new data>) which is not basic workflow
        x <- na.omit(x)
        if ( ! is.null(attr(x,"na.action"))) warning("project() called on data including NA's: it may be better to apply na.omit() on the data first.")
      }
    }
  }
  #
  projnames <- names(projectors)
  ly <- vector("list", length(projnames))
  names(ly) <- projnames
  for (st in projnames) {
    if (is.null(projectors[[st]])) { # raw statistics to be kept # for retained raw statistics, explicit NULLs are needed in the projectors list, 
      # otherwise the retained raw statistics are not found in 'projnames' and the raw stat is lost from the result.
      if (is.vector(x)) {
        ly[[st]] <- x[[st]]
      } else {
        ly[[st]] <- x[,st]
      }
    } else {
      if (st %in% checknames) {
        message(paste0("'",st,"' appears to name both a projector and a raw statistics. This suggests either that 'projectors' are misnomed"))
        message("   (they should receive names different from those of original summary statistics);")
        message("  or that 'project...' is applied recursively ('x' should not be the result of applying 'project...')")
        stop("From 'project.default': check names of 'x' and 'projectors'. I exit.")
      }
      #
      if (inherits(x,"EDF")) { ## For old workflow, not a reftable but a list of EDFs for given param values: the was tag assigned by the parent project.default function (see code if (inherits(x,"list"))...)
        # EDF => stats names are already there
        ly[[st]] <- .predictWrap(projectors[[st]], newdata=as.data.frame(x)) 
      } else {
        ly[[st]] <- .predictWrap(projectors[[st]], newdata=x, use_oob=use_oob, is_trainset=is_trainset,
                     methodArgs=methodArgs, ext_projdata=ext_projdata) # ordered as names(projectors), not ls(projectors)
      }
      
    }
    
  }
  
  if (is.vector(x)) {
    ly <- unlist(ly)
  } else { ## Result expected to be a matrix (SLik case), or data.frame (SLikj case)
    ly <- do.call(cbind, ly)
    if (inherits(x,"data.frame")) { ## projection of reftable
      inferredVars <- names(attr(x,"LOWER")) # where latentVars must have been included using declare_latent() 
      ly <- cbind(x[,inferredVars,drop=FALSE],ly) ## automatic conversion of ly to data.frame
      if ( inherits(x, "reftable") && ! inherits(ly, "reftable")) class(ly) <- c("reftable", class(ly)) 
    }
    attr(ly,"Simulate") <- attr(x,"Simulate")
    attr(ly,"control.Simulate") <- attr(x,"control.Simulate")
    # attr(ly,"Simulate_input") <- attr(x,"Simulate_input")
    attr(ly,"packages") <- attr(x,"packages")
    attr(ly,"env") <- attr(x,"env") ## an environment!
    attr(ly,"workflow_env") <- attr(x,"workflow_env") # copy without modif 
    attr(ly,"projectors") <- projectors # an environment
    attr(ly,"LOWER") <- attr(x,"LOWER") 
    attr(ly,"UPPER") <- attr(x,"UPPER") 
  }
  attr(ly,"raw_data") <- x 
  # We have here 'raw_data' (data from the point of view of the projection) attached to the return value
  # which is either a reftable OR the *data* (stat.obs)  of the Infusion inference problem.
  # For stat.obs, 'raw_data' will remain attached to it, while for reftable, they will be moved to the SLik_j's $reftable_raw
  # For primitive-workflow SLik's the reftable_raw equivalent is attr(logLs, 'raw_data').
  # Whatever the class of the objects and then non-redundancy of stat.obs vs reftable, 
  ly ## same type as input (as doc'ed)
}

# user-level version of .project_reftable_raw: only distinction is that of keeping the latentVars attr
project.default <- function (x, projectors, use_oob=Infusion.getOption("use_oob"), is_trainset=FALSE, 
                                                methodArgs=list(), ext_projdata, # in a programming context it's better to take them from colTypes
                                                ...) {
  mc <- match.call(expand.dots=TRUE) 
  mc[[1L]] <- get(".project_reftable_raw", asNamespace("Infusion"), inherits=FALSE) 
  resu <- eval(mc,envir = parent.frame())
  #
  attr(resu,"latentVars") <- attr(x,"latentVars")
  resu
}
  

.update_reft_raw <- function(olddata, newdata) {
  mostAttrs <- names(attributes(olddata))
  mostAttrs <- setdiff(mostAttrs,c("class","dim","dimnames","names","row.names"))
  newdata <- rbind(olddata,newdata)
  for (attrname in mostAttrs) {
    attr(newdata,attrname) <- attr(olddata,attrname)
  }
  if (is.null(cumul_n <- attr(olddata, "cumul_n"))) cumul_n <- c(0L, nrow(olddata))
  attr(newdata, "cumul_n") <- c(cumul_n, nrow(newdata))
  return(newdata)
}

get_projector <- function(...) project.character(...) 
get_projection <- function(...) project.default(...) 

neuralNet <- function(formula,data) {
  if (isNamespaceLoaded("doSNOW")) {
    unloadNamespace("doSNOW")
    on.exit(do.call("loadNamespace",list(package="doSNOW"))) ## ...but not the outer one
  }
  .do_call_wrap("train", list(form=formula, data=data, method='nnet', linout=TRUE, trace = FALSE), pack="caret")
}


# "plot_proj" <- function(object, ...) UseMethod("plot_proj") ## makes it easy to develop new inference methods

plot_proj <- function(object, # SLik_j ; range object not handled, see comment in commented code.
                      parm=NULL, proj, 
                      new_rawdata=NULL, 
                      use_oob=Infusion.getOption("use_oob"), 
                      is_trainset=FALSE, 
                      xlab=NULL, ylab=NULL, 
                      ...) {
  
  if (inherits(object,"ranger")) {
    warning("Experimental use of plot_proj() on ranger object: 'new_rawdata' needed and must be the training set.")
    parm <- paste(object$call$formula[[2]])
    x <- new_rawdata[,parm]
    y <- object$predictions
  } else {
    projectors <- object$projectors
    if (is.null(parm)) {
      if (missing(proj)) {
        projnames <- ls(projectors)
        np <- length(projnames)
        if (np>3L) {
          mai <- c(0.45,0.5,0.1,0.1)
        } else mai <- NULL
        opar <- .fittedPars_plot_pars(np=np, plotpar=list(), ylab=ylab, mai=mai)
        resu <- setNames(vector("list", np), projnames)
        for (st in projnames) {
          resu[[st]] <- plot_proj(object, # SLik_j ; range object not handled, see comment in commented code.
                    proj=st, 
                    new_rawdata=new_rawdata, 
                    use_oob=use_oob, 
                    is_trainset=is_trainset, 
                    xlab=xlab, ylab=ylab, 
                    ...)
        }
        par(opar)
        return(invisible(resu))
      } else {
        projector <- projectors[[proj]] # there must be proj
        parm <- attr(projector,"project_call")$x
      }
    } else {
      projnames <- ls(projectors)
      pars <- sapply(projnames, function(st) attr(projectors[[st]],"project_call")$x) 
      proj <- names(which(pars==parm))
      projector <- projectors[[proj]]
    }
    # x and y chosen so that E[y] ~ x ('regression of true value on predicted one') which is how RF are trained 
    # (conversely one does NOT see x as the mean of the y values)  
    timeStamps <- object$timeStamps
    is_timeStamps_suspect <-  (! is.null(proj_timeStamp <- timeStamps$projectors) && 
                              ! is.null(logLs_timeStamp <- timeStamps$proj_stats) &&
                              proj_timeStamp>logLs_timeStamp )
    if ( is.null(new_rawdata)) {
      if (is_timeStamps_suspect) {
        warning(paste("Projectors were updated after projected statistics were computed:\n",
                      "plot will not represent properties of any fit object.\n",
                      "See Details of ?plot_proj for more information."))
      }
      logLs <- object$logLs
      y_parm <- logLs[,parm] # parm value
      x_pred <- logLs[,proj] # predictions
      oob_positions <- object$proj_trainset
    } else {
      if (is_timeStamps_suspect) {
        message(paste("Projectors were updated after projected statistics were computed:\n",
                      "plot will represent properties of updated projector, not of original one.\n",
                      "See Details of ?plot_proj for more information."))
      }
      y_parm <- new_rawdata[,parm]
      x_pred <- .predictWrap(projector, newdata=new_rawdata,
                        use_oob=use_oob, is_trainset=is_trainset, 
                        ext_projdata=object$projdata)
      oob_positions <- attr(y,"oob_positions")
      attr(y,"oob_positions") <- NULL
    }
  }
  if (is.null(xlab)) xlab <- paste("predicted values of", parm)
  if (is.null(ylab)) ylab <- paste("true", parm) 
  
  plot(x=x_pred,
       y=y_parm, 
       xlab=xlab, ylab=ylab,...) 
  if ( ! is.null(oob_positions)) points(x_pred[oob_positions], y_parm[oob_positions], col="blue")
  abline(0,1)
  invisible(list(x=x_pred,y=y_parm, xlab=xlab, ylab=ylab))
}

plot_importance <- function (object, parm, proj, n.var = 30L, xlim=NULL, xlab = "Variable Importance", ylab = "", main="",
                             ...) {
  if (inherits(object,"ranger")) {
    projector <- object
  } else if (inherits(object,"SLik_j")) {
    if (is.null(parm)) {
      projector <- object$projectors[[proj]] # there must be proj
    } else {
      projnames <- ls(object$projectors)
      pars <- sapply(projnames, function(st) paste(object$projectors[[st]]$call$formula[[2]]))
      proj <- names(which(pars==parm))
      projector <- object$projectors[[proj]]
    }
  }
  imp <- projector$variable.importance
  n.var <- min(n.var, length(imp))
  ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
  if (is.null(xlim)) xlim <-  c(0, max(imp))
  dotchart(imp[ord], xlim = xlim, xlab = xlab, ylab = ylab, main=main, ...)
  invisible(imp)
}
