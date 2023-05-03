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


project.character <- local({
  RF_warned <- FALSE
  nThr_warned <- FALSE
  function(x,
           stats,
           data,
           trainingsize= eval(Infusion.getOption("trainingsize")),
           train_cP_size= eval(Infusion.getOption("train_cP_size")),
           method, methodArgs=list(), 
           #nb_cores= Infusion.getOption("nb_cores"), # but see methodArgs$num.threads for ranger
           verbose=TRUE,
           ...) {
    if (x %in% stats) stop(paste0("Given parameter name '",x,"' is one of the statistics' names."))
    if (missing(method)) {
      method <- "ranger"
      #if (requireNamespace(method,quietly=TRUE)) {
        if ( ! RF_warned ) {
          message("Selecting 'ranger' as default method ")
          RF_warned <<- TRUE
        }
      # } else {
      #   if ( ! RF_warned ) {
      #     message(paste("  If the 'ranger' package were installed, fast projection would be possible.\n",
      #                   "  Instead, a slow method (REML) will be used."))
      #     RF_warned <<- TRUE
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
        if (essfit$spaMM.version<"2.4.26") {
          corrPars1 <- essfit$corrPars[["1"]]
        } else corrPars1 <- get_ranPars(essfit,which="corrPars")[["1"]]
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
      methodArgs_for_reproject <- methodArgs
      #
      if ("formula" %in% formalNames && is.null(methodArgs$formula)) methodArgs$formula <- raw_form 
      if ("data" %in% formalNames && is.null(methodArgs$data)) methodArgs$data <- totsim[trainsample,] # for ranger, notably
      if (method_string=="ranger") {
        if (is.null(methodArgs$importance)) methodArgs$importance <- "permutation"
        if (is.null(methodArgs$num.threads)) methodArgs$num.threads <- 
            max(1L,Infusion.getOption("nb_cores")) # ranger's default is to use all cores!!!  => NULL is distinct from 1L for ranger
        if ( ! nThr_warned && methodArgs$num.threads==1L && nrow(methodArgs$data)>2000L ) {
          message(paste("Parallelisation might be useful for ranger() calls. See e.g. 'cluster_args' argument of refine().\n"))
          nThr_warned <<- TRUE
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
          plot(predict(projector, totsim[,stats,drop=FALSE])[,1],yy,xlab=paste("predicted",x),ylab=paste("true",x))
        } else if (inherits(projector,"keras.engine.training.Model")) {
          plot(predict(projector, x=methodArgs$x),yy,xlab=paste("predicted",x),ylab=paste("true",x))
        } else if (inherits(projector,"ranger")) {
          plot(projector$predictions,yy,xlab=paste("oob-predicted",x),ylab=paste("true",x)) 
        } else if (inherits(projector,"randomForest")) {
          plot(predict(projector),yy,xlab=paste("oob-predicted",x),ylab=paste("true",x))
        } else plot(predict(projector,newdata=totsim[trainsample,stats,drop=FALSE]),yy,xlab=paste("predicted",x),ylab=paste("true",x))
        abline(0,1)
      }
    }
    mc <- match.call()
    if (is.data.frame(data)) mc$data <- data[FALSE,] ## <0 rows> data frame that keeps the colnames info that will be useful for reprojection
    # but data is not always a data frame (found while testing obsolete examples, specifically nnet)
    # but we need to secure x = parName, stats = statNames too
    #call_env <- list2env(lapply(as.list(match.call(expand.dots=TRUE))[-1],eval)) # does not work: it may store x=local value of "somevar" but the promise remains x= somevar where somevar is not a variable in the list  
    #environment(mc) <- call_env
    mc$x <- x
    mc$stats <- stats
    mc$method <- method_string # not its interpretation as a function!
    mc$trainingsize <- trainingsize
    mc$train_cP_size <- train_cP_size
    if (method_string=="ranger") mc$methodArgs <- methodArgs_for_reproject # no data, in particular, otherwise future data would not be taken into account. 
    attr(projector,"project_call") <- mc
    return(projector) 
  }
})

# wrapper for handling names 
.predictWrap <- function(oneprojector,newdata, use_oob=Infusion.getOption("use_oob"), 
                         # Currently only the first two arguments are used so the following default is always obeyed and the ... never used.
                         # nb_cores=Infusion.getOption("nb_cores"), 
                         is_trainset=FALSE, methodArgs=NULL,
                         ...) {
  if (inherits(newdata,"numeric")) { ## not data.frame...
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
    # On the training data, this 1st pred is overwritten (inelegantly)
    # not-oob predictions are retained for:
    #  * possibly part of the reftable if it is not fully used for training.
    #  * the stat.obs; and 
    #  * new simulations from goftest
    if (is_trainset && use_oob) {
      pred <- oneprojector$predictions
    } else {
      num.threads <- methodArgs$num.threads
      if (is.null(num.threads)) num.threads <- max(1L,Infusion.getOption("nb_cores")) # ranger's default is to use all cores!!! 
                                              # => NULL is distinct from 1L for ranger, and must be avoided.
      pred <- predict(oneprojector,data=newdata,
                      num.threads=num.threads, 
                      ...)$predictions 
      
      if (use_oob &&  ! is.null(dim(newdata))) {
        #message("oob used")
        x <- oneprojector$call$data[,oneprojector$forest$independent.variable.names, drop=FALSE] # data.frame `==`
        #(without the drop arg, a 1-col df is reduced to a vector => pb for toy example projecting a 1-dim stat)
        # removed awful old dist()-based code that memory-failed for large matrices...
        # Next line is slow hence 'use_oob' and 'is_trainset' controls were implemented.
        posinold <- match(data.frame(t(newdata[,colnames(x),drop=FALSE])),data.frame(t(x)) ) # https://stackoverflow.com/questions/12697122/in-r-match-function-for-rows-or-columns-of-matrix
        newinold <- na.omit(posinold)
        if (length(newinold)) pred[ ! is.na(posinold)] <- oneprojector$predictions[na.omit(posinold)] # using out-of-bag predictions.
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
  }
    
}

# x attributes used late in code, hence x should not be modified
project.default <- function (x,projectors, use_oob=Infusion.getOption("use_oob"), is_trainset=FALSE, 
                             methodArgs=list(), ...) {
  #
  if (inherits(projectors,"list")) {
    projectors <- list2env(projectors)
  } else if (! inherits(projectors,"environment")) stop("'projectors' argument must be a environment or a list.")
  #
  if (inherits(x,"list")) { ##for list of EDFs where pars are not in the table (ie old workflow, not a reftable)
    ly <- lapply(x, function(lt) {
      class(lt) <- c(class(lt),"EDF")
      tmp <- project(lt,projectors=projectors)
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
  }
  if (any(checknames %in% names(projectors))) {
    message("Some names of 'x' and 'projectors' match. This suggests either that 'projectors' are misnomed")
    message("   (they should receive names different from those of original summary statistics);")
    message("  or that 'project...' is applied recursively ('x' should not be the result of applying 'project...')")
    stop("From 'project.default': check names of 'x' and 'projectors'. I exit.")
  }
  #
  if (inherits(x,"EDF")) { ## For old workflow, not a reftable but a list of EDFs for given param values: the was tag assigned by the parent project.default function (see code if (inherits(x,"list"))...)
    # EDF => stats names are already there
    ly <- lapply(projectors,.predictWrap, newdata=as.data.frame(x)) 
  } else {
    if ( ( ! is.null(dim(x))) && is.null(attr(x,"na.action"))) { # is.null(...) possible in direct call of project(<projector>, x=<new data>) which is not basic workflow
      x <- na.omit(x)
      if ( ! is.null(attr(x,"na.action"))) warning("project() called on data including NA's: it may be better to apply na.omit() on the data first.")
    }
    # ly <- vector("list",length(projectors))
    # for (projit in seq_along(projectors)) ly[[projit]] <- .predictWrap(projectors[[projit]],newdata=x) # fails on projectors being an *environment*
    ly <- lapply(projectors, .predictWrap, newdata=x, use_oob=use_oob, is_trainset=is_trainset,
                 methodArgs=methodArgs) # ordered as names(projectors), not ls(projectors)
  }
  ly <- do.call(cbind,ly) ## binding is over projectors
  if (is.vector(x)) {
    ly <- as.numeric(ly)
    names(ly) <- names(projectors)
  } else { ## expects ly to be a matrix (SLik case), or data.frame (SLikj case)
    colnames(ly) <- names(projectors)
    if (inherits(x,"data.frame")) { ## projection of reftable
      ly <- cbind(x[,names(attr(x,"LOWER")),drop=FALSE],ly) ## automatic conversion of ly to data.frame
    }
    # if ( ! is.null(parms <- attr(x,"parms"))) { ## assumes that the new data 'x' has not be modified
    #   ly <- cbind(x[,parms,drop=FALSE],ly)
    #   attr(ly,"parms") <- parms
    #   attr(ly,"stats") <- names(projectors)
    # }
  }
  attr(ly,"Simulate") <- attr(x,"Simulate")
  attr(ly,"control.Simulate") <- attr(x,"control.Simulate")
  # attr(ly,"Simulate_input") <- attr(x,"Simulate_input")
  attr(ly,"packages") <- attr(x,"packages")
  attr(ly,"env") <- attr(x,"env") ## an environment!
  attr(ly,"workflow_env") <- attr(x,"workflow_env") # copy without modif 
  attr(ly,"projectors") <- projectors # an environment
  attr(ly,"raw_data") <- x ## F I X M E but then we duplicate info about its attributes
  ly ## same type as input (as doc'ed)
}

.update_raw_data <- function(olddata, newdata) {
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

plot_proj <- function(object, parm, proj, xlab=parm, ylab=proj, ...) {
  plot(x=attr(object$logLs,"raw_data")[,parm],
       y=predict(object$projectors[[proj]], 
                 data=attr(object$logLs,"raw_data"))$predictions, xlab=xlab, ylab=ylab,...) 
}
