project <- function(x, ...) UseMethod("project")

.trainsize_fn <- function(method,stats) {
  if (method %in% c("REML")) {
    400
  } else NULL
}

.knotnbr_fn <- function(method,stats) {
  if (method %in% c("REML")) {
    Infusion.getOption("projKnotNbr")
  } else if (method == "keras") {
    Inf
  } else {floor(1000*log2(length(stats)+1))}
}

project.character <- local({
  RF_warned <- FALSE
  function(x,
           stats,
           data,
           trainingsize= eval(Infusion.getOption("projTrainingSize")),
           knotnbr= eval(Infusion.getOption("knotnbr")),
           method, methodArgs=list(),
           verbose=TRUE,
           ...) {
    mc <- match.call()
    if (x %in% stats) stop(paste0("Given parameter name '",x,"' is one of the statistics' names."))
    if (missing(method)) {
      method <- "ranger"
      if (requireNamespace(method,quietly=TRUE)) {
        if ( ! RF_warned ) {
          message("Selecting 'ranger' as default method ")
          RF_warned <<- TRUE
        }
      } else {
        if ( ! RF_warned ) {
          message(paste("  If the 'ranger' package were installed, fast projection would be possible.\n",
                        "  Instead, a slow method (REML) will be used."))
          RF_warned <<- TRUE
        }
        method <- "REML"
      }
    }
    if (!is.character(method)) {stop("'method' should be a character string, \n in particular a function name rather than a function")}
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
    rownames(totsim) <- make.names(rownames(totsim),unique = TRUE)
    nr <- nrow(totsim)
    if (nr>knotnbr) {
      message("projection data reduced according to 'knotnbr' argument.")
      randsample <- sample(nr,knotnbr)
    } else randsample <- seq(nr)
    if (is.character(method) && method %in% c("GCV","REML")) {
      form <- as.formula(paste(x," ~ 1 + Matern(1|",paste(stats,collapse="+"),")"))
      # utiliser totsim dans ce qui suit => plantage severe
      if (method=="REML") {
        ## hummmm prend ~ 3 minutes pour 300 points
        ## pas de replicat de X available for projection -> estim phi/lambda pb et pas de GCV possible...
        if ((nrs <- length(randsample))>trainingsize) {
          trainsample <- randsample[sample(nrs,trainingsize)] ## row indices of totsim...
        } else trainsample <- randsample
        if (eval(Infusion.getOption("fitmeCondition"))) {
          essfit <- fitme(form,data=totsim[trainsample,], fixed=list(nu=4), 
                          method=method, init=list(rho=rep(NA,length(stats)))) 
        } else essfit <- corrHLfit(form,data=totsim[trainsample,],init.corrHLfit=list(rho=rep(NA,length(stats))),ranFix=list(nu=4))  
        if (essfit$spaMM.version<"2.4.26") {
          corrPars1 <- essfit$corrPars[["1"]]
        } else corrPars1 <- get_ranPars(essfit,which="corrPars")[["1"]]
        ranfix <- c(corrPars1,list(lambda=essfit$lambda,phi=essfit$phi))
      } else if (method=="GCV") {## no difference between trainingSize and knotNbr
        forGCV <- prepareData(data=totsim[randsample,],ParameterNames=stats,
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
      if (eval(Infusion.getOption("fitmeCondition"))) {
        projector <- fitme(form,data=totsim[randsample,], fixed=ranfix, method="REML") 
      } else projector <- corrHLfit(form,data=totsim[randsample,],ranFix=ranfix)  
    } else if (method=="keras") {
      method <- get("fit.keras.engine.training.Model", asNamespace("keras"))
      build_model <- methodArgs$build_model
      if (is.null(build_model)) build_model <- .default_build_model
      formalNames <- names(formals(method))
      methodArgs[setdiff(names(methodArgs),formalNames)] <- NULL ## so that only arguments of 'method' remain in methodArgs
      if ("x" %in% formalNames && is.null(methodArgs$x)) methodArgs$x <- as.matrix(totsim[randsample,stats,drop=FALSE]) 
      if ("y" %in% formalNames && is.null(methodArgs$y)) methodArgs$y <- totsim[randsample,x] 
      methodArgs$object <- do.call(build_model, args=list(input_shape=dim(methodArgs$x)[[2]]))
      if (is.null(methodArgs$epochs)) methodArgs$epochs <- 20L 
      abyss <- do.call("fit",methodArgs) # the LHS is not the projector
      projector <- methodArgs$object
      attr(projector,"stats") <- stats 
    } else {
      if (is.character(method)) {
        if (method=="ranger") { # availability already checked
          ## Given this is a Suggests'ed package, the necessary functions cannot be imported-from in the NAMESPACE  
          method <- get("ranger", asNamespace("ranger")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
        # } else if (method=="randomForest") { # availability already checked
        #   ## Given this is a Suggests'ed package, the necessary functions cannot be imported-from in the NAMESPACE  
        #   method <- get("randomForest.default", asNamespace("randomForest")) ## https://stackoverflow.com/questions/10022436/do-call-in-combination-with
        } else method <- eval(parse(text=method)) ## ranger must o through this
      }
      formalNames <- names(formals(method))  
      methodArgs[setdiff(names(methodArgs),formalNames)] <- NULL ## so that only arguments of 'method' remain in methodArgs
      if ("x" %in% formalNames && is.null(methodArgs$x)) methodArgs$x <- totsim[randsample,stats,drop=FALSE] 
      if ("y" %in% formalNames && is.null(methodArgs$y)) methodArgs$y <- totsim[randsample,x] ## drop is important... 
      if ("data" %in% formalNames && is.null(methodArgs$data)) methodArgs$data <- totsim[randsample,]
      if ("formula" %in% formalNames && is.null(methodArgs$formula)) methodArgs$formula <- 
        as.formula(paste(x," ~ ",paste(stats,collapse="+")))
      dotargs <- match.call(expand.dots = FALSE)$... ## produce a pairlist of (essentially) promises. No quote() needed
      for (st in names(dotargs)) methodArgs[[st]] <- dotargs[[st]]
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
        yy <- totsim[randsample,x] # but the calls could be used
        OKplot <- TRUE
      } else if (inherits(projector,"randomForest")) {
        yy <- totsim[randsample,x] ## (which is not returned in the randomForest object) # but the calls could be used
        OKplot <- TRUE
      } else if (inherits(projector,"keras.engine.training.Model")) {
        yy <- methodArgs$y
        OKplot <- TRUE
      }
      if (OKplot) { ## should try to plot oob.
        if (inherits(projector,"keras.engine.training.Model")) {
          plot(predict(projector, x=methodArgs$x),yy,xlab=paste("predicted",x),ylab=paste("true",x))
        } else if (inherits(projector,"ranger")) {
          plot(projector$predictions,yy,xlab=paste("oob-predicted",x),ylab=paste("true",x)) 
        } else if (inherits(projector,"randomForest")) {
          plot(predict(projector),yy,xlab=paste("oob-predicted",x),ylab=paste("true",x))
        } else plot(predict(projector,newdata=totsim[randsample,stats,drop=FALSE]),yy,xlab=paste("predicted",x),ylab=paste("true",x))
        abline(0,1)
      }
    }
    mc$data <- eval(mc$data)
    attr(projector,"project_call") <- mc
    return(projector) 
  }
})

# wrapper for handling names 
.predictWrap <- function(oneprojector,newdata, oob=Infusion.getOption("oob"), ...) {
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
  } ## else colnames should already be set
  if (inherits(oneprojector,"keras.engine.training.Model")) {
    if (is.data.frame(newdata)) {
      stats <- attr(oneprojector,"stats") 
      newdata <- newdata[,stats,drop=FALSE]
      newdata <- as.matrix(newdata)
    }
    predict(oneprojector,x=newdata,...)
  } else if (inherits(oneprojector,"ranger")) { 
    pred <- predict(oneprojector,data=newdata,...)$predictions 
    if (oob &&  ! is.null(dim(newdata))) {
      #message("oob used")
      x <- oneprojector$call$data[,oneprojector$forest$independent.variable.names] # data.frame
      dno <- dist(x,newdata[,colnames(x),drop=FALSE])
      newinold <- apply(dno==0L,2L,any)
      if (any(newinold)) {
        posinold <- apply(dno[,newinold,drop=FALSE]==0L,2L,function(v) which(v)[1L]) # in cas of proba masses which() will report several matches
        oobpred <- oneprojector$predictions
        pred[which(newinold)] <- oobpred[posinold] 
      } 
    }
    return(pred)
  } else if (inherits(oneprojector,"randomForest") && oob) { 
    pred <- predict(oneprojector,newdata=newdata,...) 
    if ( ! is.null(dim(newdata))) {
      x <- oneprojector$call$x
      dno <- dist(x,newdata[,colnames(x)])
      newinold <- apply(dno==0L,2L,any)
      if (any(newinold)) {
        posinold <- apply(dno[,newinold,drop=FALSE]==0L,2L,function(v) which(v)[1L]) # in cas of proba masses which() will report several matches
        oobpred <- predict(oneprojector,...) # ie oneprojector$predicted
        pred[which(newinold)] <- oobpred[posinold] 
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
project.default <- function (x,projectors,...) {
  #
  if (inherits(projectors,"list")) {
    projectors <- list2env(projectors)
  } else if (! inherits(projectors,"environment")) stop("'projectors' argument must be a environment or a list.")
  #
  if (inherits(x,"list")) { ##for list of EDFs where pars are not in the table 
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
  if (inherits(x,"EDF")) { ## tag assigned by this project.default function (see code if (inherits(x,"list"))...)
    # EDF => stats names are already there
    ly <- lapply(projectors,predict,newdata=as.data.frame(x)) 
  } else {
    #if ( ! is.null(stats)) {
    #  ly <- lapply(projectors,locpredict,newdata=x[,stats,drop=FALSE]) ## keep x cols unchanged for later use
    #} else 
    ly <- lapply(projectors, .predictWrap, newdata=x) # ordered as names(projectors), not ls(projectors)
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
  attr(ly,"packages") <- attr(x,"packages")
  attr(ly,"env") <- attr(x,"env") ## an environment!
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
