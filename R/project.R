project <- function(x, ...) UseMethod("project")


project.character <- function(x,
                         stats,
                         data,
                         trainingsize= if (method=="REML") {Infusion.getOption("projTrainingSize")} else {NULL},
                         knotnbr= if (method %in% c("REML","GCV")) {
                           Infusion.getOption("projKnotNbr")
                           } else {floor(1000*log2(length(stats)+1))},
                         method="REML",methodArgs=list(),
                         verbose=TRUE,
                         ...) {
  if (!is.character(method)) {stop("'method' should be a character string, \n in particular a function name rather than a function")}
  if (is.data.frame(data)) {
    totsim <- data
  } else if (inherits(data,"list")) {
    pardata <- lapply(data,function(ll) {
      cbind(ll,attr(ll,"par")[x],row.names=NULL)
    }) ## puts the parameter values into the data
    totsim <- do.call(rbind,pardata)
  } else {
    stop("'data' are neither a data frame nor a list of data frames.")
  }
  nr <- nrow(totsim)
  if (is.character(method) && method %in% c("GCV","REML")) {
    form <- as.formula(paste(x," ~ 1 + Matern(1|",paste(stats,collapse="+"),")"))
    # utiliser totsim dans ce qui suit => plantage severe
    if (method=="REML") {
      randsim <- totsim[sample(nr,min(nr,trainingsize)),]
      ## hummmm prend ~ 3 minutes pour 300 points
      ## pas de replicat de X available for projection -> estim phi/lambda pb et pas de GCV possible...
      if (eval(Infusion.getOption("fitmeCondition"))) {
        essfit <- fitme(form,data=randsim, fixed=list(nu=4), 
                         method=method, init=list(rho=rep(NA,length(stats)))) 
      } else essfit <- corrHLfit(form,data=randsim,init.corrHLfit=list(rho=rep(NA,length(stats))),ranFix=list(nu=4))  
      ranfix <- c(essfit$corrPars,list(lambda=essfit$lambda,phi=essfit$phi))
      randsim <- totsim[sample(nrow(totsim),min(nr,knotnbr)),]
    } else if (method=="GCV") {## no difference between trainingSize and knotNbr
      randsim <- totsim[sample(nrow(totsim),min(nr,knotnbr)),]
      forGCV <- prepareData(data=randsim,ParameterNames=stats,
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
      projector <- fitme(form,data=randsim, fixed=ranfix, method="REML") 
    } else projector <- corrHLfit(form,data=randsim,ranFix=ranfix)  
  } else {
    randsim <- totsim[sample(nrow(totsim),min(nr,knotnbr)),]
    if (is.character(method)) method <- eval(parse(text=method))
    formalNames <- names(formals(method))  
    methodArgs[setdiff(names(methodArgs),formalNames)] <- NULL ## so that only arguments of 'method' remain in methodArgs
    if ("x" %in% formalNames && is.null(methodArgs$x)) methodArgs$x <- randsim[,stats,drop=FALSE] 
    if ("y" %in% formalNames && is.null(methodArgs$y)) methodArgs$y <- randsim[,x] ## drop is important... 
    if ("data" %in% formalNames && is.null(methodArgs$data)) methodArgs$data <- randsim
    if ("formula" %in% formalNames && is.null(methodArgs$formula)) methodArgs$formula <- 
      as.formula(paste(x," ~ ",paste(stats,collapse="+")))
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
    } else if (inherits(projector,"randomForest")) {
      yy <- randsim[,x] ## (which is not returned in the randomForest object)
      OKplot <- TRUE
    }
    if (OKplot) {
      plot(predict(projector),yy,xlab=paste("predicted",x),ylab=paste("true",x))
      abline(0,1)
    }
  }
  return(projector) ## note the two returns
}

# wrapper for handling names 
predictWrap <- function(oneprojector,newdata,...) {
  if (inherits(newdata,"numeric")) {
    stats <- attr(oneprojector,"stats") ## should be non-null if projector.character used a formula
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
  predict(oneprojector,newdata=newdata,...)
}

# x attributes used late in code, hence x should not be modified
project.default <- function (x,projectors,...) {
  #
  if (! inherits(projectors,"list")) stop("'projectors' argument must be a list")
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
    attr(ly,"projectors") <- as.list(match.call())$projectors
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
    ly <- lapply(projectors,predictWrap,newdata=x)    
  }
  ly <- do.call(cbind,ly) ## binding is over projectors
  if (is.vector(x)) {
    ly <- as.numeric(ly)
    names(ly) <- names(projectors)
  } else { ## expects ly to be a matrix
    colnames(ly) <- names(projectors)
    if ( ! is.null(parms <- attr(x,"parms"))) { ## assumes that the new data 'x' has not be modified
      ly <- cbind(x[,parms,drop=FALSE],ly)
      attr(ly,"parms") <- parms
      attr(ly,"stats") <- names(projectors)
    }
  }
  attr(ly,"Simulate") <- attr(x,"Simulate")
  attr(ly,"projectors") <- as.list(match.call())$projectors
  ly ## basically a matrix
}

neuralNet <- function(formula,data) {
  if (requireNamespace("caret",quietly=TRUE)) {
    caret::train(formula, data, method='nnet',
          linout=TRUE,
          trace = FALSE)
  } else {
    stop("'caret' package not available.")
  }
}