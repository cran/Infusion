# Used only for the never-operational smoothing by Poisson GLMM
.seekSimplex <- function(point,constraintsList) {
  it <- 1L
  while(it < length(constraintsList)+1L) {
    if (isPointInCHull(point,constraints=constraintsList[[it]])) {
      return(it)
    } else it <- it+1L
  }
  return(NA) ## failed to locate the point in the tesselation
}

## wrapper for rsimplex, works on several focal points, handles case where focal is out of triangulation
#jitterPoints <- function(focalpoints,vT,constraintsList,expand=1,u=NULL) { ## FR->FR not used removed 04/2016
#...

.findConstraints <- function(vT) {
  sT <- vT$simplicesTable
  knots <- vT$vertices
  constraintsList <- lapply(seq(nrow(sT)),function(it) {
    resetCHull(knots[sT[it,],,drop=FALSE],formats=c("constraints"))[c("a","b")]
  }) ## bizarrement lent... => saved as attribute
  return(constraintsList)
}

## Used only for the never-operational smoothing by Poisson GLMM (=> only for the primitive workflow)
#
## * Constructs a binning 
##   using the vertices of the convex hull + a random sample of the points to be binned
##   this grid is further "expanded" 
##   FR->FR rather use rhull ?
##   or use first volTr
## * First tests if focal point is in convex envelope before perfoming all these complex operations. 
multi_binning <- function(m,subsize=trunc(nrow(m)^(Infusion.getOption("binningExponent"))),expand=5/100,focal=NULL) {
  m <- m[,names(focal),drop=FALSE] ## convenient way of selection cols
  if (ncol(m)==0L) stop("'focal' elements must have names, matching colnames of 'm'")
  chull <- resetCHull(m,formats=c("vertices","constraints")) ## full triangulation also covers the *convex* hull
  ## test whether focal point is in convex hull
  if ( ! is.null(focal) && ! isPointInCHull(focal,constraints=chull[c("a","b")])) {
    return(NULL)
  } ## ELSE
  center <- colMeans(m)
  knots <- m[sample(nrow(m),subsize),,drop=FALSE] ## take random points
  envelope <- chull$vertices ## exterior vertices = convex hull
  knots <- as.matrix(rbind(knots,envelope)) ## add convex envelope
  dif <- apply(knots,1,'-',center)
  if (is.null(dim(dif))) {
    dim(dif) <- c(length(dif),1)
  } else dif <- t(dif) 
  knots <- knots + dif*expand ## expand## FR->FR adapt to the average relative volume of the cells ?
  dif <- apply(envelope,1,'-',center)
  if (is.null(dim(dif))) {
    dim(dif) <- c(length(dif),1)
  } else dif <- t(dif) 
  expandEnv <- envelope + dif*Infusion.getOption("zeromargin") ## expand more
  knots <- as.matrix(rbind(knots,expandEnv)) ## add convex envelope
  knots <-rbind(knots,focal) ## to increase precision here
  # pb when focal is near the edge... and ugly solution since findConstraints is slow
  ## add vertices of simplex where focal is ! $vertices not ordered as m !
  #   vT <- volTriangulation(m)  
  #   constraintsList <- .findConstraints(vT)
  #   whereIsStatObs <- .seekSimplex(focal,constraintsList)
  #   knots <- rbind(knots,vT$vertices[vT$simplicesTable[whereIsStatObs,],]) 
  #
  vT <- volTriangulation(as.matrix(knots)) ## j'ai rajouté les bary dans la sortie. Y penser pour migraineLike 
  ## locate points from m into the simplices
  constraintsList <- .findConstraints(vT) ## FR->FR slow
  foundSimplices <- apply(m,1, .seekSimplex, constraintsList=constraintsList)
  messy <- table(foundSimplices) ## zero counts restored below  (while 'tabulate' keep them only for internal categories !)
  counts <- integer(length(constraintsList))
  counts[as.numeric(names(messy))] <- messy ## restores zero counts
  volumes <- vT$vol
  barys <- vT$bary
  sumCounts <- sum(counts)
  resu <- data.frame(cbind(barys,binFactor=volumes*sumCounts,count=counts))   
  attr(resu,"stats") <- colnames(barys)
  attr(resu,"counts") <- "count"
  attr(resu,"binFactor") <- "binFactor"
  resu
}

.smoothEDF_s <- function(data,pars=c(), verbose=FALSE,...) {
  #browser()
  stats <- attr(data,"stats")
  counts <- attr(data,"counts")
  
  if (length(pars)>0) { ## pas vraiment implementé
    form <- paste(counts,"~ 0 + Matern(1|",paste(c(pars,stats),collapse=" + "),")")    
  } else form <- paste(counts,"~ 0 + Matern(1|",paste(stats,collapse="+"),")")
  if (is.null(data$trend)) {
    form <- as.formula(paste(form,"+offset(log(",attr(data,"binFactor"),"))")) ## 
  } else form <- as.formula(paste(form,"+offset(log(",attr(data,"binFactor"),")+log(trend))")) ## 
  init <- list(rho=rep(1,length(stats)+length(pars)))
  arglist <- list(formula=form,init=init,
                  HLmethod=.Infusion.data$options$HLmethod,family=poisson(),data=data,fixed=list(nu=4)) ## HLmethod makes a difference at low response
  ## non-default control of fitme
  dotlist <- list(...)
  formalsNames <- names(formals(fitme))
  controlNames <- intersect(names(dotlist),formalsNames)
  ## fixed must be treated differently from other formalsNames to keep default fixed$nu if no replacement
  controlNames <- setdiff(controlNames,"fixed") ##  
  ## same idea for init
  controlNames <- setdiff(controlNames,"init")  
  arglist[controlNames] <- dotlist[controlNames] 
  arglist$init[names(dotlist$init)] <- dotlist$init
  arglist$fixed[names(dotlist$fixed)] <- dotlist$fixed
  arglist$init[names(dotlist$fixed)] <- NULL
  ## tentative subsampling
  nr <- nrow(data)
  ns <- max(20,nr^(2/3))
  subidx <- sample(nr,ns)
  arglist$data <- data[subidx,]
  fit <- do.call(fitme,arglist)
  #if(verbose) cat("Smoothing the data, may be slow...\n")
  arglist$data <- data
  corrPars1 <- get_ranPars(fit,which="corrPars")[["1"]]
  ranfix <- c(corrPars1,list(lambda=fit$lambda))
  arglist$fixed <- ranfix
  arglist$init <- NULL
  fit <- do.call(fitme,arglist)
  ## to see the fit, redata <- data; redata$binFactor <- 1; predict(fit,newdata=redata)
  #browser()
  return(fit)  
}




## the cases where there are replicates, one with and one without a logL estimate, are among the misleading ones dealt here:  
.remove.pairswithNas <- function(logLs) {
  fittedPars <- attr(logLs,"colTypes")$fittedPars 
  logLname <- attr(logLs,"colTypes")$logLname
  ordered.logLs <- logLs[do.call(order,logLs),]
  y <- ordered.logLs[,logLname]
  diffs <- diff(as.matrix(ordered.logLs[,fittedPars,drop=FALSE]))
  maxdiffs <- apply(diffs,1,function(v) {max(abs(v))})
  # One need to makesure if two successive pairs P1, P2, the second row of P1 having NA, the P2 pair is not declared misleading 
  # hence,it is incorrect to [test and expand by c(FALSE...) | ...] separately the maxdiffs and the NAs
  #  as expension of the test on NA's only would yield T,T,T,F
  #  and expension of the test on pairs only would yield T,T,T,T => joint test would exclude the first three rows!
  isNAy <- is.na(y)
  misleading <- maxdiffs==0 & (isNAy[-length(isNAy)] | isNAy[-1]) ## replicates avec au moins un des deux ayant NA
  # in above example, this is (T,F,T) & (T,T,F) = (T,F,F) and expansion is (T,T,F,F)
  misleading <- c(FALSE,misleading) | c(misleading,FALSE) ## tous les paires des replicats dont un membre est comme ci dessus
  ## this leaves NA's not in pairs in the data ! : 
  ordered.logLs[ ! misleading,,drop=FALSE]
}

.prettysignif <- function(x,extradigits=0) {
  ## extradigits=0 => 99(99)99->99(99)99; 999.9->999.9; 9.999->9.999; 0.xxx -> 3 chifres signif
  for(i in seq_len(length(x))) {
    xi <- x[i]
    if (is.na(xi)) x[i] <- NA
    #ELSE
    if (xi<1) {n<-3} else if (xi>1000) {n <- ceiling(log(xi,10))} else {n <- 4}
    n <- n+extradigits
    x[i] <- signif(xi,n)
  }
  x
}

# new variable name
.makenewname <- function(base,varnames) { 
  pattern <- paste(base,"[0-9]+",sep="")
  allmatches <- dir(pattern=pattern)
  allremainders <- substring(allmatches, nchar(base)+1L) 
  allnumericremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  ## 2015/03/04
  if (length(allremainders) == 0L && length(allnumericremainders) == 0L) { ## if base = allremainders => length(allnumericremainders) == 0 not sufficient
    fvname <- base
  } else {
    if (length(allnumericremainders)) {
      num <- max(allnumericremainders)+1L
    } else num <- 1L
    fvname <- paste ( base , num , sep="") 
  }
  fvname
}

.generateFileName <- function(base="tmp",ext="") { ## for a file
  pattern <- paste(base,".*",ext,sep="") # regexp .* means any character, zero or more time
  allmatches <- dir(pattern=pattern)
  allremainders <- substring(allmatches,nchar(base)+1)
  allremainders <- unlist(strsplit(allremainders,ext)) ## removes the extension from the remainder 
  allremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  if (length(allremainders) == 0) {
    num <- 0
  } else num <- max(allremainders)+1
  validFileName <-paste ( base , num , ext,sep="") 
  return(validFileName)
}

.lookup <- function(arglist, which=NULL, try_in =NULL) {
  if ( ! is.null(try_in)) {
    resu <- .lookup(arglist[[try_in]], which=which, try_in=NULL)
    if ( ! is.null(resu)) return(resu)
  }
  if ( ! is.null(which)) {
    arglist[[which]]
  } else arglist
}

# extractor handling SLik/Slik_j difference and changes in Slik_j structure (from v2.1.28 onwards)
.get_reft_raw <- function(object) {
  if (inherits(object,"SLik")) {
    rr <- attr(object$logLs,"raw_data")
  } else if (is.null(rr <- object$reftable_raw)) rr <- object$raw_data
  #  object$reftable_raw should be null when there are no projector, and rr <- object$raw_data only for back compat
  rr
}


