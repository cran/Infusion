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
  init.corrHLfit <- list(rho=rep(1,length(stats)+length(pars)))
  arglist <- list(formula=form,init.corrHLfit=init.corrHLfit,
                  HLmethod=.Infusion.data$options$HLmethod,family=poisson(),data=data,ranFix=list(nu=4)) ## HLmethod makes a difference at low response
  ## non-default control of corrHLfit
  dotlist <- list(...)
  formalsNames <- names(formals(corrHLfit))
  controlNames <- intersect(names(dotlist),formalsNames)
  ## ranFix must be trated differently from other formalsNames to keep default ranFix$nu if no replacement
  controlNames <- setdiff(controlNames,"ranFix") ##  
  ## same idea for init.corrHLfit
  controlNames <- setdiff(controlNames,"init.corrHLfit")  
  arglist[controlNames] <- dotlist[controlNames] 
  arglist$init.corrHLfit[names(dotlist$init.corrHLfit)] <- dotlist$init.corrHLfit
  arglist$ranFix[names(dotlist$ranFix)] <- dotlist$ranFix
  arglist$init.corrHLfit[names(dotlist$ranFix)] <- NULL
  ## tentative subsampling
  nr <- nrow(data)
  ns <- max(20,nr^(2/3))
  subidx <- sample(nr,ns)
  arglist$data <- data[subidx,]
  fit <- do.call(corrHLfit,arglist)
  #if(verbose) cat("Smoothing the data, may be slow...\n")
  arglist$data <- data
  corrPars <- fit$corrPars[["1"]]
  if (is.null(corrPars)) corrPars <- fit$corrPars ## F I X M E transitional code 
  ranfix <- c(corrPars,list(lambda=fit$lambda))
  arglist$ranFix <- ranfix
  arglist$`init.corrHLfit` <- NULL
  fit <- do.call(corrHLfit,arglist)
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
  if (is.na(x)) return(NA)
  #ELSE
  if (x<1) {n<-3} else if (x>1000) {n <- ceiling(log(x,10))} else {n <- 4}
  n <- n+extradigits
  signif(x,n)
}

.allCIs <- function(object,level=0.95, verbose=TRUE,method="LR",...) {
  CIs <- list()
  for (st in object$colTypes$fittedPars) {CIs[[st]] <- confint(object,st,level=level, verbose=verbose,method=method,...)}
  lowers <- lapply(CIs,function(li) {li$lowerpar})
  if ( ! is.null(lowers)) names(lowers) <-  paste("low.",names(lowers),sep="")
  uppers <- lapply(CIs,function(li) {li$upperpar})
  if ( ! is.null(uppers)) names(uppers) <-  paste("up.",names(uppers),sep="")
  # the list elements are either numeric vectors or asingle NA... 
  ordre <- order(c(seq_len(length(lowers)),seq_len(length(uppers))))
  bounds <- c(lowers,uppers)[ordre]
  checkvec <- function(vec) {if (identical(vec,NA)) {return(NULL)} else {return(vec)} }
  whichNAs <- unlist(lapply(bounds,function(vec) {identical(vec,NA)}))
  missingBounds <- names(bounds[whichNAs])
  bounds <- do.call(rbind,bounds[ ! whichNAs])
  return(list(CIs=CIs,bounds=bounds,missingBounds=missingBounds,level=level)) 
}

.makenewname <- function(base,varnames) { ## post CRAN 1.4.1
  varnames <- varnames[which(substring(varnames,1,nchar(base))==base)] 
  allremainders <- substring(varnames,nchar(base)+1) 
  allnumericremainders <- as.numeric(allremainders[which( ! is.na(as.numeric(allremainders )))  ]) ## as.numeric("...")
  ## 2015/03/04
  if (length(allremainders) == 0L && length(allnumericremainders) == 0L) { ## if base = allremainders => length(allnumericremainders) == 0 not sufficient
    fvname <- base
  } else {
    if (length(allnumericremainders)) {
      num <- max(allnumericremainders)+1L
    } else num <- 1L
    fvname <-paste ( base , num , sep="") 
  }
  fvname
}
