reparam_reftable <- function(fitobject, 
                             to, # all pars...
                             reparamfn, # two arguments: 1st arg is a reftable, second is 'to' or ... 
                             LOWER=NULL, UPPER=NULL,
                             raw=FALSE, 
                             reftable_attrs=NULL # final control of attributes
                             ) { # may contain only specifi values) {
  # I shouldwork on fixed+fittedPars...
  allPars <- c(fitobject$colTypes$fittedPars,
               names(fitobject$colTypes$fixedPars))
  if (raw) {
    oldreft <- get_from(fitobject,"reftable_raw") 
  } else {
    oldreft <- get_from(fitobject,"reftable")
    allPars <- intersect(allPars, colnames(oldreft))
  }   
  np <- length(allPars)
  if (length(to)!=np) stop("length(to) != length(allPars)")
  mostAttrs <- names(attributes(oldreft))
  mostAttrs <- setdiff(mostAttrs,c("class","dim","dimnames","names","row.names",
                                   "allPars","LOWER","UPPER"))
  reft <- reparamfn(oldreft[,allPars], to=to)
  reft <- cbind(reft, oldreft[,setdiff(colnames(oldreft), allPars)])
  for (attrname in mostAttrs) attr(reft,attrname) <- attr(oldreft,attrname)
  attr(reft,"allPars") <- to # infer_SLik_joint expects this. colTypes will be re-evaluated by it.
  oldLOWER <- fitobject$LOWER
  oldUPPER <- fitobject$UPPER
  newUPPER <- newLOWER <- setNames(rep(NA_real_,np),to)
  for (st in to) {
    lowval <- oldLOWER[st]
    if (is.null(lowval) || is.na(lowval)) lowval <- LOWER[st]
    if (is.null(lowval) || is.na(lowval)) lowval <- min(reft[,st])
    newLOWER[st] <- lowval
    
    upval <- oldUPPER[st]
    if (is.null(upval) || is.na(upval)) upval <- UPPER[st]
    if (is.null(upval) || is.na(upval)) upval <- max(reft[,st])
    newUPPER[st] <- upval
  }
  attr(reft,"LOWER") <- newLOWER[to]
  attr(reft,"UPPER") <- newUPPER[to]

  for (st in names(reftable_attrs)) attr(reft,st) <- reftable_attrs[[st]]

  class(reft) <- c("reftable",class(reft))
  reft
} 

reparam_fit <- function(fitobject, 
                         to, 
                         reparamfn, # two arguments: 1st arg is a reftable, second is 'to' or ... 
                         LOWER=NULL, UPPER=NULL, # may contain only specific values
                         nbCluster="max",
                         constr_crits=get_from(fitobject,"constr_crits"),
                        raw=FALSE, 
                        #statNames=names(get_from(fitobject,"raw_data")),
                        reftable_attrs=NULL,
                         ...) {
  reft <- reparam_reftable(fitobject=fitobject, 
                           to=to, 
                           reparamfn=reparamfn, raw=raw,
                           LOWER=LOWER, UPPER=UPPER, reftable_attrs=reftable_attrs)
  if ( ! is.null(constr_crits)) {
    ## Checking syntactic correctness only
    cat("Checking feasability of 'constr_crits' on transformed parameter space... ")
    eval(constr_crits, envir=as.list(reft[1, ]))
    cat("OK")
  }
  if (raw) {
    isVar_Pars <- apply(reft[,to,drop=FALSE], 2L, function(v) length(unique(range(v)))>1L)
    fittedPars <- names(which(isVar_Pars))
    pfittedPars <- paste0("p",fittedPars)
    statNames <- names(get_from(fitobject,"raw_data"))
    projectors <- setNames(vector("list", length(fittedPars)), pfittedPars)
    for (it in seq_along(fittedPars)) {
      projectors[[pfittedPars[it]]] <- project(fittedPars[it],stats=statNames,
                                               data=reft,
                                               verbose= FALSE)
    }
    reft <- project(reft, projectors = projectors,is_trainset=TRUE, 
                    use_oob=TRUE)
    pSobs <- project(get_from(fitobject,"raw_data"), projectors = projectors,
                     is_trainset=FALSE, use_oob=FALSE)
    refit <- infer_SLik_joint(reft,
                              stat.obs = pSobs,
                              nbCluster=nbCluster,
                              constr_crits=constr_crits)
  } else  {
    refit <- infer_SLik_joint(reft, 
                              stat.obs = attr(reft,"stat.obs"), 
                              nbCluster=nbCluster,
                              constr_crits=constr_crits)
  }
  MSL(refit, ...)
}

