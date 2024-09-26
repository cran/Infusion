\name{reparam_fit}
\alias{reparam_fit}
\alias{reparam_reftable}
\title{
Conversion to new parameter spaces
}
\description{
Functions to facilitate inferences using alternative parametrizations of the same model, reusing an existing reference table. \code{reparam_reftable} produces a reference table in the new parametrization.  \code{reparam_fit} does the same internally, and runs \code{infer_SLik_joint} and \code{MSL} on the new reference table. 

\code{reparam_fit} is \bold{experimental} and may have various limitations. %To circumvent such a limitation, more extended coding is required, as illustrated in the \code{test-reparam} example of the \code{testthat} subdirectory of the package sources.
}
\usage{
reparam_fit(fitobject, to, reparamfn, 
             LOWER=NULL, UPPER=NULL, nbCluster="max", 
             constr_crits = get_from(fitobject, "constr_crits"),
             raw=FALSE, 
             %statNames=names(get_from(fitobject,"raw_data")),
             reftable_attrs=NULL,
             ...)
reparam_reftable(fitobject, to, reparamfn, 
             LOWER=NULL, UPPER=NULL, raw=FALSE, reftable_attrs=NULL)
}
\arguments{
  \item{fitobject}{
an object of class \code{SLik_j}.
}
  \item{to}{
Character vector: names of all parameters in the new parametrization.
}
  \item{reparamfn}{
A function which can convert a data frame/matrix/vector of \dQuote{old} parameters to an object of the same class in the new parametrization. It must have two arguments: its first argument must hold the \dQuote{old} data frame/matrix/vector; and the second argument is either \code{to} (holding the \code{to} argument given to \code{reparam_reftable}) if \code{reparamfn} needs this information, or \code{\ldots}.
}
  \item{LOWER, UPPER}{
Optional named vectors of bounds. They may be incomplete, containing values only for new parameters, and not necessarily for all of them (missing information is deduced from the observed ranges in the reference table).  
}
  \item{nbCluster}{
Passed to \code{\link{infer_SLik_joint}}.
}
  \item{constr_crits}{
\code{\link{constr_crits}} applicable in the new parametrization. The suitability of such constraints is checked on the transformed reference table. When this argument is ignored, its default value is taken from the input object and therefore refers to the old parametrization. The check may then highlight the need for providing constraints redefined in reference to the new parametrization.   
}
  \item{raw}{Boolean; if TRUE, the object is re-built starting from the raw 
  reference table. In particular the projections are re-computed.
}
\item{reftable_attrs}{
A \code{list} whose elements are set as attributes to the re-parametrized reference table (see Value in \code{\link{add_reftable}}). Elements not provided by this argument will be copied from the input reference table. A typical use is to provide as \code{Simulate} function in the new parametrization (see Examples). 
}
\item{\ldots}{
Passed to \code{\link{MSL}}.
}
}
\value{
\code{reparam_reftable} returns a reference table with attributes, suitable as input for \code{\link{infer_SLik_joint}}.
\code{reparam_fit} returns the return value of an \code{\link{MSL}} call.

The information about projections retained in these objects come from original \code{fitobject}.  
}
\examples{
\dontrun{

## Toy simulation function 
# (inspired by elementary population-genetic scenario)

hezsim <- function(logNe=parvec["logNe"],
                   logmu=parvec["logmu"],parvec) {
  Ne <- 10^logNe
  mu <- 10^logmu
  Es <- Ne*mu
  Vs <- 1/log(1+Ne) 
  genom_s <- rgamma(5, shape=Es/Vs,scale=Vs) # 5 summary statistics
  names(genom_s) <- paste0("stat",seq(5))
  genom_s
} 


{ ## Analysis with 'canonical' parameters
  #
  ## simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- hezsim(logNe=4,logmu=-4) 
  #
  parsp <- init_reftable(lower=c(logNe=1,logmu=-5), 
                         upper=c(logNe=6,logmu=-2))
  init_reft_size <- nrow(parsp)
  simuls <- add_reftable(Simulate=hezsim, parsTable=parsp)
  
  { 
    plogNe <- project("logNe", data=simuls, stats=paste0("stat",seq(5)))
    plogmu <- project("logmu", data=simuls, stats=paste0("stat",seq(5)))

    dprojectors <- list(plogNe=plogNe,plogmu=plogmu)
    
    projSimuls <- project(simuls,projectors=dprojectors,verbose=FALSE)
    projSobs <- project(Sobs,projectors=dprojectors)
  }
  
  { ## Estimation: 
    ddensv <- infer_SLik_joint(projSimuls,stat.obs=projSobs)
    dslik_j <- MSL(ddensv, eval_RMSEs=FALSE) ## find the maximum of the log-likelihood surface
    refined_dslik_j <- refine(dslik_j, eval_RMSEs=FALSE, CIs=FALSE)
  }
}

{  ## Reparametrization to composite parameters

  locreparamfn <- function(object, ...) {
    logTh <- object[["logmu"]]+object[["logNe"]]
    if (inherits(object,"data.frame")) { # *data.frame case always needed.*
      data.frame(logTh=logTh,
                 logNe=object[["logNe"]])
    } else if (is.matrix(object)) {
      cbind(logTh=logTh,
            logNe=object[["logNe"]])
    } else c(logTh=logTh,
             logNe=object[["logNe"]])
  }
  
  { ## without re-projection
     rps <- reparam_fit(refined_dslik_j, to=c("logTh","logNe"),
                      reparamfn = locreparamfn)
     plot(rps)
  }

  
  { ## with re-projection [necessary to allow refine()'s]
  
    # For refine() a new simulation will be needed, with new input parameters: 
    hezsim2 <- function(logNe=parvec["logNe"],logTh=parvec["logTh"],parvec) {
      hezsim(logNe=logNe,logmu=logTh-logNe)
    } 
    
    rps <- reparam_fit(refined_dslik_j, to=c("logTh","logNe"),
                       reparamfn = locreparamfn, 
                       raw=TRUE,                   # to allow re-projection
                       reftable_attrs=list(Simulate=hezsim2))
    plot(rps)
    refine(rps)

  }
  
}
}
}
