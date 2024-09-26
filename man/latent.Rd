\name{declare_latent}
\alias{declare_latent}
\alias{pplatent}
\alias{latint}
\title{
Modeling and predicting latent variables
}
\description{
Latent variables are unobserved variables which, for given model parameters, are random. Since they are unobserved, they cannot appear in data nor used to infer parameters. However, they can be predicted if their joint distribution with the data is learned from the reference table. Thus for inference about latent variables, these should be returned along summary statistics by the simulation function generating the samples for the reference table, but they should be declared as such so that later inference steps can distinguish them from both parameters and summary statistics. The \code{declare_latent} function is used for that purpose. 

The \code{pplatent} function can be used to point-predict latent values, by their inferred mean or median given the (projected) summary statistics and fitted parameter values. 

The \code{latint} function provides prediction intervals for the latent variable, accounting for uncertainty in parameter estimates, using a bootstrap method (see Details).     
}
\usage{
declare_latent(reftable, latentVars)
pplatent(object, type="mean",
         newDP=NULL,
         sumstats= t(get_from(object,"stat.obs")), 
         pars=t(object$MSL$MSLE), 
         size=1000L, ...)
latint(object, nsim=199L, levels=c(0.025,0.975), 
       sumstats= t(get_from(object,"stat.obs")), 
       Simulate, control.Simulate=get_from(object,"control.Simulate"), 
       bootSamples=NULL,
       ...)
}
\arguments{
  \item{reftable}{
A reference table of simulation as returned by \code{\link{add_reftable}}
}
  \item{latentVars}{
A vector of names of variables to be treated as latent variables.
}
\item{object}{
  An object of class \code{SLik_j}.
}
\item{type}{
  Character: the only handled non-default value is \code{"median"}.
}
\item{newDP, sumstats, pars}{
  Matrices of data and/or (\emph{projected}) summary statistics, with one column for each variable. \code{newDP} should contain both, and if NULL, is constructed from the next two arguments, \code{sumstats} holding statistics and \code{pars} holding parameters.
}
\item{levels}{
  Numeric vector: one-sided confidence levels (cumulative probabilities at which quantiles of the predictive distribution will be returned).
}
\item{nsim}{
  Integer: number of bootstrap replicates.
}
\item{size}{
  Integer: number of draws for estimating the median.
}
  \item{Simulate}{
  May be used to provide the simulation function if it is not stored in the \code{object}; usage is then as for \code{\link{add_reftable}}. If it is set to NULL, bootstrap samples will instead be drawn from the gaussian mixture fit, but this shoudl be avoided for good performance. When the argument is missing (default), the \code{Simulate} information is sought in the \code{object}, and if absent, the gaussian mixture fit is used with a warning. 
}
  \item{control.Simulate}{
    A list of arguments of the \code{Simulate} function (see\code{\link{add_simulation}}). The default value should be used unless you understand enough of its structure to modify it wisely (e.g., it may contain the path of an executable used to perform the fit on one machine and a different path may be specified to compute the prediction interval on another machine).  
  }
  \item{bootSamples}{
A data frame. The bootstrap samples may be provided by this argument, otherwise they will be automatically simulated by the function provided by the \code{Simulate} argument.  The boot samples may for example be obtained by calling \code{simulate(}<\emph{Slik_j object}>\code{, SGP=TRUE)}. 
  }
  \item{\ldots}{
  For \code{pplatent}: Not currently used. For \code{latint}: may be used to pass arguments \code{verbose}, \code{nb_cores}, \code{packages}, \code{env}, \code{control.Simulate}, \code{cluster_args}, and \code{cl_seed}, with usage as described for \code{\link{add_reftable}}, to the bootstrap sample simulation step. Parallelisation controls will also be used for the other steps of the bootstrap procedure.
}
}
\details{
The \code{latint} function aims to provide intervals for the latent variable \eqn{V}, with controlled coverage for given parameter values. The parametric bootstrap method of Lawless & Fredette (2005) can be adapted to the present problem. First, new data \eqn{D^*} and new latent values \eqn{V^*} are jointly simulated, given the summary-ML estimates (this uses the sample-generating process simulator, but these simulations are not added to the reference table). Second, their original method requires the evaluation for each new \eqn{(D^*,V^*)} of \eqn{F_V(V=V^*|D^*;\theta(D^*))}, the value of the cumulative distribution function \eqn{F_V} of \eqn{V} evaluated at \eqn{V=V^*} values, given \eqn{D^*} and given new parameter estimates \eqn{\theta(D^*)}. % Their eq 14 in particular.
It also requires the quantile function of \eqn{V} for the original data and parameter estimates. \code{latint} uses the fit of the multivariate gaussian mixture model to the reference table, stored in the fit object, for fast parameter refitting and for fast estimation of these functions on each bootstrap sample \eqn{(D^*,V^*)}. 

\code{latint(., levels=c(0.5))} will return the median of the predictive distribution, 
conceptually distinct from the plug-in prediction by \code{pplatent(slik_j, type="median")}   

Quantile computations in \code{pplatent} and \code{latint} are approximate and might be modified in some future version, but were sufficient to obtain reasonable results in simulations. 
}
\value{
\code{declare_latent} returns the input \code{reftable} with modified attributes. \code{pplatent} returns a vector of predicted latent values.

\code{pplatent} returns a single numeric value or a vector.
\code{latint} returns a list of matrices. Each matrix has with one column per element of \code{levels}, and one row per row of \code{sumstats}. There is one such matrix for each latent variable.  
}
\references{
Lawless J. F., Fredette M. (2005) Frequentist prediction intervals and predictive distributions. Biometrika 92: 529–542, <doi:10.1093/biomet/92.3.529>
}
\examples{
\dontrun{
  ##### A toy example motivated by some inference problem for genomic data #####
  # A model with three parameters logNe, Vs and Es is fitted.
  # Data from 100 loci are here summarized by three genome-wide summary statistics 
  # (slogNe, sVs and sEs), and one locus-specific statistic that will provide 
  # information about a locus-specific latent variable.
  
  ## Simulation function 
  genomloc <- function(logNe=parvec["logNe"],Es=parvec["Es"],Vs=parvec["Vs"], 
                       latent=TRUE, # returns the latent value by default
                       parvec) {
    slogNe <- rnorm(1,logNe, sd=0.3)
    genom_s <- rgamma(99, shape=Es/Vs,scale=Vs) # all loci except the focal one
    sEs <- mean(genom_s)
    sVs <- var(genom_s)
    latentv <- rgamma(1, shape=Es/Vs,scale=Vs) # locus-specific latent variable to predict
    sloc <- rnorm(1, mean=latentv-sEs,sd=latentv/3) # locus-specific statistic
    resu <- list(slogNe=slogNe,sEs=sEs,sVs=sVs, sloc=sloc)
    if (latent) resu$latentv <- latentv
    unlist(resu)
  } 
  #
  ## simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- genomloc(logNe=4,Es=0.05, Vs=0.1,latent=FALSE) ## no latent value
  #
  workflow_design <- get_workflow_design(npar=3L, n_proj_stats=4L, n_latent=1L)
  parsp <- init_reftable(lower=c(logNe=2,Es=0.001,Vs=0.001), 
                         upper=c(logNe=6,Es=0.2,Vs=0.2),
                         nUnique=workflow_design$init_reft_size)
  simuls <- add_reftable(Simulate=genomloc, parsTable=parsp)
  
  simuls <- declare_latent(simuls,"latentv") 

  ## Projections are not necessary here since the number of statistics is minimal,
  # but will be discussed later.
  { ############ Without projections
    { ## Usual workflow for estimation: 
      densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
      slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
      slik_j <- refine(slik_j,maxit=2,update_projectors=TRUE)
      # plot1Dprof(slik_j) ## 1D profiles show parameter inference is OK
    }
    { ## inference about latent values:
      pplatent(slik_j)
      pplatent(slik_j, type="median")
      latint(slik_j, nsim=999, levels=c(0.025,0.5,0.975))
    }
    { ## Assessing prediction of latent variable:
      # Builds testing set:
      test_simuls <- t(replicate(1000, genomloc(logNe=4,Es=0.05, Vs=0.1)))
      test_data <- test_simuls[,-5]
      # Point prediction:
      pred <- pplatent(slik_j, sumstats = test_data)
      
      plot(test_simuls[,"latentv"], pred); abline(0,1) # prediction vs. true latent values
    }
  }
  
  { ########## Beyond standard use of projections for estimation of parameter values, 
    # projections can also be used when several individual-level statistics inform about 
    # the latent variable, to reduce them to a single summary statistic.
    # Projection will then be needed at the prediction step too.
    
    { # projection with latent variable as response:
      platent <- (project("latentv", data=simuls, stats=c("slogNe","sEs","sVs","sloc")))
      # (This example only serves to show the syntax since no dimention reduction occurs)
      
      dprojectors <- list(SLOC=platent,slogNe=NULL,sEs=NULL, sVs=NULL)
      
      # => As soon as one projection is used, The 'projectors' argument must include 
      # all projectors used for the inference, whether for parameters or for latent variables. 
      # NULL projectors should then be declared for raw statistics retained 
      # in the projected reference table.
      
      # Apply projections on simulated statistics and 'data':
      projSimuls <- project(simuls,projectors=dprojectors,verbose=FALSE)
      projSobs <- project(Sobs,projectors=dprojectors)
    }
    
    { ## Estimation: 
      ddensv <- infer_SLik_joint(projSimuls,stat.obs=projSobs)
      dslik_j <- MSL(ddensv) ## find the maximum of the log-likelihood surface
      dslik_j <- refine(dslik_j,maxit=2,update_projectors=TRUE)
      # plot1Dprof(dslik_j)
    }
    
    { ## Assessing prediction of latent variable: do not forget to project!
      
      test_simuls <- t(replicate(1000, genomloc(logNe=4,Es=0.05, Vs=0.1)))
      test_data <- test_simuls[,-5] # removing column of latent variable
      ptest_data <- project(test_data,projectors=dprojectors,verbose=FALSE) # Here!
      pred <- pplatent(dslik_j, sumstats = ptest_data)
      
      plot(test_simuls[,"latentv"], pred); abline(0,1)
    }
  }

}
}
