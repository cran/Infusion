\name{options}
\alias{Infusion.options}
\alias{Infusion.getOption}
\alias{parallel}
\title{Infusion options settings}
\description{
  Allow the user to set and examine a variety of \emph{options}
  which affect operations of the Infusion package. 
  However, typically these should not be modified, and if they are, not more than once in a data analysis.
}
\usage{
Infusion.options(\dots)

Infusion.getOption(x)

}
\arguments{
  \item{x}{a character string holding an option name.}
  \item{\dots}{A named value, or several of them, or a single unnamed argument which is a named list). 
  The following values, with their defaults, are used in \code{Infusion}:

  \describe{
\item{\code{version}}{A \code{\link{package_version}} object or something that can be converted to such an object by \code{\link{as.package_version}}. The default is the installed version of \pkg{Infusion}. Alternative values may be used by some functions (such as \code{\link{get_workflow_design}}) to use control values of the algorithms from other versions the the package. 
;   }
\item{\code{mixturing}}{character string: package or function to be used for mixture modelling. Recognized packages are \code{"Rmixmod"} (the default) and \code{"mclust"};   }
\item{\code{train_cP_size}:}{Expression for \code{train_cP_size} argument of \code{\link{project.character}}.}
\item{\code{trainingsize}:}{Expression for \code{trainingsize} argument of \code{\link{project.character}}.}
\item{\code{projKnotNbr = 1000}:}{default value of \code{trainingsize} argument of \code{\link{project.character}} for REML (as implied by default expression for \code{trainingsize}).}
\item{\code{logLname = "logL"}:}{default value of \code{logLname} argument of \code{\link{infer_logLs}}. The name given to the inferred log likelihoods in all analyses.}
\item{\code{LRthreshold= - qchisq(0.999,df=1)/2}:}{A value used internally by \code{\link{sample_volume}} to sample points 
      in the upper region of the likelihood surface, as defined by the given likelihood ratio threshold.
}
\item{\code{precision = 0.25}:}{default value of \code{precision} argument of \code{\link{refine}}. Targets RMSE of log L and log LR estimates.}
\item{\code{nRealizations=1000}:}{default value of \code{nRealizations} argument of \code{\link{add_simulation}}. Number of realizations for each empirical distribution.}
\item{\code{mixmodGaussianModel="Gaussian_pk_Lk_Ck"}:}{default models used in clustering by \code{Rmixmod}. Run \code{Rmixmod::mixmodGaussianModel()} for a list of possible models, and see the statistical documentation (Mixmod Team 2016) for explanations about them.}
\item{\code{global_strategy_args}:}{list of arguments for \code{Rmixmod::mixmodStrategy()}.}
\item{\code{nbClu_pow_rule_fn = function(nr) {....}} and\cr
      \code{maxnbCluster = function(projdata, ...) {....} }:}{
functions controlling the value of \code{nbCluster} used in clustering by \code{Rmixmod} or other mixture modelling backends. See \code{\link{get_nbCluster_range}} for details.      
}
\item{\code{example_maxtime=2.5}:}{Used in the documentation to control whether the longer examples should be run. 
The approximate running time of given examples (or some very rough approximation for it) on one author's laptop is compared to this value.} 
\item{\code{nb_cores}}{Number of cores for parallel computations (see Details for implementation of these).}
\item{\code{gof_nstats_fn}}{See \code{\link{goftest}}.}
} 
and possibly other undocumented values for development purposes.
  }

}
\details{
  \code{Infusion} can perform parallel computations if several cores are available and requested though \code{Infusion.options(nb_cores=.)}. If the \code{doSNOW} back-end is attached (by explicit request from the user), it will be used; otherwise, \code{pbapply} will be used. Both provide progress bars, but \code{doSNOW} may provide more efficient load-balancing. The character shown in the progress bar is \code{'P'} for parallel via \code{doSNOW} backend, \code{'p'} for parallel via \code{pbapply} functions, and \code{'s'} for serial via \code{pbapply} functions. In addition, \code{add_simulation} can parallelise at two levels: at an outer level over parameter point, or at an inner level over simulation replicates for each parameter point. The progress bar of the outer computation is shown, but the character shown in the progress bar is \code{'N'} if the inner computation is parallel via the \code{doSNOW} backend, and \code{'n'} if it is parallel via \code{pbapply} functions. So, one should see either \code{'P'} or \code{'N'} when using \code{doSNOW}.
}
\references{
Mixmod Team (2016). Mixmod Statistical Documentation. Université de Franche-Comté,
Besançon, France. Version: February 10, 2016 retrieved from \url{https://www.mixmod.org}.
}
\value{
  For \code{Infusion.getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{Infusion.options()}, a list of all set options.  For
  \code{Infusion.options(name)}, a list of length one containing the set value,
  or \code{NULL} if it is unset.  For uses setting one or more options,
  a list with the previous values of the options changed (returned
  invisibly).
}
\examples{
  Infusion.options()
  Infusion.getOption("LRthreshold")
  \dontrun{
  Infusion.options(LRthreshold=- qchisq(0.99,df=1)/2)
  }
}
