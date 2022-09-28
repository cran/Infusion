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
  \item{\dots}{A named value or a list of named values. The following values, with their defaults, 
    are used in \code{Infusion}:

  \describe{
\item{\code{mixturing}}{character string: package or function to be used for mixture modelling. Recognized packages are \code{"Rmixmod"} (the default) and \code{"mclust"};   }
\item{\code{train_cP_size}:}{Expression for \code{train_cP_size} argument of \code{\link{project.character}}.}
\item{\code{trainingsize}:}{Expression for \code{trainingsize} argument of \code{\link{project.character}}.}
\item{\code{projKnotNbr = 1000}:}{default value of \code{trainingsize} argument of \code{\link{project.character}} for REML (as implied by default expression for \code{trainingsize}).}
\item{\code{logLname = "logL"}:}{default value of \code{logLname} argument of \code{\link{infer_logLs}}. The name given to the inferred log likelihoods in all analyses.}
\item{\code{LRthreshold= - qchisq(0.999,df=1)/2}:}{A value used internally by \code{\link{sample_volume}} to sample points 
      in the upper region of the likelihood surface, as defined by the given likelihood ratio threshold.
}
\item{\code{precision = 0.1}:}{default value of \code{precision} argument of \code{\link{refine}}. Targets RMSE of log L and log LR estimates.}
\item{\code{nRealizations=1000}:}{default value of \code{nRealizations} argument of \code{\link{add_simulation}}. Number of realizations for each empirical distribution.}
\item{\code{mixmodGaussianModel="Gaussian_pk_Lk_Dk_A_Dk"}:}{default models used in clustering by \code{Rmixmod}. Run \code{Rmixmod::mixmodGaussianModel()} for a list of possible models, and see the statistical documentation (Mixmod Team 2016) for explanations about them.}
\item{\code{seq_nbCluster= function(projdata, nr=nrow(projdata)) {seq(ceiling(nr^0.3))}}:}{function to control the value of \code{nbCluster} used in clustering by \code{Rmixmod} (see Details for discussion of this default).}
\item{\code{maxnbCluster = function(projdata) {...} }:}{function to control the maximum number of clusters (see Details).}
\item{\code{example_maxtime=2.5}:}{Used in the documentation to control whether the longer examples should be run. 
The approximate running time of given examples (or some very rough approximation for it) on one author's laptop is compared to this value.} 
\item{\code{nb_cores}}{Number of cores for parallel computations (see Details for implementation of these).}
\item{\code{gof_nstats_fn}}{See \code{\link{goftest}}.}
} 
and possibly other undocumented values for development purposes.
  }

}
\details{
  %Up to version 1.3.0, the default upper value of the \code{nbCluster} range was the value recommended (\eqn{n^0.3} as function fo the number of points \eqn{n}) in the \code{mixmod} statistical documentation (Mixmod Team, 2016). However, simulations have suggested that this is conservative (i.e., it underfitted the data, resulting in systematic biases in inferences of likelihood surfaces), hence \eqn{2 n^0.3} is used now. 
  The set of the number of clusters tried (\code{nbCluster} argument in \code{Rmixmod}) is controlled by two options: \code{seq_nbCluster} and \code{maxnbCluster}. The second is used to correct the first, using the dimensions of the \code{projdata} locally used for clustering, which typically differs from the dimensions of the user-level \code{data} (if projections have been applied, in particular). The default upper value of the \code{nbCluster} range is the value recommended in the \code{mixmod} statistical documentation (Mixmod Team, 2016). But this default may be suitable only for low-dimensional data, hence the need for correcting it by\code{maxnbCluster}.
  
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
