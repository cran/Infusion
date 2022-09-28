\name{MSL}
\alias{MSL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Maximum likelihood from an inferred likelihood surface
}
\description{
This computes the maximum of an object of class \code{SLik} representing an inferred (summary) likelihood surface}
\usage{
MSL(object, CIs = TRUE, level = 0.95, verbose = interactive(),
    eval_RMSEs = TRUE, cluster_args=list(),init=NULL, prior_logL=NULL,
    ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{SLik_j} as produced by \code{\link{infer_SLik_joint}} (or, in the primitive workflow,  of class \code{SLik} as produced by \code{\link{infer_surface.logLs}}).
}
  \item{CIs}{
  If \code{TRUE}, construct one-dimensional confidence intervals for all parameters.
}
  \item{level}{
  Intended coverage probability of the confidence intervals.
}
  \item{verbose}{
  Whether to display some information about progress and results.
}
  \item{eval_RMSEs}{
  Logical: whether to evaluate prediction uncertainty for likelihoods/ likelihood ratios/ parameters. 
  %By default TRUE for \code{SLik} objects, and best kept so, as it is necessary for the automated iterative method. Default may be FALSE for other classes of objects, in particular \code{SLik_j}.  
}
\item{cluster_args}{A list of arguments, passed to \code{\link[parallel]{makeCluster}}, to control parallel computation of RMSEs. Beware that parallel computation of RMSEs tends to be memory-intensive. The list may contain a non-null \code{spec} element, in which case the \code{nb_cores} global \pkg{Infusion} option is ignored. Do \bold{*not*} use a structured list with an \code{RMSE} element as is possible for \code{refine} (see Details of \code{\link{refine}} documentation).}

  \item{init}{Initial value for the optimiser. Better ignored.
  }
  \item{prior_logL}{(effective only for up-to-date workflow using gaussian mixture modelling of a joint distribution of parameters and statistics) a function that returns a vector of prior log-likelihood values, which is then added to the likelihood deduced from the summary likelihood analysis. The function's single argument must handle a matrix similar to the \code{newdata} argument of \code{\link{predict.SLik_j}}. 
  }
  \item{\dots}{Further arguments passed from or to other methods. 
  }
}
\value{
  The \code{object} is returned invisibly, with the following added members, each of which being (as from version 1.5.0) an environment:
  \describe{
    \item{\code{MSL}}{containing variables \code{MSLE} and \code{maxlogL} that match the \code{par} and \code{value} returned by an \code{optim} call. Also contain the \code{hessian} of summary likelihood at its maximum.}
    \item{\code{RMSEs}}{containing, as variable \code{RMSEs}, the root mean square errors of the log-likelihood at its inferred maximum and of the log-likelihood ratios at the CI bounds.}
    \item{\code{par_RMSEs}}{containing, as variable \code{par_RMSEs}, root mean square errors of the CI bounds.}
  }
  To ensure backward-compatibility of code to possible future changes in the structure of the objects, the extractor function \code{\link{get_from}} should be used to extract the \code{RMSEs} and \code{par_RMSEs} variables from their respective environments, and more generally to extract any element from the objects.
}
\details{
If Kriging has been used to construct the likelihood surface, \code{RMSEs} are computed using approximate formulas for prediction (co-)variances in linear mixed midels (see Details in \code{\link[spaMM]{predict}}). Otherwise, a more computer-intensive bootstrap method is used.   
\code{par_RMSEs} are computed from \code{RMSEs} and from the numerical gradient of profile log-likelihood at each CI bound. Only \code{RMSEs}, not \code{par_RMSEs}, are compared to \code{precision}.
}
\examples{
## see main documentation page for the package
}
