\name{refine}
\alias{refine}
\alias{refine.default}
\alias{refine.SLik}
\alias{refine.SLik_j}
\alias{refine.SLikp}
\title{
  Refine estimates iteratively.    
}
\description{
% FR->FR must beupdated for SLik_j
  This is a generic function with currently methods for \code{SLik}, \code{SLik_j} and \code{SLikp} objects (as produced by \code{\link{MSL}}). Depending on the value of its \code{newsimuls} argument, and on whether the function used to generate empirical distributions can be called by R, it (1) defines new parameters points and/or (2) infers their summary likelihood or tail probabilities for each parameter point independently, adds the inferred values results as input for refined inference of likelihood or P-value response surface, and provides new point estimates and confidence intervals.
}
\usage{
\method{refine}{SLik}(object, method=NULL, ...)
%\method{refine}{SLik_j}(object,...)
%\method{refine}{SLikp}(object,...)
\method{refine}{default}(object, surfaceData, Simulate =
            attr(surfaceData,"Simulate"), maxit = 1, n = NULL, 
            useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), 
            newsimuls = NULL, useCI = TRUE, level = 0.95, 
            verbose = list(most=interactive(),final=NULL,movie=FALSE,proj=FALSE),
            precision = Infusion.getOption("precision"),
            nb_cores = NULL, packages=attr(object$logLs,"packages"), 
            env=attr(object$logLs,"env"), method, 
            eval_RMSEs=TRUE, update_projectors = FALSE,
            cluster_args=list(),
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} object
}
  \item{surfaceData}{
A data.frame with attributes, usually taken from the \code{object} and thus \bold{not} specified by user, usable as input for \code{\link{infer_surface}}.  
}
  \item{Simulate}{
    Character string: name of the function used to simulate sample. The only meaningful non-default value is \code{NULL}, in which case \code{refine} may return (if \code{newsimuls} is also \code{NULL}) a data frame of parameter points on which to run a simulation function.
  }
  \item{maxit}{
  Maximum number of iterative refinements (see also \code{precision} argument)
}
  \item{n}{
A number of parameter points (excluding replicates and confidence interval points), whose likelihood should be computed 
(see \code{n} argument of \code{\link{sample_volume}})   
}
  \item{useEI}{
  Cf this argument in \code{\link{rparam}}
}
  \item{newsimuls}{
  For the \code{SLik_j} method, a matrix or data frame, with the same parameters and summary statistics as the \code{data} of the original \code{\link{infer_SLik_joint}} call.
  
  For other methods, a \code{list} of simulation of distributions of summary statistics, in the same format as for \code{link{add_simulation}}. 
  If no such list is provided (i.e., if \code{newsimuls} remains \code{NULL}), the \code{attr(object$logLs,"Simulate")} function is used (it is inherited from the \code{Simulate} 
  argument of \code{\link{add_simulation}} through the initial sequence of calls of functions \code{add_simulation},
  \code{infer_logLs} or \code{infer_tailp}, and \code{infer_surface}). If no such function is available, then this function returns parameters for which new distribution should be provided by the user.
}
  \item{useCI}{
  whether to include parameter points near the inferred confidence interval points in the set of points which likelihood should be computed 
}
  \item{level}{
  Intended coverage of confidence intervals 
}
  \item{verbose}{ A list as shown by the default, or simply a vector of booleans. \code{verbose$most} controls whether to display information about progress and results, except plots; \code{$final} controls whether to \code{plot()} the final \code{object} to show the final likelihood surface. Default is to plot it only in an interactive session and if fewer than three parameters are estimated; \code{$movie} controls whether to \code{plot()} the updated \code{object} in each iteration; \code{verbose$proj} controls the \code{verbose} argument of \code{\link{project.character}}. If \code{verbose} is a vector of booleans, they are matched to as many elements from \code{"most","movie","final","proj"}, in that order. 
  }
  \item{precision}{
    Requested local precision of surface estimation, in terms of prediction standard errors (RMSEs) of both the maximum summary log-likelihood and the likelihood ratio at any CI bound available. Iterations will stop when either \code{maxit} is reached, or if the RMSEs have been computed for the object (see \code{eval_RMSEs} argument) and this precision is reached for the RMSEs.
    A given precision on the CI bounds themselves might seem more interesting, but is not well specified by a single precision parameter if the parameters are on widely different scales.
  }
  \item{cluster_args}{A list of arguments, passed to \code{\link[parallel]{makeCluster}}. If thus affects parallel computations in the functions to which it may be passed: \code{add_simulation}, \code{MSL} (for bootstrap computations) or \code{infer_logLs}. It may contain a non-null \code{spec} element, in which case the distinct \code{nb_cores} argument of \code{refine} is ignored. If \code{spec} is a \bold{named} integer, it has an additional effect on \code{\link{add_simulation}}.}
  \item{nb_cores}{Shortcut for \code{cluster_args$spec}.}
  \item{packages}{NULL or a list with possible elements \code{add_simulation} and \code{logL_method}, passed respectively as the \code{packages} arguments of \code{add_simulation} and \code{infer_logLs}, wherein they are the additional packages to be loaded on child processes. The default value keeps pre-\code{refine} values over iterations.}
  \item{env}{An environment, passed as the \code{env} argument to \code{add_simulation}. The default value keeps the pre-\code{refine} value over iterations.}
  \item{method}{(A vector of) suggested method(s) for estimation of smoothing parameters (see \code{method} argument of \code{\link{infer_surface}}). The ith element of the vector is
  used in the ith iteration, if available; otherwise the last element is used. This argument is not always heeded, in that REML may be used if the suggested method is GCV but it appears to perform poorly. The default for \code{SLikp}, \code{SLik_j}, and \code{SLikp} objects are \code{"REML"}, \code{"mixmodCluster"}, and \code{"PQL"}, respectively. 
  }
  \item{eval_RMSEs}{passed to \code{\link{MSL}}}
  \item{update_projectors}{Boolean; whether to update the projectors at each iteration.}
  \item{\dots}{
further arguments passed to or from other methods. \code{refine} passes these arguments to the \code{plot} method suitable for the \code{object}. 
}
}
\details{
New parameter points are sampled as follows: the algorithm aims to sample uniformly the space of parameters contained in the confidence regions defined by the \code{level} argument, and to surround it by a region sampled proportionally to likelihood. In each iteration the algorithm aims to add as many points (say \emph{n}) as computed in the first iteration, so that after \emph{k} iterations of \code{refine}, there are \eqn{n * (k+1)} points in the simulation table. However, when not enough points satisfy certain criteria, only \emph{n/5} points may be added in an iteration, this being compensated in further iterations. For example, if \eqn{n=600}, the table may include only 720 points after the first refine, but 1800 after the second.
}
\value{
  An updated \code{SLik} or \code{SLik_j} object.
}
\note{See workflow examples in \code{\link{example_reftable}} and \code{\link{example_raw}}}
\examples{
  ## see Note.
}
