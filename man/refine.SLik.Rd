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
            newsimuls = NULL, trypoints=NULL, CIs = useCI, useCI = TRUE, level = 0.95, 
            verbose = list(most=interactive(),final=NULL,movie=FALSE,proj=FALSE),
            precision = Infusion.getOption("precision"),
            nb_cores = NULL, packages=attr(object$logLs,"packages"), 
            env=attr(object$logLs,"env"), method,  using = object$using, 
            eval_RMSEs=TRUE, update_projectors = FALSE,
            cluster_args=list(),
            cl_seed=.update_seed(object),
            nbCluster=quote(refine_nbCluster(nr=nrow(data))),
            ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} or \code{SLik_j} object
}
  \item{surfaceData}{
A data.frame with attributes, usually taken from the \code{object} and thus \bold{not} specified by user, usable as input for \code{\link{infer_surface}}.  
}
  \item{Simulate}{
    Character string: name of the function used to simulate samples. The only meaningful non-default value is \code{NULL}, in which case \code{refine} may return (if \code{newsimuls} is also \code{NULL}) a data frame of parameter points on which to run a simulation function.
  }
  \item{maxit}{
  Maximum number of iterative refinements (see also \code{precision} argument)
}
  \item{n}{
NULL or numeric, for a number of parameter points (excluding replicates and confidence interval points in the primitive workflow), whose likelihood should be computed 
(see \code{n} argument of \code{\link{sample_volume}}). This argument is typically not heeded in the first refinement iteration (only one fifth as many points may be produced), but will be closely approached in later ones (so four refinement iterations with \code{n=1000} is expected to produce 3200 new points). If \code{n} is left NULL, the number of points of the initial reference table is used as a reference, but with a somewhat different effect: four refinement iterations starting from a reference table of \code{1000} ones iis expected to produce 4000 new points (though again, possibly only 200 in the first refinement iteration).      
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
\item{trypoints}{
 A data frame of parameters on which the simulation function \code{attr(object$logLs,"Simulate")} should be called to extend the reference table. Only for programming by expert users, because poorly thought input \code{trypoints} could severely affect the inferences.     
}
  \item{CIs}{
  Boolean: whether to infer bounds of (one-dimensional, profile) confidence intervals. Their computation is not quite reliable in parameter spaces of large dimensions, so they should not be trusted per se, yet they may be useful for the definition of new parameter points.   
}
  \item{useCI}{
  whether to include parameter points near the inferred confidence interval points in the set of points whose likelihood should be computed. Effective only if \code{CIs} was TRUE.
}
  \item{level}{
  Intended coverage of confidence intervals 
}
  \item{verbose}{ A list as shown by the default, or simply a vector of booleans. \code{verbose$most} controls whether to display information about progress and results, except plots; \code{$final} controls whether to \code{plot()} the final \code{object} to show the final likelihood surface. Default is to plot it only in an interactive session and if fewer than three parameters are estimated; \code{$movie} controls whether to \code{plot()} the updated \code{object} in each iteration; \code{verbose$proj} controls the \code{verbose} argument of \code{\link{project.character}}. If \code{verbose} is an unnamed vector of booleans, they are matched to as many elements from \code{"most","movie","final","proj"}, in that order. 
  }
  \item{precision}{
    Requested local precision of surface estimation, in terms of prediction standard errors (RMSEs) of both the maximum summary log-likelihood and the likelihood ratio at any CI bound available. Iterations will stop when either \code{maxit} is reached, or if the RMSEs have been computed for the object (see \code{eval_RMSEs} argument) and this precision is reached for the RMSEs.
    A given precision on the CI bounds themselves might seem more interesting, but is not well specified by a single precision parameter if the parameters are on widely different scales.
  }
  \item{nb_cores}{Shortcut for \code{cluster_args$spec} for sample simulation.}
  \item{cluster_args}{A list of arguments for \code{\link[parallel]{makeCluster}}, in addition to \code{makeCluster}'s \code{spec} argument which is in most cases best specified by the \code{nb_cores} argument. Cluster arguments allow independent control of parallel computations for the different steps of a \code{refine} iteration (see Details; as a rough but effective summary, use only \code{nb_cores} when the simulations support it, and only \code{cluster_args=list(project=list(num.threads=<.>))} when they do not).}
  \item{packages}{NULL or a list with possible elements \code{add_simulation} and \code{logL_method}, passed respectively as the \code{packages} arguments of \code{add_simulation} and \code{infer_logLs}, wherein they are the additional packages to be loaded on child processes. The default value keeps pre-\code{refine} values over iterations.}
  \item{env}{An environment, passed as the \code{env} argument to \code{add_simulation}. The default value keeps the pre-\code{refine} value over iterations.}
  \item{using}{Passed to \code{\link{infer_SLik_joint}}: a charcter string used to control the joint-density estimation method, as documented for that function. Default is to use to same method as in the the first iteration, but this argument allows a change of method.}
  \item{method}{(A vector of) suggested method(s) for estimation of smoothing parameters (see \code{method} argument of \code{\link{infer_surface}}), and therefore controlling the primitive workflow (see \code{using} instead for controlling the up-to-date workflow). The ith element of the vector is
  used in the ith iteration, if available; otherwise the last element is used. This argument is not always heeded, in that REML may be used if the suggested method is GCV but it appears to perform poorly. The default for \code{SLikp} and \code{SLikp} objects are \code{"REML"} and \code{"PQL"}, respectively. 
  }
  \item{eval_RMSEs}{passed to \code{\link{MSL}}}
  \item{update_projectors}{Boolean; whether to update the projectors at each iteration.}
    \item{cl_seed}{NULL or integer, passed to \code{add_simulation}. The default code uses an internal function, \code{.update_seed}, to update it from a previous iteration.}
    \item{nbCluster}{Passed to \code{\link{infer_SLik_joint}}. The \code{data} in the expression for the default value refers to the \code{data} argument of the latter function.}
  \item{\dots}{
further arguments passed to or from other methods. \code{refine} passes these arguments to the \code{plot} method suitable for the \code{object}. 
}
}
\details{
New parameter points are sampled as follows: the algorithm aims to sample uniformly the space of parameters contained in the confidence regions defined by the \code{level} argument, and to surround it by a region sampled proportionally to likelihood. In each iteration the algorithm aims to add as many points (say \emph{n}) as computed in the first iteration, so that after \emph{k} iterations of \code{refine}, there are \eqn{n * (k+1)} points in the simulation table. However, when not enough points satisfy certain criteria, only \emph{n/5} points may be added in an iteration, this being compensated in further iterations. For example, if \eqn{n=600}, the table may include only 720 points after the first refine, but 1800 after the second.

Independent control of parallelisation may be needed in the different steps, e.g. if the simulations are not easily parallelised whereas the projection method natively handles parallelisation. In the up-to-date workflow with default \code{ranger} projection method, prarallelisation controls may be passed to \code{add_reftable} for sample simulations, to \code{project} methods when projections are updated, and to \code{MSL} for RMSE computations (alternatively for the primitive workflow, \code{add_simulation}, \code{infer_logLs} and \code{MSL} are called). \code{nb_cores}, if given and not overcome by other options, will control simulation and projection steps (but not RMSE computation): \code{nb_cores} gives the number of parallel processes for sample simulation, with additional \code{makeCluster} arguments taken from \code{cluster_args}, but RMSE computations are performed serially. Further independent control is possible as follows:\cr
\code{cluster_args=list(project=list(num.threads=<.>))} allows control of the \code{num.threads} argument of \pkg{ranger} functions;\cr 
\code{cluster_args=list(RMSE=list(spec=<number of 'children'>))} can be used to force parallel computation of RMSEs;\cr
\code{cluster_args=list(spec=<.>, <other makeCluster arguments>))} would instead apply the same arguments to both reference table and RMSE computation, overcoming the default effect of \code{nb_cores} in both of them; finally\cr 
\code{cluster_args=list(reftable=list(<makeCluster arguments>),RMSEs=list(<makeCluster arguments>))} allows full independent control of parallelisation for the two computations.
}
\value{
  \code{refine} returns an updated \code{SLik} or \code{SLik_j} object.
}
\note{See workflow examples in (by order of decreasing relevance) \code{\link{example_reftable}}, \code{\link{example_raw_proj}} and \code{\link{example_raw}}.}
\examples{
  ## see Note for links to examples.
}
