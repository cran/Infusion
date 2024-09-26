\name{refine}
\alias{refine}
\alias{refine.default}
\alias{refine.SLik}
\alias{refine.SLik_j}
\alias{refine.SLikp}
\alias{reproject}
\alias{recluster}
\title{
  Refine estimates iteratively    
}
\description{
  \code{refine} is a generic function with methods for objects of the classes produced by \code{\link{MSL}}. In the up-to-date workflow, it can automatically (1) define new parameters points, (2) add simulations to the reference table for these points, (3) optionally recompute projections, (4) update the inference of the likelihood surface, and (5) provides new point estimates, confidence intervals, and other results of an \code{\link{MSL}} call. It can repeat these steps iteratively as controlled by its \code{workflow_design}. 
Although it has many control arguments, few of them may be needed in any application. In particular it is designed to use reasonable default controls for the number of iterations, the number of points added in each iteration, and whether to update projections or not, when given only the current fit object as input.
  
\code{reproject} and \code{recluster} are wrappers for \code{refine(..., ntot=0L)}, updating the object after either recomputing the projections or only re-performing the multivariate gaussian mixture clustering.
}
\usage{
\method{refine}{SLik}(object, method=NULL, ...)
%\method{refine}{SLik_j}(object,...)
%\method{refine}{SLikp}(object,...)
\method{refine}{default}(
    object, 
    ##       reference table simulations  
    Simulate = attr(surfaceData,"Simulate"),
    control.Simulate = attr(surfaceData,"control.Simulate"),
    newsimuls = NULL,
    ##       CIs
    CIs = workflow_design$reftable_sizes[useCI], 
    useCI = prod(dim(object$logLs))<12000L, level = 0.95,
    ##       workflow design
    workflow_design = get_workflow_design(
        npar=length(fittedPars), n_proj_stats = length(statNames),
        n_latent=length(latentVars)), 
    maxit, ntot= maxit*.get_size_first_iter(object), n=NULL,
    ##       termination conditions 
    precision = Infusion.getOption("precision"),
    eval_RMSEs = workflow_design$reftable_sizes,
    ##       verbosity
    verbose = list(notable=TRUE, most=interactive(),final=NULL, movie=FALSE,
                 proj=FALSE, rparam=NULL, progress_bars=interactive()),
    ##       projection controls
    update_projectors = NULL,
    methodArgs = list(),
    ##       Likelihood surface modeling (up-to-date workflow)
    using = object$using, 
    nbCluster = quote(refine_nbCluster(nr=nrow(data))),
    ##       parallelisation
    cluster_args = list(), nb_cores=NULL, env=get_from(object,"env"), 
    packages = get_from(object,"packages"), cl_seed=.update_seed(object),
    ##       obscure stuff
    target_LR = NULL,  
    ##       not explicitly needed in up-to-date workflow
    trypoints = NULL,
    surfaceData, 
    method, 
    useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE),
    rparamFn = Infusion.getOption("rparamFn"),
    ## 
    ... 
)
            
reproject(object, eval_RMSEs = NULL, CIs = NULL, ...)
recluster(object, eval_RMSEs = NULL, CIs = NULL, update_projectors=FALSE, ...)
}
\arguments{
  \item{object}{
an \code{SLik} or \code{SLik_j} object
}
\bold{## reference table simulations}
  \item{Simulate}{
    Character string: name of the function used to simulate samples. As it is typically stored in the object this argument does not need to be explicitly given; otherwise this should be the same function provided to \code{\link{add_reftable}}, whose documentation details the design requirements. The only meaningful non-default value is \code{NULL}, in which case \code{refine} may return (if \code{newsimuls} is also \code{NULL}) a data frame of parameter points on which to run a simulation function.
  }
  \item{control.Simulate}{
    A list of arguments of the \code{Simulate} function (see \code{\link{add_simulation}}). The default value should be used unless you understand enough of its structure to modify it wisely (e.g., it may contain the path of an executable on one machine and a different path may be specified to refine a fit on another machine).  
  }
  \item{newsimuls}{
  For the \code{SLik_j} method, a matrix or data frame, with the same parameters and summary statistics as the \code{data} of the original \code{\link{infer_SLik_joint}} call.
  
  For other methods, a \code{list} of simulation of distributions of summary statistics, in the same format as for \code{link{add_simulation}}. 
  If no such list is provided (i.e., if \code{newsimuls} remains \code{NULL}), the function extracted by \code{get_from(object,"Simulate")} is used (it is inherited from the \code{Simulate} 
  argument of \code{\link{add_simulation}} through the initial sequence of calls of functions \code{add_simulation},
  \code{infer_logLs} or \code{infer_tailp}, and \code{infer_surface}). If no such function is available, then this function returns parameters for which new distribution should be provided by the user.
}
\bold{## CIs}
  \item{CIs}{
  Boolean, or boolean vector, or numeric (preferably integer) vector: controls to infer bounds of (one-dimensional, profile) confidence intervals. The numeric vector form allows to specify reference table size(s) for which CIs should be computed when these sizes are first reached. TRUE or FALSE will force or inhibit computation in all iterations. Finally (and probably less useful), a boolean vector such as \code{CIs=c(TRUE,FALSE,TRUE)} requests computation of CIs when the number of points cumulatively added reaches the target number of points for the first, third, and any subsequent iterations up to \code{maxit} (this may differ in certain cases from the first, third, and so on, iterations: see Details).  

  The default for \code{refine} is described in the Details. The default for \code{reproject} is to update the CIs if there are computed ones within the input \code{object}.
}
  \item{useCI}{
  whether to perform RMSE computations for inferred confidence interval points.
}
  \item{level}{
  Intended coverage of confidence intervals 
}
\bold{## workflow design}
  \item{workflow_design}{
    A list structured as the return value of \code{\link{get_workflow_design}}. The default value makes reference to elements of the input object's \code{colTypes} element.}
  \item{maxit}{
  Maximum number of iterative refinements (see also \code{precision} argument).
}
  \item{ntot}{
    NULL or numeric: control of the total number of simulated samples (one for each new parameter point) to be added to the reference table over the \code{maxit} iterations.  See Details for the rules used to determine the number of points added in each iteration. Reasonable default values are defined for \code{ntot} and \code{maxit} (see Details), so that beginners (and ideally, even more advanced users) do not have to find good values.

\code{ntot=0L} may be used to re-generate the projectors or the clustering without augmenting the reference table.
  }
  \item{n}{
NULL or numeric, for a number of parameter points (excluding replicates and confidence interval points in the primitive workflow), whose likelihood should be computed in each iteration 
(see \code{n} argument of \code{\link{sample_volume}}). Slightly less intutive alternative to \code{ntot} specification, as there is at least one iteration where the actual number of added points is not the nominal \code{n} (see Details).   
\code{n=0L} will have the same effect as \code{ntot=0L}.
}
\bold{## termination conditions}
  \item{precision}{
    Requested local precision of surface estimation, in terms of prediction standard errors (RMSEs) of both the maximum summary log-likelihood and the likelihood ratio at any CI bound available. Iterations will stop when either \code{maxit} is reached, or if the RMSEs have been computed for the object (see \code{eval_RMSEs} argument) and this precision is reached for the RMSEs.
    A given precision on the CI bounds themselves might seem more interesting, but is not well specified by a single precision parameter if the parameters are on widely different scales.
  }
  \item{eval_RMSEs}{Same usage as for \code{CIs}; controls the \code{eval_RMSEs} argument of \code{\link{MSL}} in each iteration. See Details for the default. 
  The default for \code{reproject} is to update the RMSEs if there are computed ones within the input \code{object}.}
\bold{## verbosity}
  \item{verbose}{ A list as shown by the default, or simply a vector of booleans. \code{verbose$most} controls whether to display information about progress and results, except plots; \code{$final} controls whether to \code{plot()} the final \code{object} to show the final likelihood surface. Default is to plot it only in an interactive session and if fewer than three parameters are estimated; \code{$movie} controls whether to \code{plot()} the updated \code{object} in each iteration; \code{verbose$proj} controls the \code{verbose} argument of \code{\link{project.character}}; \code{verbose$rparam} controls (cryptic) information about generation of new parameter points; \code{verbose$progress_bars} controls display of some progress bars. If \code{verbose} is an unnamed vector of booleans, they are interpreted as as-many first elements of the\code{verbose} vector, in the order shown by the default. 
  }
\bold{## projection controls}
  \item{update_projectors}{Same usage as for \code{CIs}; this controls in which iterations the projectors are updated. The default \code{NULL} value is strongly recommended. See Details for further explanations.}
  \item{methodArgs}{A list of arguments for the projection method. By default the \code{methodArgs} of the original \code{\link{project.character}} calls are reused over iteration, but elements of the new \code{methodArgs} list will be used to update the original \code{methodArgs}. Note that the updated list becomes the new default for further iterations.}
\bold{## Likelihood surface modeling}
  \item{using}{Passed to \code{\link{infer_SLik_joint}}: a character string used to control the joint-density estimation method, as documented for that function (see \code{method} instead for equivalent control in primitive workflow). Default is to use to same method as in the the first iteration, but this argument allows a change of method.}
    \item{nbCluster}{Passed to \code{\link{infer_SLik_joint}}. The \code{data} in the expression for the default value refers to the \code{data} argument of the latter function.}
\bold{## parallelisation}
  \item{cluster_args}{A list of arguments for \code{\link[parallel]{makeCluster}}, in addition to \code{makeCluster}'s \code{spec} argument which is in most cases best specified by the \code{nb_cores} argument. Cluster arguments allow independent control of parallel computations for the different steps of a \code{refine} iteration (see Details; as a rough but effective summary, use only \code{nb_cores} when the simulations support it, and see the \code{methodArgs} argument if independent control of parallelisation of the projection procedure is needed).}
    \item{nb_cores}{Integer: shortcut for specifying \code{cluster_args$spec} for sample simulation.}
  \item{packages}{NULL or a list with possible elements \code{add_simulation} and \code{logL_method} (the latter for the primitive workflow). These elements should be formatted as the \code{packages} arguments of \code{add_simulation} and \code{infer_logLs}, respectively, wherein they are the additional packages to be loaded on child processes. The effect of the default value of this argument is to pass over successive \code{refine} calls the value stored in the input fit object (itself determined by the latest use of the \code{packages} argument in, e.g., \code{add_simulation} or in previous \code{refine}s).}
  \item{env}{An environment, passed as the \code{env} argument to \code{add_simulation}. The default value keeps the pre-\code{refine} value over iterations.}
    \item{cl_seed}{NULL or integer, passed to \code{add_simulation}. The default code uses an internal function, \code{.update_seed}, to update it from a previous iteration.}
\bold{## others}
    \item{target_LR}{Likelihood ratio threshold used to control the sampling of new points and the selection of points for projections. Do not change it unless you known what you are doing.}
  \item{method}{For the primitive workflow: (a vector of) suggested method(s) for estimation of smoothing parameters (see \code{method} argument of \code{\link{infer_surface}}). The ith element of the vector is
  used in the ith iteration, if available; otherwise the last element is used. This argument is not always heeded, in that REML may be used if the suggested method is GCV but it appears to perform poorly. The default for \code{SLikp} and \code{SLikp} objects are \code{"REML"} and \code{"PQL"}, respectively. 
  }
\item{trypoints}{
 A data frame of parameters on which the simulation function \code{get_from(object,"Simulate")} should be called to extend the reference table. Only for programming by expert users, because poorly thought input \code{trypoints} could severely affect the inferences.     
}
  \item{useEI}{
  for the primitive workflow only: cf this argument in \code{\link{rparam}}.
}
      \item{surfaceData}{
 for the primitive workflow only: a data.frame with attributes, usually taken from the \code{object} and thus \bold{not} specified by user, usable as input for \code{\link{infer_surface}}.  
}
\item{rparamFn}{
  Function used to sample new parameter values.
}
  \item{\dots}{
further arguments passed to or from other methods. \code{refine} passes these arguments to the \code{plot} method suitable for the \code{object}. 
}
}
\details{
\code{ * }\bold{Controls of exploration of parameter space}: New parameter points are sampled so as to fill the space of parameters contained in the confidence regions defined by the \code{level} argument, and to surround it by a region sampled proportionally to likelihood. 

Each \code{refine} call performs several iterations, these iterations stopping when \code{ntot} points have been added to the simulation table. The target number of points potentially added in each iteration is controlled by the \code{ntot} and \code{maxit} arguments as described below, but fewer points may be actually added in each iteration, and more than \code{maxit} iterations may be needed to add the \code{ntot} points, if in a given iteration too few \dQuote{good} candidate points are generated according to the internal rules for sampling the parameter region with high likelihood.  In that case, the next iteration tries to keep up with the missing points by adding more points than the target number, but if not enough points have been added after \code{maxit} iterations, further iterations will be run. 

CIs and RMSEs may be computed in any iteration but the default values of \code{eval_RMSEs} and \code{CIs} are chosen so as to avoid performing these computations too often, particularly when they are expected to be slow. The default implies that the RMSE for the maximum logL will be computed at the end each block of iterations that defines a refine (itself defined to reach to reference table sizes specified by the \code{workflow_design} and its default value). If the reference table is not too large (see default value of \code{useCI} for the precise condition), RMSEs of the logL are also computed at the inferred bounds of profile-based confidence intervals for each parameter. 

Although the \code{update_projectors} argument allow similar control of the iterations where projections are updated, it is advised to keep it NULL (default value), so that whether projectors are updated in a given iteration is controlled by default internal rules. Setting it to \code{TRUE} would induce updating whenever any of the target reference table sizes implied by the \code{workflow_design$subblock_sizes} is reached. The default \code{NULL}, as the same effect subject to additional conditions: updating may not be performed when the training set is considered too similar to the one used to compute pre-existing projections, or when the train set includes more samples than the limit define by the global package option \code{upd_proj_subrows_thr}

Default values of \code{ntot} and \code{maxit} are controlled by the value of the \code{workflow_design}, which itself has the shown default value, and are distinct for the first vs. subsequent \code{refine}s. 
The target number of points in each iteration is also controlled differently for the first vs. subsequent \code{refine}s. This design is motivated by the fact that the likelihood surface is typically poorly inferred in the first refine so that the parameter points sampled then tend to be less relevant than those that can be sampled in later iterations. In the first \code{refine} call, the target number of points increases roughly as powers of two over iterations, to reach \code{ntot} cumulatively after \code{maxit} iterations. The default \code{ntot} is twice the size of the initial reference table, and the default \code{maxit} is 5. The \code{\link{example_reftable}} Example illustrates this, where the initial reference table holds 200 simulations, and the default target number of points to be added in 5 iterations by the first \code{refine} call are 25, 25, 50, 100 and 200.  In later \code{refine} calls, the target number is \code{ntot/maxit} in each iteration.

\code{ * }Independent \bold{control of parallelisation} may be needed in the different steps, e.g. if the simulations are not easily parallelised whereas the projection method natively handles parallelisation. In the up-to-date workflow with default \code{ranger} projection method, distinct parallelisation controls may be passed to \code{add_reftable} for sample simulations, to \code{project} methods when projections are updated, and to \code{MSL} for RMSE computations (alternatively for the primitive workflow, \code{add_simulation}, \code{infer_logLs} and \code{MSL} are called). 
The most explicit way of specifying distinct controls is by a list structured as \preformatted{
cluster_args=list(reftable=list(<makeCluster arguments>),
                  RMSEs=list(<makeCluster arguments>))
} 
A \code{project=list(num.threads=<.>)} element can be added to this list, providing control of the \code{num.threads} argument of \pkg{ranger} functions. However, this is retained mainly for back compatibility as the \code{methodArgs} argument can now be used to specify the \code{num.threads}.

Simpler arguments may be used and will be interpreted as follows: \code{nb_cores}, if given and not overriden by a \code{spec} argument in \code{cluster_args} (or in sublists of it), will control simulation and projection steps (but not RMSE computation): that is, \code{nb_cores} then gives the number of parallel processes for sample simulation, with additional \code{makeCluster} arguments taken from \code{cluster_args}, but RMSE computations are performed serially. On the other hand, a \code{spec} argument in
\code{cluster_args=list(spec=<.>, <other makeCluster arguments>))} will instead apply the same arguments to both reference table and RMSE computation, overcoming the default effect of \code{nb_cores} in both of them.
}
\value{
  \code{refine} returns an updated \code{SLik} or \code{SLik_j} object, unless both \code{newsimuls} and \code{Simulate} arguments are NULL, in which case a data frame of parameter points is returned.
}
\note{See workflow examples in (by order of decreasing relevance) \code{\link{example_reftable}}, \code{\link{example_raw_proj}} and \code{\link{example_raw}}.

See \code{\link{get_workflow_design}}, the function that controls 
the default value of the \code{workflow_design} argument, and can be used to provide non-default controls.
}
\examples{
  ## see Note for links to examples.
}
