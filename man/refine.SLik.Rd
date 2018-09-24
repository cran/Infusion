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
  This is a generic function with currently methods for \code{SLik} and \code{SLikp} objects (as produced by \code{\link{MSL}}). Depending on the value of its \code{newsimuls} argument, and on whether the function used to generate empirical distributions can be called by R, it (1) defines new parameters points and/or (2) infers their summary likelihood or tail probabilities for each parameter point independently, adds the inferred values results as input for refined inference of likelihood or P-value response surface, and provides new point estimates and confidence intervals.
}
\usage{
\method{refine}{SLik}(object, method=NULL, ...)
%\method{refine}{SLik_j}(object,...)
%\method{refine}{SLikp}(object,...)
\method{refine}{default}(object, surfaceData, maxit = 1, n = NULL, 
            useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), 
            newsimuls = NULL, useCI = TRUE, level = 0.95, 
            verbose = list(most=interactive(),movie=FALSE),
            precision = Infusion.getOption("precision"),
            nb_cores = NULL, packages=attr(object$logLs,"packages"), 
            env=attr(object$logLs,"env"), method, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} object
}
  \item{surfaceData}{
A data.frame with attributes, usually taken from the \code{object} and thus \bold{not} specified by user, usable as input for \code{\link{infer_surface}}.  
}
  \item{maxit}{
  Maximum number of iterative refinements (see also \code{precision} argument)
}
  \item{n}{
A number of parameter points (excluding replicates and confidence interval points), which likelihood should be computed 
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
  \item{verbose}{ A list as shown by the default, or simply a vector of booleans, indicating respectively
  whether to display (1) some information about progress and results; (2) a plot of a the likelihood surface, if a one- or two-parameters surface is inferred. \code{verbose$movie} controls whether such plots are produced for all iterations. \code{verbose$most} controls other verbosity, including whether a final plot is produced. 
  }
  \item{precision}{
    Requested local precision of surface estimation, in terms of prediction variance of both the maximum summary log-likelihood and the likelihood ratio at any CI bound available. Iterations will stop when either \code{maxit} is reached, or this precision is reached for these MSEs.
    A given precision on the CI bounds themselves might seem more interesting, but is not well specified by a single precision parameter if the parameters are on widely different scales.
  }
  \item{nb_cores}{integer (possibly named, see \code{\link{add_simulation}}), passed to \code{add_simulation} and \code{infer_logLs}.}
  \item{packages}{NULL or a list with possible elements \code{add_simulation} and \code{logL_method}, passed respectively as the \code{packages} arguments of \code{add_simulation} and \code{infer_logLs}, wherein they are the additional packages to be loaded on child processes. The default value keeps pre-\code{refine} values over iterations.}
  \item{env}{An environment, passed as the \code{env} argument to \code{add_simulation}. The default value keeps the pre-\code{refine} value over iterations.}
  \item{method}{(A vector of) suggested method(s) for estimation of smoothing parameters (see \code{method} argument of \code{\link{infer_surface}}). The ith element of the vector is
  used in the ith iteration, if available; otherwise the last element is used. This argument is not always heeded, in that REML may be used if the suggested method is GCV but it appears to perform poorly. The default for \code{SLikp}, \code{SLik_j}, and \code{SLikp} objects are \code{"REML"}, \code{"mixmodCluster"}, and \code{"PQL"}, respectively. 
  }
  \item{\dots}{
further arguments passed to or from other methods. \code{refine} passes these arguments to the \code{plot} method suitable for the \code{object}. 
}
}
\value{
  An updated \code{SLik} or \code{SLik_j} object.
}
\examples{
## see main documentation page for the package
}
