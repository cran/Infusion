\name{rparam}
\alias{rparam}
\alias{sample_volume}
\title{
  Sample the parameter space
}
\description{
    These functions take an \code{SLik} object (as produced by \code{\link{MSL}}) and samples its parameter space in (hopefully) clever ways, not yet well documented. \code{rparam} calls \code{sample_volume} to define points targeting the likelihood maximum and the bounds of confidence intervals, with \code{n} for these different targets dependent on the mean square error of prediction of likelihood at the maximum and at CI bounds.
}
\usage{
rparam(object, n= 1, useEI = list(max=TRUE,profileCI=TRUE,rawCI=FALSE), 
       useCI = TRUE, verbose = interactive(), tryn=30*n,  
       level = 0.95, CIweight=Infusion.getOption("CIweight"))

sample_volume(object, n = 6, useEI, vertices=NULL,
              dlr = NULL, verbose = interactive(), 
              fixed = NULL, tryn= 30*n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} object
}
  \item{n}{
The number of parameter points to be produced 
}
  \item{useEI}{ List of booleans, each determining whether to use an \dQuote{expected improvement} (EI) criterion (e.g. Bingham et al., 2014) to select candidate parameter points to better ascertain a particular focal point. The elements \code{max}, \code{profileCI} and \code{rawCI} determine this for three types of focal points, respectively the MSL estimate, profile CI bounds, and full-dimensional bounds. When EI is used, \code{n} points with best EI are selected among \code{tryn} points randomlygenerated in some neighborhood of the focal point. 
}
\item{vertices}{Points are sampled within a convex hull defined by \code{vertices}. By default, these vertices are taken from \code{object$fit$data}.  
}
  \item{useCI}{
  Whether to define points targeting the bounds of confidence intervals for the parameters. An expected improvement criterion is also used here.
}
  \item{level}{
    If \code{useCI} is \code{TRUE} but confidence intervals are not available from the \code{object}, such intervals are computed with coverage \code{level}. 
  }
  \item{dlr}{
  A (log)likelihood ratio threshold used to select points in the upper region of the likelihood surface. Default value is 
  given by \code{Infusion.getOption("LRthreshold")}   
  }
  \item{verbose}{
  Whether to display some information about selection of points, or not
  }
  \item{fixed}{
    A list or named vector, of which each element is of the form \code{<parameter name>=<value>}, defining a one-dimensional constraint in parameter space.
    Points will be sampled in the intersection of the volume defined by the \code{object} and of such constraint(s).
  }
  \item{tryn}{
  See \code{useEI} argument. 
  }
  \item{CIweight}{
  For development purposes, not documented. 
  }
}
\value{
a data frame of parameter points. Only parameters variable in the \code{SLik} object are considered. 
}
\references{
D. Bingham, P. Ranjan, and W.J. Welch (2014) Design of Computer Experiments for Optimization, Estimation of Function Contours, and Related Objectives, pp. 109-124 in Statistics in Action: A Canadian Outlook (J.F. Lawless, ed.). Chapman and Hall/CRC.
}
\examples{
if (Infusion.getOption("example_maxtime")>10) {
 data(densv)
 summliksurf <- infer_surface(densv) ## infer a log-likelihood surface
 sample_volume(summliksurf)
}
}
