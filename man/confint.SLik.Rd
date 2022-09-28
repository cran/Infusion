\name{confint.SLik}
\alias{confint.SLik}
\alias{confint.SLik_j}
\alias{confint.SLikp}
\alias{confint}
\title{
  Compute confidence intervals by (profile) summary likelihood 
}
\description{
  This takes an \code{SLik} object (as produced by \code{\link{MSL}}) and deduces confidence bounds for each parameter, using a (profile, if relevant) likelihood ratio method.
}
\usage{
\method{confint}{SLik}(object, parm,
                       level=0.95, verbose=interactive(),
                       fixed=NULL,which=c(TRUE,TRUE),...)
%\method{confint}{SLik_j}(object, parm,
%                       level=0.95, verbose=interactive(),
%                       fixed=NULL,which=c(TRUE,TRUE),...)
%\method{confint}{SLikp}(object, parm,
%                       level=0.95, verbose=interactive(),
%                       fixed=NULL,which=c(TRUE,TRUE),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} or \code{SLik_j} object
}
  \item{parm}{
  The parameter which confidence bounds are to be computed
}
  \item{level}{
  The desired coverage of the interval
 }
  \item{verbose}{
  Whether to print some information or not
 }
  \item{fixed}{When this is \code{NULL} the computed interval is a profile confidence interval over all parameters excluding \code{parm}.
     \code{fixed} allows one to set fixed values to some of these parameters.  
}
  \item{which}{
  A pair of booleans, controlling whether to compute respectively the lower and the upper CI bounds.
 }
  \item{\dots}{
further arguments passed to or from other methods (currently not used).
}
}
\value{
A list with sublists for each parameter, each sublist containing of three vectors: the bounds of the one-dimensional confidence interval; the \dQuote{full} 
(only parameters variable in the \code{SLik} object are considered) parameter point for the lower bound, and the full parameter point for the upper bound
}
\examples{
## see main documentation page for the package
}
