\name{plot1Dprof}
\alias{plot1Dprof}
\alias{plot2Dprof}
\title{
Plot likelihood profiles
}
\description{
These functions plot 1D and 2D profiles from an SLik object
}
\usage{
plot1Dprof(object, pars=object$colTypes$fittedPars, type="logLR",   
           gridSteps=21, xlabs=list(), ylab, scales=NULL,
           plotpar=list(pch=20))
plot2Dprof(object, pars=object$colTypes$fittedPars, type="logLR",  
           gridSteps=17, xylabs=list(), main, scales=NULL,
           plotpar=list(pch=20), margefrac = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An SLik object }
  \item{pars}{ The parameters for which profiles will be computed. For 2D plots, all pairs of parameters in \code{pars} are considered }
  \item{type}{ \code{logLR} or \code{LR} or \code{logL} for the (log) likelihood ratio or the log likelihood}
  \item{gridSteps}{ The number of values (in each dimension for 2D plots) which likelihood  should be computed. For 1D plots, 
    \code{gridSteps=0} will call \code{curve} rather than a plot of points 
  }
  \item{xlabs}{ A \emph{list} of alternative axis labels. The names  of the list elements should be elements of \code{pars} (see Examples)}
  \item{xylabs}{ Same as \code{xlabs} but affecting both axes in 2D plots}
  \item{ylab}{ Same as \code{ylab} argument of \code{plot}. Default depends on \code{type} argument.}
  \item{main}{ Same as \code{main} argument of \code{plot}. Default depends on \code{type} argument.}
  \item{scales}{ A named character vector, which controls ticks and tick labels on axes, so that these can be expressed as (say) the exponential of the parameter inferred in the SLik object. 
    For example if the likelihood of \code{logPop} = log(population size) was thus inferred, \code{scales=c(logPop="log")} will give population size values on the axis (but will retain a log scale for this parameter).  
  }
  \item{plotpar}{ arguments for \code{par()} such as font sizes, etc.}
  \item{margefrac}{ For development purposes, not documented.}
}
\value{
No return value.
}
\examples{
if (Infusion.getOption("example_maxtime")>40) {
 data(densv)
 slik <- infer_surface(densv) ## infer a log-likelihood surface
 slik <- MSL(slik) ## find the maximum of the log-likelihood surface
 plot1Dprof(slik,pars="s2",gridSteps=40,xlabs=list(s2=expression(paste(sigma^2))))
}
}
