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
           plotpar=list(pch=20), 
           control=list(min=-7.568353, shadow_col="grey70"))
plot2Dprof(object, pars=object$colTypes$fittedPars, type="logLR",  
           gridSteps=17, xylabs=list(), main, scales=NULL,
           plotpar=list(pch=20), margefrac = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An SLik object }
  \item{pars}{ The parameters for which profiles will be computed. For 2D plots, all pairs of parameters in \code{pars} are considered }
  \item{type}{Character: \code{"logL"} to plot the log-likelihood profile; \code{"logLR"} (or \code{"LR"} for the not-log version) to plot the log-likelihood-ratio profile (the default); or \code{"zoom"} or \code{"dual"} for variants of \code{"logLR"} (see details). }
  \item{gridSteps}{ The number of values (in each dimension for 2D plots) which likelihood  should be computed. For 1D plots, 
    \code{gridSteps=0} is now obsolete.
  }
  \item{xlabs}{ A \emph{list} of alternative axis labels. The names  of the list elements should be elements of \code{pars} (see Examples)}
  \item{xylabs}{ Same as \code{xlabs} but affecting both axes in 2D plots}
  \item{ylab}{ Same as \code{ylab} argument of \code{plot}. Default depends on \code{type} argument.}
  \item{main}{ Same as \code{main} argument of \code{plot}. Default depends on \code{type} argument.}
  \item{scales}{ A named character vector, which controls ticks and tick labels on axes, so that these can be expressed as (say) the exponential of the parameter inferred in the SLik object. 
    For example if the likelihood of \code{logPop} = log(population size) was thus inferred, \code{scales=c(logPop="log")} will give population size values on the axis (but will retain a log scale for this parameter). Possible values of each element of the vector are \code{"identity"} (default), \code{"log"}, and \code{"log10"},   
  }
  \item{plotpar}{ arguments for \code{par()} such as font sizes, etc.}
  \item{control}{Control of \code{"zoom"} or \code{"dual"} plots (see Details).}
  \item{margefrac}{ For development purposes, not documented.}
}
\details{
A \code{"zoom"} plot shows only the top part of the profile, defined as points whose y-values are above a threshold minus-log-likelihood ratio \code{control$min}, whose default is -7.568353, the 0.9999 p-value threshold.

A \code{"dual"} plot displays both the zoom, and a shadow graph showing the full profile. The dual plot is shown only when requested and if there are values above and below \code{control$min}. The shadow curve color is given by \code{control$shadow_col}.
}
\value{
Both functions return a list, which currently has a single element \code{MSL_updated} which is a boolean indicating whether the summary-likelihood maximum (but not the intervals) has been recomputed.
}
\examples{
if (Infusion.getOption("example_maxtime")>40) {
 data(densv)
 slik <- infer_surface(densv) ## infer a log-likelihood surface
 slik <- MSL(slik) ## find the maximum of the log-likelihood surface
 plot1Dprof(slik,pars="s2",gridSteps=40,xlabs=list(s2=expression(paste(sigma^2)))) 
}
}
