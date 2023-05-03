\name{plot1Dprof}
\alias{plot1Dprof}
\alias{plot2Dprof}
\title{
Plot likelihood profiles
}
\description{
These functions plot 1D and 2D profiles from a summary-likelihood object.

High quality 2D plots may be slow to compute, and there may be many of them in high-dimensional parameter spaces, so parallelization of the computation of each profile point has been implemented for them. Usual caveats apply: there is an time cost of launching processes on a cluster, particularly on socket clusters, possibly offsetting the benefits of parallelization when each profile point is fast to evaluate. Further, summary-likelihood objects are typically big (memory-wise), which may constrain the number of concurrent processes.    
}
\usage{
plot1Dprof(object, pars=object$colTypes$fittedPars, type="logLR",   
           gridSteps=21, xlabs=list(), ylab, scales=NULL,
           plotpar=list(pch=20), 
           control=list(min=-7.568353, shadow_col="grey70"),
           decorations = function(par) NULL, ...)
plot2Dprof(object, pars=object$colTypes$fittedPars, type="logLR",  
           gridSteps=17, xylabs=list(), main, scales=NULL,
           plotpar=list(pch=20), margefrac = 0,
           decorations = function(par1,par2) NULL, 
           filled.contour.fn = "spaMM.filled.contour", 
           cluster_args=NULL, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{SLik} or \code{SLik_j} object}
  \item{pars}{ Control of parameters for which profiles will be computed. If \code{pars} is specified as a vector of names, profiles are plotted for each parameter, or (2D case) for all pairs of distinct parameters. Finer control is possible in the 2D case (see Details).}
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
  \item{plotpar}{Arguments for \code{par()} such as font sizes, etc.}
  \item{control}{Control of \code{"zoom"} or \code{"dual"} plots (see Details).}
  \item{decorations}{ A function with formals parameters as shown by the default (being parameters names represented as character strings), implementing graphic directives added to the plot.}
  \item{margefrac}{ For development purposes, not documented.}
  \item{filled.contour.fn}{Name of a possible alternative to \code{graphics::filled.contour} to be used for rendering the plot.}
  \item{cluster_args}{NULL, or a list in which case a cluster may be created and used. The list elements must match the arguments \code{spec} and \code{type} of \code{parallel::makeCluster}. A socket cluster is created unless \code{type="FORK"} (on operating systems that support fork clusters).}
  \item{\dots}{ Further arguments passed by another function. Currently these arguments are ignored.
  }
}
\details{
In the 2D case, the \code{pars} may be specified as a two-column natrix, in which case profiles are generated for all pairs of distinct parameters specified by rows of the matrix. It may also be specified as a two-element list, where each element is a vector of parameter names. In that case, profiles are generated for all pairs of distinct parameters combining one element of each vector.

A \code{"zoom"} plot shows only the top part of the profile, defined as points whose y-values are above a threshold minus-log-likelihood ratio \code{control$min}, whose default is -7.568353, the 0.9999 p-value threshold.

A \code{"dual"} plot displays both the zoom, and a shadow graph showing the full profile. The dual plot is shown only when requested and if there are values above and below \code{control$min}. The shadow curve color is given by \code{control$shadow_col}.
}
\value{
Both functions return a list, which currently has a single element \code{MSL_updated} which is a boolean indicating whether the summary-likelihood maximum (but not the intervals) has been recomputed.
}
\examples{

if (Infusion.getOption("example_maxtime")>20) { 
  #### Toy bivariate gaussian model, three parameters, no projections
  #
  myrnorm2 <- function(mu1,mu2,s2,sample.size) {
    sam1 <- rnorm(n=sample.size,mean=mu1,sd=sqrt(s2))
    sam2 <- rnorm(n=sample.size,mean=mu2,sd=sqrt(s2))
    s <- c(sam1,sam2)
    e_mu <- mean(s)
    e_s2 <- var(s)
    c(mean=e_mu,var=e_s2,kurt=sum((s-e_mu)^4)/e_s2^2)
  } 
  #
  ## pseudo-sample, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- myrnorm2(mu1=4,mu2=2,s2=1,sample.size=40) ## 
  #
  ## build reference table
  parsp <- init_reftable(lower=c(mu1=2.8,mu2=1,s2=0.2), 
                         upper=c(mu1=5.2,mu2=3,s2=3),
                         nUnique=600)
  parsp <- cbind(parsp,sample.size=40)
  simuls <- add_reftable(Simulate="myrnorm2",par.grid=parsp)
  
  ## Inferring the summary-likelihood surface...
  densv <- infer_SLik_joint(simuls,stat.obs=Sobs)
  slik_j <- MSL(densv) ## find the maximum of the log-likelihood surface
  
  ## plots
  plot2Dprof(slik_j,gridSteps=21,
             ## alternative syntaxes for non-default 'pars':
             # pars = c("mu1","mu2"), # => all combinations of given elements
             # pars = list("s2",c("mu1","mu2")), # => combinations via expand.grid()
             # pars = matrix(c("mu1","mu2","s2","mu1"), ncol=2), # => each row of matrix
             xylabs=list(
               mu1=expression(paste(mu[1])),
               mu2=expression(paste(mu[2])),
               s2=expression(paste(sigma^2))
             )) 
  # One could also add (e.g.) 
  #            cluster_args=list(spec=4, type="FORK"), 
  # when longer computations are requested.
}

if (Infusion.getOption("example_maxtime")>40) {
 #### Older example with primitive workflow 
 data(densv)
 slik <- infer_surface(densv) ## infer a log-likelihood surface
 slik <- MSL(slik) ## find the maximum of the log-likelihood surface
 plot1Dprof(slik,pars="s2",gridSteps=40,xlabs=list(s2=expression(paste(sigma^2)))) 
}

}
