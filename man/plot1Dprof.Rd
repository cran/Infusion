\name{plot1Dprof}
\alias{plot1Dprof}
\alias{plot2Dprof}
\title{
Plot likelihood profiles
}
\description{
These functions plot 1D and 2D profiles from a summary-likelihood object.

If you feel the 1D profiles are ugly, see the Details. Confidence intervals may still be correct in such cases.

High quality 2D plots may be slow to compute, and there may be many of them in high-dimensional parameter spaces, so parallelization of the computation of each profile point has been implemented for them. Usual caveats apply: there is an time cost of launching processes on a cluster, particularly on socket clusters, possibly offsetting the benefits of parallelization when each profile point is fast to evaluate. Further, summary-likelihood objects are typically big (memory-wise), notably when theyinclude many projections, and this may constrain the number of processes that can run in parallel.

Parallelization is also implemented for 1D profiles, but over the parameter for which profiles are computed, rather than over points in a profile. So it is effective only if profiles are computed for several parameters.  

In default 2D plots, some areas may be left blank, for distinct reasons: the function values may be too low (as controlled by the \code{min_logLR} argument), or the likelihood profile may be maximized at parameter values which do not satisfy constraints defined by the \code{constr_crits} of the \code{\link{infer_SLik_joint}} function. Low profile values, even when shown, are not accurate anyway, since the inference workflow aims at inferring the top of likelihood \dQuote{hill}  relevant for the computation of confidence regions.
}
\usage{
plot1Dprof(object, pars=object$colTypes$fittedPars, fixed=NULL, 
           type="logLR", gridSteps=21, xlabs=list(), ylab, scales=NULL,
           plotpar=list(pch=20), 
           control=list(min=-7.568353, shadow_col="grey70"),
           decorations = function(par) NULL, 
           profiles=NULL, add=NULL, safe=TRUE,
           cluster_args=NULL, do_plot=TRUE, CIlevels=NULL,
           lower=NULL, upper=NULL,
           verbose=TRUE,
           ...)
plot2Dprof(object, pars=object$colTypes$fittedPars, fixed=NULL, 
           type="logLR", gridSteps=17, xylabs=list(), main, scales=NULL,
           plotpar=list(pch=20), margefrac = 0,
           decorations = function(par1,par2) NULL, 
           filled.contour.fn = "spaMM.filled.contour", cluster_args=NULL, 
           min_logLR = qchisq(0.95,df=length(object$colTypes$fittedPars))/2 +3,
           lower=NULL, upper=NULL, color.palette=NULL, profiles=NULL,
           ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An \code{SLik} or \code{SLik_j} object}
  \item{pars}{ Control of parameters for which profiles will be computed. If \code{pars} is specified as a vector of names, profiles are plotted for each parameter, or (2D case) for all pairs of distinct parameters. 
  
  Finer control is possible in the 2D case: the \code{pars} may be specified as a two-column natrix, in which case profiles are generated for all pairs of distinct parameters specified by rows of the matrix. It may also be specified as a two-element list, where each element is a vector of parameter names. In that case, profiles are generated for all pairs of distinct parameters combining one element of each vector.}
  \item{fixed}{A named vector of parameter values to keep fixed in the profile computation: these parameters will be excluded from the \code{pars} as well as from the paramters over which maximization occurs.}
  \item{type}{
  Character: possible values are \code{"logL"} to plot the log-likelihood profile; \code{"logLR"} (or \code{"LR"} for the not-log version) to plot the log-likelihood-ratio profile; and, for 1D profiles, \code{"zoom"}, \code{"ranges"} and \code{"dual"} which are variants of the \code{"logLR"} plot controlling the plot ranges in different ways (see Details). }
  \item{gridSteps}{ The number of values (in each dimension for 2D plots) for which likelihood  should be computed. For 1D plots, 
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
  \item{decorations}{ A function implementing graphic directives added to the plot (anything that is not a function is converted to the default function). Its formal parameters are as shown by its trivial default value, the \code{par}[.] arguments being parameters names (as a vector of character strings), to allow the additional plot elements to depend on the parameter of each subplot (see Plot 3 in Examples).}
  \item{margefrac}{ For development purposes, not documented.}
  \item{safe}{ For development purposes, not documented.} % FALSE -> nlminb with gradient of 'unsafe' logL is used
  \item{filled.contour.fn}{Name of a possible alternative to \code{graphics::filled.contour} to be used for rendering the plot.}
  \item{cluster_args}{NULL, or a list in which case a cluster may be created and used. The list elements must match the arguments \code{spec} and \code{type} of \code{parallel::makeCluster}. A socket cluster is created unless \code{type="FORK"} (on operating systems that support fork clusters). For \code{plot1Dprof}, parallelisation is over the parametrs for which profiles are computed; For \code{plot2Dprof},it is over the grid of points for each profile.}
  \item{profiles}{The \code{profiles} element of the return value of a previous call of the function. The point coordinates it provides for a given parameter or pair of parameters will be used, instead of recomputing the profile.}
  \item{add}{The \code{profiles} element of the return value of a previous \code{plot1Dprof()} call. The point coordinates it provides for a given parameter will be added to the profile produced by other arguments, using the graphic directives specified by the \dots\ to distingush them.}
  \item{min_logLR}{Numeric; the minimal value of the log-likelihood ratio to be shown on 2D plots of type \code{logLR}, parameter regions having lower ratios being left blank.
  }
  \item{do_plot}{Boolean: whether to actually plot the profiles or only return computations for it.}
  \item{CIlevels}{NULL, or a vector of confidence levels whose bounds will be added to the plot. Suggested values are \code{c(0.99, 0.95,0.90,0.75)}}
  \item{lower, upper}{Numeric vectors of bounds for the parameters.}
  \item{verbose}{Boolean: whether to print information about the progress of the computation.}
  \item{color.palette}{
  Either NULL or a function that can replace the default color function used by \code{plot2Dprof}. The function must have a single argument, giving the number of color levels. 
}
  \item{\dots}{ Further arguments passed by another function. Currently these arguments are ignored, except when handling the \code{add} argument, where they are passed to \code{graphics::lines}.
  }
}
\details{
\bold{Possible issues in computing the profiles:}
Computation of profiles is complicated by local maximization issues, which may result in highly visible artefacts in the plots. \code{plot1Dprof} tries to reduce their impact by computing the profile points starting from parameter values closest to the identified likelihood-maximizing value, and by using a method for selecting initial values for maximization which is partially controlled by the result of computation of the previous profile point. A second sequential computation of points profile may additionally be performed, this time starting from parameter values closest to the identified 95\%\ confidence intervals rather than from the MLEs, when such intervals are available in the fit object or if they are requested by the \code{CIlevels} argument.  

\bold{Graphic details of the plots:} 
The 2D profile plots will typically include a contour level for the 95\% 2D confidence region, labelled \dQuote{2D CI} on the scale.

A \code{"zoom"} plot shows only the top part of the profile, defined as points whose y-values are above a threshold minus-log-likelihood ratio \code{control$min}, whose default is -7.568353, the 0.9999 p-value threshold.

A \code{"ranges"} (for \code{plot1Dprof}) is similar to \code{"zoom"} but maintains a fixed y range, facilitating comparison of profiles with and among different \code{plot1Dprof} calls. Few of the originally computed points may appear in the retained x range, in which case it may be useful to generate additional points ensured to appear in the plot, using the \code{CIlevels} argument.

A \code{"dual"} plot displays both the zoom, and a shadow graph showing the full profile. The dual plot is shown only when requested and if there are values above and below \code{control$min}. The shadow curve color is given by \code{control$shadow_col}.

}
\value{
Both functions return a list, invisibly for \code{plot1Dprof}. The list has elements\cr
\code{ * MSL_updated} which is a boolean indicating whether the summary-likelihood maximum has been recomputed (if it is TRUE, a message is printed);\cr
\code{ * profiles}, itself a list which stores information about each profile. The format of the information per profile is not yet stable (subject to changes without notice), % ____F I X M E____
but consistent with the one handled by the \code{profiles} and \code{add} argument). Its elements currently include\cr
\code{    - }the x,y or x,y,z coordinates of the putative plot;\cr 
\code{    - }the full coordinates of the profile points in a \code{profpts} matrix (1D plot) or a 3D array (2D plot; first dimension are the fitted parameters, second and third are x and y steps). 

\code{plot1Dprof} may have the side effect of adding confidence interval information to the fit object.
}
\examples{

if (Infusion.getOption("example_maxtime")>20) { # 2D plots relatively slow
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
  ## simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- myrnorm2(mu1=4,mu2=2,s2=1,sample.size=40) ## 
  #
  ## build reference table
  parsp <- init_reftable(lower=c(mu1=2.8,mu2=1,s2=0.2), 
                         upper=c(mu1=5.2,mu2=3,s2=3))
  parsp <- cbind(parsp,sample.size=40)
  simuls <- add_reftable(Simulate="myrnorm2", parsTable=parsp)
  
  ## Inferring the summary-likelihood surface...
  densvj <- infer_SLik_joint(simuls,stat.obs=Sobs)
  slik_j <- MSL(densvj) ## find the maximum of the log-likelihood surface
  
  ### plots
  # Plot 1: a 1D profile:
  prof1 <- plot1Dprof(slik_j,pars="s2",gridSteps=40,xlabs=list(s2=expression(paste(sigma^2))))
  
  # Using 'add' for comparison of successive profiles:
  slik_2 <- refine(slik_j, n=600)
  # Plot 2: comparing 1D profiles of different iterations:
  prof2 <- plot1Dprof(slik_2,"s2", gridSteps=40,xlabs=list(s2=expression(paste(sigma^2))), 
                      add=prof1$profiles, col="grey30") 
                      
  # Plot 3: using 'decorations' and 'profiles' for recycling previous computation for s2:
  DGpars <- c(mu1=4,mu2=2,s2=1,sample.size=40) # data-generating parameters
  plot1Dprof(slik_2, gridSteps=40,xlabs=list(s2=expression(paste(sigma^2))), 
             profiles=prof2$profiles,
   decorations=function(par) {
      points(y=-7,x=DGpars[par],pch=20,cex=2,col="red");
      points(y=-0.5,x=slik_2$MSL$MSLE[par],pch=20,cex=2,col="blue")
   }
  ) 

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

if (Infusion.getOption("example_maxtime")>5) {
 #### Older example with primitive workflow 
 data(densv)
 slik <- infer_surface(densv) ## infer a log-likelihood surface
 slik <- MSL(slik) ## find the maximum of the log-likelihood surface
 prof1 <- plot1Dprof(slik,pars="s2",gridSteps=40,xlabs=list(s2=expression(paste(sigma^2)))) 
}

}
