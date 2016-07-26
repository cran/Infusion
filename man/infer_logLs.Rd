\name{infer_logLs}
\alias{infer_logLs}
\alias{infer_tailp}
\alias{infer_logL_by_GLMM}
\alias{infer_logL_by_Rmixmod}
\alias{infer_logL_by_mclust}
\alias{infer_logL_by_Hlscv.diag}
\title{
Infer log Likelihoods using simulated distributions of summary statistics
}
\description{
  For each simulated distribution of summary statistics, \code{infer_logLs} infers a probability density function, and the density of the observed values of the summary statistics is deduced. By default, inference of each density is performed by \code{infer_logL_by_Rmixmod}, which fits a distribution of summary statistics using procedures from the \code{\link{Rmixmod}} package.        
}
\usage{
infer_logLs(object, stat.obs, 
            logLname = Infusion.getOption("logLname"), 
            verbose = list(most=interactive(), 
                           final=FALSE), 
            method="infer_logL_by_Rmixmod",
            ...)
infer_tailp(object, refDensity, stat.obs,
                tailNames=Infusion.getOption("tailNames"),
                verbose=interactive(), method=NULL,...)
infer_logL_by_GLMM(EDF,stat.obs,logLname,verbose)
infer_logL_by_Rmixmod(EDF,stat.obs,logLname,verbose)
infer_logL_by_mclust(EDF,stat.obs,logLname,verbose)
infer_logL_by_Hlscv.diag(EDF,stat.obs,logLname,verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  A list of simulated distributions (the return object of \code{\link{add_simulation}})
}
\item{EDF}{
  An empirical distribution, with a required \code{par} attribute (an element of the \code{object} list).
}
\item{stat.obs}{
  Named numeric vector of observed values of summary statistics.
}
 \item{logLname}{
 The name to be given to the log Likelihood in the return object, or the root of the latter name in case of conflict with other names in this object.
 }
  \item{tailNames}{
  Names of \dQuote{positives} and \dQuote{negatives} in the binomial response for the inference of tail probabilities. 
}
  \item{refDensity}{
  An object representing a reference density (such as an \code{\link{HLfit}} fit object or other objects with a similar \code{predict} method) which, together with the density inferred from each empirical density, defines a likelihood ratio used to define a rejection region.   
}
  \item{verbose}{ A list as shown by the default, or simply a vector of booleans, indicating respectively
  whether to display (1) some information about progress; (2) a final summary of the results after all elements of \code{simuls} have been processed. If a count of 'outlier'(s) is reported, this typically means that \code{stat.obs} is not within the envelope of a simulated distribution (or whatever other meaning the user attaches to an \code{FALSE isValid} code: see Details)
  }
  \item{method}{
  A function for density estimation. See Description for the default behaviour and Details for the constraints on input and output of the function. 
}
  \item{\dots}{
further arguments passed to or from other methods (currently not used).
}
}
\value{
For \code{infer_logLs}, a data frame containing parameter values and their log likelihoods, and additional information such as attributes providing information about the parameter names and statistics names (not detailed here). These attributes are essential for further inferences.

See Details for the required value of the \code{method}s called by \code{infer_logLs}.
}
\details{
By default, density estimation is based on \code{Rmixmod} methods. Other available methods are not routinely used and not all of \code{Infusion} features may work with them. The function \code{\link{mixmodCluster}} is called, with arguments \code{nbCluster=Infusion.getOption("nbCluster")} and \code{mixmodGaussianModel=Infusion.getOption("mixmodGaussianModel")}. If \code{Infusion.getOption("nbCluster")} specifies a sequence of values, then several clusterings are computed and AIC is used to select among them.  

\code{infer_logL_by_GLMM}, \code{infer_logL_by_Rmixmod}, \code{infer_logL_by_mclust}, and \code{infer_logL_by_Hlscv.diag} are examples of the method that may be provided for density estimation. Other \code{method}s may be provided with the same arguments. Their return value must include the element \code{logL}, an estimate of the log-density of \code{stat.obs}, and the element \code{isValid} with values \code{FALSE}/\code{TRUE} (or 0/1). The standard format for the return value is \code{unlist(c(attr(EDF,"par"),logL,isValid=isValid))}.

\code{isValid} is primarily intended to indicate whether the log likelihood of \code{stat.obs} inferred by a given density estimation method was suitable input for inference of the likelihood surface. \code{isValid} has two effects: to distinguish points for which isValid is FALSE in the plot produced by \code{\link{plot.SLik}}; and more critically, to control the sampling of new parameter points within \code{\link{refine}} so that points for which isValid is FALSE are less likely to be sampled. 

Invalid values may for example indicate a likelihood estimated as zero (since log(0) is not suitable input), or (for density estimation methods which may infer erroneously large values when extrapolating), whether \code{stat.obs} is within the convex hull of the EDF. In user-defined \code{method}s, invalid inferred logL should be replaced by some alternative low estimate, as all methods included in the package do.    

The source code of \code{infer_logL_by_Hlscv.diag} illustrates how to test whether \code{stat.obs} is within the convex hull of the EDF, using functions \code{resetCHull} and \code{isPointInCHull} (exported from the \code{blackbox} package).

\code{infer_logL_by_Rmixmod} calls \code{\link[Rmixmod]{mixmodCluster}}, \code{infer_logL_by_mclust} calls \code{\link[mclust]{densityMclust}}, 
\code{infer_logL_by_Hlscv.diag} calls \code{\link[ks]{kde}}, and \code{infer_logL_by_GLMM} fits a binned distribution of summary statistics using a Poisson GLMM with autocorrelated random effects, where the binning is based on a tesselation of a volume containing the whole simulated distribution. Limited experimentations so far suggest that the mixture models methods are fast and appropriate (\code{Rmixmod}, being a bit faster, is the default method); that the kernel smoothing method is more erratic and moreover requires additional input from the user, hence is not really applicable, for distributions in dimension \emph{d}= 4 or above; and that the GLMM method is a very good density estimator for \emph{d}=2 but will challenge one's patience for \emph{d}=3 and further challenge the computer's memory for \emph{d}=4.
}
\seealso{
 See step (3) of the workflow in the Example on the main \code{\link{Infusion}} documentation page.  
}
