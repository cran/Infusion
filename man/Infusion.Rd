\name{Infusion}
% ----- \Rdversion{1.1}
\alias{Infusion}
\alias{Infusion-package}

\docType{package}

\title{Inference using simulation}

\description{
  Implements a collection of methods to perform inferences based on simulation of realizations of the model considered. In particular it implements 
  \dQuote{summary likelihood}, an approach that effectively evaluates and uses the likelihood of simulated summary statistics. 

The procedures for the \dQuote{primitive} workflow, implemented in the first published version of the package, are being maintained for back compatibility, but users are urged to use the distinct, up-to-date workflow (see Examples). The package is expected to perform best when used in combination with \pkg{Rmixmod} for multivariate Gaussian mixture modeling, although the \pkg{mclust} package can also be used (with less control and perhaps decreased performance). The inference workflow typically (but not necessarily) includes dimension-reduction steps (\dQuote{projections}) for summary statistics, and has been more finely tuned for projections performed using the \pkg{ranger} package (though \pkg{Infusion} should handle other methods). 

The up-to-date workflow builds and updates an object of class \code{"SLik_j"} which is designed to carry all the information required for pursuing the inference. Thus, successive inference steps (in particular, successive calls to the key \code{\link{refine}} function) can be carried on different computers, with the caveats noted below, by transferring the object between computers. 

The \code{"SLik_j"} object includes in particular the reference table, information about all projections (this can be memory-expensive, so it has been made possible to remove non-essential information, specifically for \pkg{ranger} results, using \code{\link{deforest_projectors}}: this information will be regenerated automatically if needed), and a sample-simulation function. If the latter function is a wrapper for an external simulation program, then this simulation program must be provided in addition to the \code{"SLik_j"} object. If an external program has to be called differently on different computers, information specific to each computer can be provided by the optional \code{control.Simulate} argument of \code{\link{refine}}.

People used to the functional-programming style common in R, where the return value entirely defines the effect of a function call, may be surprized by the distinct style of some of functions of this package. Indeed, the \code{"SLik_j"} objects include environments that may be modified by functions, independently of what these functions return. Notably, the profile plots (see \code{\link{plot1Dprof}}) and summary-likelihood ratio tests (\code{\link{SLRT}}) may update the  summary-likelihood estimates. 
However, the more significant caveat of such usage of environments, from an inferential viewpoint, is that the projectors stored in input and output fit objects of a \code{refine} call are stored in the same environment, and therefore the projectors stored in the input object are modified in light of new simulations (or by any other operation affecting them in the refine), which alters the statistical meaning of any subsequent operation reusing the projectors from the input object. To keep the unmodified version of the projectors, users may need to save a fit object on disk before a refine call. This feature could of course be modified, but is retained for the following reasons: (1) in the routine workflow, only the latest fit object and its projectors will be relevant (but when comparing performance of the inference method for, e.g., reference tables of different sizes, more care is needed); and (2) projectors can be memory-expensive objects, so keeping distinct versions of them in successive fit objects may have substantial drawbacks.  
}

\details{
  The methods implemented in \code{Infusion} by default assume that the summary statistics have densities. Special values of some statistic, having discrete probability mass, can be more or less automatically handled by the up-to-date workflow, which also handles automatically NA values of summary statistics. For the primitive workflow, both of these problems could be handled to a limited extent using the \code{boundaries} attribute of the observed summary statistics (see \code{\link{handling_NAs}} for one use of this attribute).              
}
\note{See examples \code{\link{example_reftable}} for the most complete example using up-to-date workflow, and \code{\link{example_raw_proj}} or \code{\link{example_raw}} for older workflows.}
\examples{
  ## see Note for links to examples.
}
\keyword{ package }
