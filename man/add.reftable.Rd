\name{add_reftable}
\alias{add_reftable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create or augment a list of simulated distributions of summary statistics 
}
\description{
\code{add_reftable} creates or augments a reference table of simulations, and formats the results appropriately for further use. The user does not have to think about this return format. Instead, s-he only has to think about the very simple return format of the function given as its \code{Simulate} argument. The primary role of his function is to wrap the call(s) of the function specified by \code{Simulate}. Depending on the arguments, parallel or serial computation is performed. 

When parallelization is implied, it is performed by by default a \dQuote{socket} cluster, available on all operating systems. Special care is then needed to ensure that all required packages are loaded in the called processes, and that all required variables and functions are passed therein: check the \code{packages} and \code{env} arguments. For socket clusters, \code{foreach} or \code{pbapply} is called depending whether the \code{doSNOW} package is attached (\code{doSNOW} allows more efficient load balancing than \code{pbapply}). 

Alternatively, if the simulation function cannot be called directly by the R code, simulated samples can be added using the \code{newsimuls} argument. Finally, a generic data frame of simulated samples can be reformatted as a reference table by using only the \code{reftable} argument.

\code{add_simulation} is a wrapper for \code{add_reftable}, suitable when \code{nRealizations}>1. It is now distinctly documented: the distinct features of \code{add_simulation} were conceived for the first workflow implemented in \code{Infusion} but are somewhat obsolete now.

}
\usage{
add_reftable(reftable=NULL, Simulate, parsTable=par.grid, par.grid=NULL, 
               nRealizations = 1L, newsimuls = NULL, 
               verbose = interactive(), nb_cores = NULL, packages = NULL, env = NULL,
               control.Simulate=NULL, cluster_args=list(), cl_seed=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reftable}{
Data frame: a reference table. Each row contains parameters value of a simulated realization of the data-generating process, and the simulated summary statistics.  
As parameters should be told apart from statistics by \pkg{Infusion} functions, information about parameter names should be attached to the \code{reftable} \bold{*if*} it is not available otherwise. Thus if no \code{parsTable} is provided, the \code{reftable} should have an attribute \code{"LOWER"} (a named vectors giving lower bounds for the parameters which will vary in the analysis, as in the return value of the function).  
}
  \item{Simulate}{
An *R* function, or the name (as a character string) of an *R* function used to generate summary statistics for samples form a data-generating process. When an external simulation program is called, \code{Simulate} must therefore be an R function wrapping the call to the external program. Two function APIs are handled:\cr \code{ * }If the function has a \code{parsTable} argument, it must return a data frame of summary statistics, each line of which contains the vector of summary statistics for one realization of the data-generating process. The \code{parsTable} argument of \code{add_reftable} will be passed to \code{Simulate} and lines of the output data frame must be ordered, as in the input \code{parsTable} as these two data frames will be bound together.\cr \code{ * }Otherwise, the \code{Simulate} function must have one argument for each element of the parameter vector (i.e. of each row of \code{parsTable}). It must return a vector of summary statistics with named vector member.
}
  \item{parsTable, par.grid}{
 A data frame of which each line is the vector of parameters needed by \code{Simulate} for each simulation of the data-generating process. \code{par.grid} is an alias for \code{parsTable}; the latter argument may be preferred in order not to suggest that the parameter values should form a regular grid.
}
\item{nRealizations}{
  The number of simulated samples of summary statistics, for each parameter vector (each row of \code{parsTable}). If not 1, theold wrkflow is assumed and \code{\link{add_simulation}} is called.
}
  \item{newsimuls}{
  If the function used to generate empirical distributions cannot be called by R, then \code{newsimuls} can be used to provide these distributions. See Details for the structure of this argument.
}
  \item{nb_cores}{Number of cores for parallel simulation; \code{NULL} or integer value, acting as a shortcut for \code{cluster_args$spec}. This is effective only if the simulation function is called separately for each row of \code{parsTable}. Otherwise, if the simulation function is called once one the whole \code{parsTable}, parallelisation could be controlled only through that function's own arguments.   
 }
 \item{cluster_args}{A list of arguments, passed to \code{\link[parallel]{makeCluster}}. May contain a non-null \code{spec} element, in which case the distinct \code{nb_cores} argument and the global \pkg{Infusion} option \code{nb_cores} are ignored. A typical usage would thus be \code{control_args=list(spec=<number of 'children'>)}. Additional elements \code{outfile="log.txt"} may be useful to collect output from the nodes, and \code{type="FORK"} may be used to force a fork cluster on linux(-alikes) (otherwise a socket cluster is set up as this is the default effect of \code{parallel::makeCluster}). Do \bold{*not*} use a structured list with an \code{add_reftable} element as is possible for \code{refine} (see Details of \code{\link{refine}} documentation).
}
  \item{verbose}{
  Whether to print some information or not.
 }
  \item{\ldots}{Additional arguments passed to \code{Simulate}, beyond the parameter vector. These arguments should be constant through all the simulation workflow. 
 }
 \item{control.Simulate}{A list, used as an exclusive alternative to \dQuote{\ldots} to pass additional arguments to \code{Simulate}, beyond the parameter vector. The list must contain the same elements as would otherwise go in the \dQuote{\ldots} (if \code{control.Simulate} is left NULL, a default value is constructed from the \ldots).}
   \item{packages}{For parallel evaluation: Names of additional libraries to be loaded on the cores, necessary for \code{Simulate} evaluation.}
   \item{env}{For parallel evaluation: an environment containing additional objects to be exported on the cores, necessary for \code{Simulate} evaluation.}
     \item{cl_seed}{
(all parallel contexts:) Integer, or NULL. If an integer, it is used to initialize \code{"L'Ecuyer-CMRG"} random-number generator. If \code{cl_seed} is \code{NULL}, the default generator is selected on each node, where its seed is not controlled. Providing the seed allows repeatable results for given parallelization settings, but may not allow identical results across different settings. 
}
}
\details{

The \code{newsimuls} argument should have the same structure as the return value of the function itself, except that \code{newsimuls} may include only a subset of the attributes returned by the function. It is thus a data frame; its required attributes are \code{LOWER} and \code{UPPER} which are named vectors giving bounds for the parameters which are variable in the whole analysis (note that the names identify these parameters in the case this information is not available otherwise from the arguments). The values in these vectors may be incorrect in the sense of failing to bound the parameters in the \code{newsimuls}, as the actual bounds are then corrected using parameter values in \code{newsimuls} and attributes from \code{reftable}.  

}
\value{
A data.frame (with additional attributes) is returned. 

The value has the following attributes: \code{LOWER} and \code{UPPER} which are each a vector of per-parameter minima and maxima deduced from any \code{newsimuls} argument, and optionally any of the arguments \code{Simulate, control.Simulate, packages, env, parsTable} and  \code{reftable} (all corresponding to input arguments when provided, except that the actual \code{Simulate} function is returned even if it was input as a name).
}
\examples{
## see main documentation page for the package for other typical usage
}

