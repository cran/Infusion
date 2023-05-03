\name{add_simulation}
\alias{add_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create or augment a list of simulated distributions of summary statistics 
}
\description{
\code{add_simulation} is suitable for the primitive \pkg{Infusion} workflow; otherwise, it is cleaer to call \code{\link{add_reftable}} directly. 
\code{add_simulation} creates or augments a list of simulated distributions of summary statistics, and formats the results appropriately for further use. Alternatively, if the simulation function cannot be called directly by the R code, simulated distributions can be added using the \code{newsimuls} argument, using a simple format (see \code{onedistrib} in the Examples). Finally, a generic data frame of simulations can be reformatted as a reference table by using only the \code{simulations} argument.

Depending on the arguments, parallel or serial computation is performed. When parallelization is implied, by default a \dQuote{socket} cluster, available on all operating systems. Special care is then needed to ensure that all required packages are loaded in the called processes, and that all required variables and functions are passed therein: check the \code{packages} and \code{env} arguments. For socket clusters, \code{foreach} or \code{pbapply} is called depending whether the \code{doSNOW} package is attached (\code{doSNOW} allows more efficient load balancing than \code{pbapply}). 
}
\usage{
add_simulation(simulations=NULL, Simulate, parsTable=par.grid, par.grid=NULL, 
               nRealizations=Infusion.getOption("nRealizations"),
               newsimuls=NULL, verbose=interactive(), nb_cores=NULL, 
               packages=NULL, env=NULL, control.Simulate=NULL,
               cluster_args=list(), cl_seed=NULL, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulations}{
A list of matrices each representing a simulated distribution for given parameters in a format consistent with the return format of \code{add_simulation}.
}
  \item{nRealizations}{
  The number of simulated samples of summary statistics, for each empirical distribution (each row of \code{par.grid}).
}
  \item{Simulate}{
An *R* function, or the name (as a character string) of an *R* function used to generate empirical distributions of summary statistics. When an external simulation program is called, \code{Simulate} must therefore be an R function wrapping the call to the external program. The \code{Simulate} function must have one argument for each element of the parameter vector (i.e. of each row of \code{par.grid}). It must return a vector of summary statistics with named vector members; \bold{or} a single matrix of \code{nRealizations} simulations, in which case its rows and row names must represent the summary statistics, it should have \code{nRealizations} columns, and \code{nRealizations} should be named integer of the form \dQuote{\code{c(as_one=.)}} (see Examples).
}
  \item{parsTable, par.grid}{
 A data frame of which each line is the vector of parameters needed by \code{Simulate} for each simulation of the data-generating process. \code{par.grid} is an alias for \code{parsTable}; the latter argument may be preferred in order not to suggest that the parameter values should form a regular grid.
}
  \item{newsimuls}{
  If the function used to generate empirical distributions cannot be called by R, then \code{newsimuls} can be used to provide these distributions. See Details for the structure of this argument.
}
  \item{nb_cores}{Number of cores for parallel simulation; \code{NULL} or integer value, acting as a shortcut for \code{cluster_args$spec}. The effect is complicated: see Details.
 }
 \item{cluster_args}{A list of arguments, passed to \code{\link[parallel]{makeCluster}}. May contain a non-null \code{spec} element, in which case the distinct \code{nb_cores} argument and the global \pkg{Infusion} option \code{nb_cores} are ignored. A typical usage would thus be \code{control_args=list(spec=<number of 'children'>)}. Additional elements \code{outfile="log.txt"} may be useful to collect output from the nodes, and \code{type="FORK"} may be used to force a fork cluster on linux(-alikes) (otherwise a socket cluster is set up as this is the default effect of \code{parallel::makeCluster}).
}
  \item{verbose}{
  Whether to print some information or not.
 }
  \item{\ldots}{
  Arguments passed to \code{add_reftable} (and possibly beyond, to the simulation function: see \code{nsim} argument of \code{myrnorm_tab()} in the Examples. These arguments should be constant through all the simulation workflow.
 }
 \item{control.Simulate}{A list, used as an exclusive alternative to \dQuote{\ldots} to pass additional arguments to \code{Simulate}, beyond the parameter vector. The list must contain the same elements as would go in the \dQuote{\ldots}.}
   \item{packages}{For parallel evaluation: Names of additional libraries to be loaded on the cores, necessary for \code{Simulate} evaluation.}
   \item{env}{For parallel evaluation: an environment containing additional objects to be exported on the cores, necessary for \code{Simulate} evaluation.}
     \item{cl_seed}{Integer, or NULL. Providing the seed was conceived to allow repeatable results at least for given parallelization settings, if not identical results across different parallelization contexts. However, this functionality may have been been lost as the code was adapted for the up-to-date workflow using \code{add_reftable}.
}
}
\details{

The \code{newsimuls} argument should have the same structure as the return value of the function itself, except that \code{newsimuls} may include only a subset of the attributes returned by the function. \code{newsimuls} should thus be list of matrices, each with a \code{par} attribute (see Examples). Rows of each matrix stand for simulation replicates and columns stand for the different summary statistics. 

When \code{nRealizations}>1L, if \code{nb_cores} is unnamed or has name \code{"replic"} and if the simulation function does not return a single table for all replicates (thus, if \code{nRealizations} is \bold{not} a named integer of the form \dQuote{\code{c(as_one=.)}}, parallelisation is over the different samples for each parameter value (and the seed of the random number generator is not controlled in a parallel context). For any other explicit name (e.g., \code{nb_cores=c(foo=7)}), or if \code{nRealizations} is a named integer of the form \dQuote{\code{c(as_one=.)}}, parallelisation is over the parameter values (the rows of \code{par.grid}). In all cases, the progress bar is over parameter values. See Details in \code{\link{Infusion.options}} for the subtle way these different cases are distinguished in the progress bar.

Using a FORK cluster with \code{nRealizations}>1 is warned as unreliable: in particular, anyone trying this combination should check whether other desired controls, such as random generator seed, or progress bar are effective. 
}
\value{
If \code{nRealizations}>1L, the return value is an object of class \code{EDFlist}, which is a list-with-attributes of matrices-with-attribute. Each matrix contains a simulated distribution of summary statistics for given parameters, and the \code{"par"} attribute is a 1-row data.frame of parameters. If \code{Simulate} is used, this must give all the parameters to be estimated; otherwise it must at least include all variable parameters in this \bold{or later} simulations to be appended to the simulation list. 

The value has the following attributes: \code{LOWER} and \code{UPPER} which are each a vector of per-parameter minima and maxima deduced from any \code{newsimuls} argument, and optionally any of the arguments \code{Simulate, control.Simulate, packages, env, par.grid} and  \code{simulations} (all corresponding to input arguments when provided, except that the actual \code{Simulate} function is returned even if it was input as a name).

If \code{nRealizations}=1 \code{add_reftable} is called: see its distinct return value.
}
\examples{
### Examples using init_grid and add_simulation, for primitive workflow
### Use init_reftable and add_reftable for the up-to-date workflow

# example of building a list of simulations from scratch:
myrnorm <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  return(c(mean=mean(s),var=var(s)))
}
set.seed(123)
onedistrib <- t(replicate(100,myrnorm(1,1,10))) # toy example of simulated distribution
attr(onedistrib,"par") <- c(mu=1,sigma=1,sample.size=10) ## important!
simuls <- add_simulation(NULL, Simulate="myrnorm", nRealizations=500,
                         newsimuls=list("example"=onedistrib))

# standard use: smulation over a grid of parameter values
parsp <- init_grid(lower=c(mu=2.8,s2=0.2,sample.size=40),
                   upper=c(mu=5.2,s2=3,sample.size=40))
simuls <- add_simulation(NULL, Simulate="myrnorm", nRealizations=500,
                         par.grid = parsp[1:7,])
                         
\dontrun{ # example continued: parallel versions of the same
# Slow computations, notably because cluster setup is slow.

#    ... parallel over replicates, serial over par.grid rows
# => cl_seed has no effect and can be ignored
simuls <- add_simulation(NULL, Simulate="myrnorm", nRealizations=500,
                         par.grid = parsp[1:7,], nb_cores=7)
#                         
#    ... parallel over 'par.grid' rows => cl_seed is effective
simuls <- add_simulation(NULL, Simulate="myrnorm", nRealizations=500,
                         cl_seed=123, # for repeatable results
                         par.grid = parsp[1:7,], nb_cores=c(foo=7))
}
                     
####### Example where a single 'Simulate' returns all replicates:

myrnorm_tab <- function(mu,s2,sample.size, nsim) {
  ## By default, Infusion.getOption('nRealizations') would fail on nodes!
  replicate(nsim, 
            myrnorm(mu=mu,s2=s2,sample.size=sample.size)) 
}

parsp <- init_grid(lower=c(mu=2.8,s2=0.2,sample.size=40),
                   upper=c(mu=5.2,s2=3,sample.size=40))

# 'as_one' syntax for 'Simulate' function returning a simulation table: 
simuls <- add_simulation(NULL, Simulate="myrnorm_tab",
              nRealizations=c(as_one=500),
              nsim=500, # myrnorm_tab() argument, part of the 'dots'
              par.grid=parsp)

\dontrun{ # example continued: parallel versions of the same.
# Slow cluster setup again
simuls <- add_simulation(NULL,Simulate="myrnorm_tab",par.grid=parsp,
              nb_cores=7L,
              nRealizations=c(as_one=500),
              nsim=500, # myrnorm_tab() argument again
              cl_seed=123, # for repeatable results
              # need to export other variables used by *myrnorm_tab* to the nodes:
              env=list2env(list(myrnorm=myrnorm)))
}

## see main documentation page for the package for other typical usage
}

