\name{add_simulation}
\alias{add_simulation}
\alias{add_reftable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create or augment a list of simulated distributions of summary statistics 
}
\description{
\code{add_simulation} creates or augments a list of simulated distributions of summary statistics, and formats the simulation results appropriately for further use. The user does not have to think about this return format. Instead, s-he only has to think about the very simple return format of the function given as its \code{Simulate} argument. Alternatively, if the simulation function cannot be called directly by the R code, simulated distributions can be added easily using the \code{newsimuls} argument, again using a simple format (see \code{onedistrib} in the Examples). 

\code{add_reftable} is a wrapper for \code{add_simulation}, enforcing \code{nRealizations=1}: see example \code{\link{example_reftable}}.

These functions can run simulations in a parallel environment. Special care is then needed to ensure that all required packages are loaded in the called processes, and required all variables and function are passed therein: check the \code{packages} and \code{env} arguments.
}
\usage{
add_simulation(simulations=NULL, Simulate, par.grid=NULL, 
               nRealizations = NULL, newsimuls = NULL, 
               verbose = interactive(), nb_cores = NULL, packages = NULL, env = NULL,
               control.Simulate=NULL, ...)
add_reftable(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulations}{
A list of simulations
}
  \item{nRealizations}{
  The number of simulated samples of summary statistics, for each empirical distribution (each row of \code{par.grid}). If the argument is NULL, the value is obtained by \code{Infusion.getOption}. If the argument is not NULL, \code{Infusion.options(nRealizations)} is set, but restored, on exit from \code{add_simulation}, to its initial value.
}
  \item{Simulate}{
The name (as a character string) of the *R* function used to generate empirical distributions of summary statistics (notably, for external simulation program, an R function wrapping the call to the external program must be provided). The function must have a single vector as argument, matching each row of \code{par.grid}. It must return a vector of summary statistics with named vector members; \bold{or} a single matrix of \code{nRealizations} simulations, in which case its rows and row names, must represent the summary statistics, it should have \code{nRealizations} columns, and \code{nRealizations} should be an integer names as \dQuote{\code{as_one}} (see Examples).
}
  \item{par.grid}{
 A data frame of which each line matches the single vector argument of \code{Simulate}.
}
  \item{newsimuls}{
  If the function used to generate empirical distributions cannot be called by R, then \code{newsimuls} can be used to provide these distributions. See Details for the structure of this argument.
}
  \item{nb_cores}{Number of cores for parallel simulation; \code{NULL} or integer value. If the latter is unnamed or has name \code{"replic"}, parallelisation is over the different samples for each parameter value. For any other name (e.g., \code{nb_cores=c(foo=7)}), parallelisation is over the parameter values (the rows of \code{par.grid}). In all cases, the progress bar is over parameter values. See Details in \code{\link{Infusion.options}} for the subtle way these different cases are distinguished in the progress bar.
 }
  \item{verbose}{
  Whether to print some information or not.
 }
  \item{\ldots}{
  For \code{add_reftable}: arguments passed to \code{add_simulation}. Any of the \code{add_simulation} arguments is valid, except \code{nRealizations}.    For \code{add_simulation}: additional arguments passed to \code{Simulate}, beyond the parameter vector; see \code{nsim} argument of \code{myrnorm_tab()} in the Examples. These arguments should be constant through all the simulation workflow.
 }
 \item{control.Simulate}{A list, used as an exclusive alternative to \dQuote{\ldots} to pass additional arguments to \code{Simulate}, beyond the parameter vector. The list must contain the same elements as would go in the \dQuote{\ldots}.}
   \item{packages}{For parallel evaluation: Names of additional libraries to be loaded on the cores, necessary for \code{Simulate} evaluation.}
   \item{env}{For parallel evaluation: an environment containing additional objects to be exported on the cores, necessary for \code{Simulate} evaluation.}
}
\details{
The \code{newsimuls} argument should have the same structure as the return value of \code{add_simulation} itself, except that \code{newsimuls} may include only a subset of the attributes returned by \code{add_simulation}. \bold{In the reference-table case}, it is thus a data frame; its required attributes are \code{LOWER} and \code{UPPER} which are named vectors giving bounds for the parameters which are variable in the whole analysis (note that the names identify these parameters in the case this information is not available otherwise from the arguments). The values in these vectors may be incorrect in the sense of failing to bound the parameters in the \code{newsimuls}, as the actual bounds are then corrected using parameter values in \code{newsimuls} and attributes from \code{simulations}. \bold{Otherwise}, \code{newsimuls} should be list of matrices, each with a \code{par} attribute (see Examples). Rows of each matrix stand for simulation replicates and columns stand for the different summary statistics. 
}
\value{
If only one realization is computed for each (vector-valued) parameter, a data.frame (with additional attributes) is returned. 
Otherwise, the return value is an objet of class \code{EDFlist}, which is a list-with-attributes of matrices-with-attribute. Each matrix contains a simulated distribution of summary statistics for given parameters, and the \code{"par"} attribute is a 1-row data.frame of parameters. If \code{Simulate} is used, this must give all the parameters to be estimated; otherwise it must at least include all variable parameters in this \bold{or later} simulations to be appended to the simulation list. 

The value has the following attributes: \code{Simulate, control.Simulate, packages, env} (all orresponding to input arguments), \code{LOWER} and \code{UPPER} which are each a vector of per-parameter minima and maxima deduced from any \code{newsimuls} argument, any \code{par.grid} argument, and any \code{simulations} argument. 
}
\examples{
# example of building a list of simulations from scratch:
myrnorm <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  return(c(mean=mean(s),var=var(s)))
}
set.seed(123)
onedistrib <- t(replicate(100,myrnorm(1,1,10))) # toy example of simulated distribution
attr(onedistrib,"par") <- c(mu=1,sigma=1,sample.size=10) ## important!
simuls <- add_simulation(NULL, Simulate="myrnorm",
                         newsimuls=list("example"=onedistrib))
                     
## Example where a single Simulate returns all nRealizations:
myrnorm_tab <- function(mu,s2,sample.size, nsim) {
  replicate(nsim, 
            myrnorm(mu=mu,s2=s2,sample.size=sample.size)) ## n_summ_stats * nReal 
}
op <- Infusion.options(nRealizations=c(as_one=500)) 
parsp <- init_grid(lower=c(mu=2.8,s2=0.2,sample.size=40),
                   upper=c(mu=5.2,s2=3,sample.size=40))
simuls <- add_simulation(NULL,Simulate="myrnorm_tab",
                         nsim=Infusion.getOption('nRealizations'),
                         par.grid=parsp)
Infusion.options(op)

\dontrun{
# parallel version of the same:
op <- Infusion.options(nb_cores=2, nRealizations=c(as_one=500))
myrnorm_tab <- function(mu,s2,sample.size) {
  ## By default, Infusion.getOption('nRealizations') will fail on child process!
  replicate(1000,   
            myrnorm(mu=mu,s2=s2,sample.size=sample.size)) 
}
simuls <- add_simulation(NULL,Simulate="myrnorm_tab",par.grid=parsp,
                         # need to export function called by Simulate:
                         env=list2env(list(myrnorm=myrnorm)))
Infusion.options(op)
}

## see main documentation page for the package for other typical usage
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
