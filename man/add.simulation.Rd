\name{add_simulation}
\alias{add_simulation}
\alias{add_reftable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create or augment a list of simulated distributions of summary statistics 
}
\description{
\code{add_simulation} creates or augments a list of simulated distributions of summary statistics, and formats the simulation results appropriately for further use. The user does not have to think about this return format. Instead, s-he only has to think about the very simple return format of the function given as its \code{Simulate} argument. Alternatively, if the simulation function cannot be called directly by the R code, simulated distributions can be added easily using the \code{newsimuls} argument, again using a simple format (see \code{onedistrib} in the Examples). \code{add_reftable} is a wrapper for \code{add_simulation}, enforcing \code{nRealizations=1}.
}
\usage{
add_simulation(simulations=NULL, Simulate, par.grid=NULL, 
               nRealizations = NULL,
               newsimuls = NULL, verbose = interactive())
add_reftable(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{simulations}{
A list of simulations
}
  \item{Simulate}{
The name (as a character string) of the *R* function used to generate empirical distributions of summary statistics (notably, an for external simulation program, an R function wrapping the call to the external program must be provided). The function must have a single vector as argument, matching each row of \code{par.grid}, and must return a vector of summary statistics with named vector members.
}
  \item{par.grid}{
 A data frame of which each line matches the single vector argument of \code{Simulate}.
}
  \item{nRealizations}{
  The number of simulated samples of summary statistics, for each empirical distribution (each row of \code{par.grid}). If the argument is NULL, the value is obtained by \code{Infusion.getOption}. If the argument is not NULL, \code{Infusion.options(nRealizations)} is set, but restored, on exit from \code{add_simulation}, to its initial value.
}
  \item{newsimuls}{
  If the function used to generate empirical distributions cannot be called by R, then \code{newsimuls} can be used to provide these distributions. This argument should be list of matrices, each with a \code{par} attribute (see Examples). Rows of each matrix stand for simulation replicates and rows for the different summary statistics. This is the same format as the return value of \code{add_simulation}, except that the latter adds further attributes.
}
  \item{verbose}{
  Whether to print some information or not
 }
  \item{...}{
  Arguments passed to \code{add_simulation}. Any of the above aguments is valid, except \code{nRealizations}.
 }
}
\value{
If only one realization is computed for each (vector-valued) parameter, a data.frame is returned. 
Otherwise, the return value is an objet of class \code{EDFlist}, which is a list of matrices-with-attribute. Each matrix contains a simulated distribution of summary statistics for given parameters, and the \code{"par"} attribute is a vector of parameters. If \code{Simulate} is used, this must give the full parameters required for a call of the simulation function; otherwise it must at least include all variable parameters in this \bold{or later} simulations to be appended to the simulation list. 
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
add_simulation(NULL, Simulate="myrnorm",
                     newsimuls=list("example"=onedistrib))

## see main documentation page for the package for other typical usage
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
