\name{goftest}
\alias{goftest}
\alias{summary.goftest}
\alias{print.goftest}
\title{
Assessing goodness of fit of inference using simulation
}
\description{
A goodness-of-fit test is performed in the case projected statistics have been used for inference. Otherwise some plots of limited interest are produced.

\code{summary} and \code{print} methods for results of \code{goftest} call \code{str} to display the structure of this result.
}
\usage{
goftest(object, nsim = 99L, method = "", stats=NULL, plot. = TRUE, nb_cores = NULL, 
        Simulate = get_from(object,"Simulate"), 
        control.Simulate=get_from(object,"control.Simulate"),
        packages = get_from(object,"packages"), 
        env = get_from(object,"env"), verbose = interactive(),
        cl_seed=.update_seed(object), get_gof_stats=.get_gof_stats)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} or \code{SLik_j} object.
}
  \item{nsim}{
  Number of draws of summary statistics.
}
  \item{method}{
For development purposes, not documented.
}
\item{stats}{
Character vector, or NULL: the set of summary statistics to be used to construct the test. If NULL, the union, across all projections, of the raw summary statistics used for projections is potentially used for goodness of fit; however, if this set is too large for gaussian mixture modelling, a subset of variable may be selected. How they are selected is not yet fully settled (see Details).  
}
\item{plot.}{
Control diagnostic plots. \code{plot.} can be of logical, character or numeric type. If \code{plot.} is \code{FALSE}, no plot is produced. If \code{plot.} is \code{TRUE} (the default), a data frame of up to 8 goodness-of-fit statistics (the statistics denoted \emph{u} in Details) is plotted. If more than eight raw summary statistics (denoted \emph{s} in Details) were used, then only the first eight \emph{u} are retained (see Details for the ordering of the \emph{u}s here). 
If \code{plot.} is a \bold{numeric vector}, then \emph{u}\code{[plot.]} are retained (possibly more than 8 statistics, as in the next case). If \code{plot.} is a \bold{character vector}, then it is used to match the names of the \emph{u} statistics (not of \emph{s}) to be retained in the plot; the names of \emph{u} are built from names of \emph{s} by wrapping the latter within \code{"Res(".")"} (see axes labels of default plots for examples of valid names).    
}
  \item{nb_cores, Simulate, packages, env, verbose}{
See same-named \code{\link{add_simulation}} arguments.
}
  \item{control.Simulate}{
    A list of arguments of the \code{Simulate} function (see\code{\link{add_simulation}}). The default value should generally be used, unless e.g. it contains the path of an executable on one machine and a different path must be specified on another machine.  
}
\item{cl_seed}{NULL or integer (see \code{\link{refine}} for Details).}
\item{get_gof_stats}{function for selecting raw statistics (see Details).}
}
\details{
\bold{Testing goodness-of-fit:} The test is somewhat heuristic but appears to give reasonable results (the Example shows how this can be verified). It assumes that all summary statistics are reduced to projections predicting all model parameters. It is then conceived as if any projection \emph{p} predicting a parameter were a sufficient statistic for this parameter, given the information contained in the summary statistics \bold{s} (this is certainly the ideal objective of machine-learning regression methods). Then a statistic \emph{u} independent (under the fitted model) from all projections should be a suitable statistic for testing goodness of fit: if the model is correctly specified, the quantile of observed \emph{u}, in the distribution of \emph{u} under the fitted model, should be uniformly distributed over repeated sampling under the data-generating process. The procedure constructs statistics uncorrelated to all \bold{p} (over repeated sampling under the fitted model) and proceeds as if they were independent from \emph{p} (rather than simply uncorrelated). A number (depending on the size of the reference table) of statistics \emph{u} uncorrelated to \emph{p} are then defined. Each such statistic is obtained as the residual of the regression of a given raw summary statistic to all projections, where the regression input is a simulation table of \code{nsim} replicates of \bold{s} under the fitted model, and of their projections \bold{p} (using the \dQuote{projectors} constructed from the full reference table). The latter regression involves one more, small-\code{nsim}, approximation (as it is the sample correlation that is zeroed) but using the residuals is crucially better than using the original summary statistics (as some ABC software may do). An additional feature of the procedure is to construct a single test statistic \emph{t} from joint residuals \bold{u}, by estimating their joint distribution (using Gaussian mixture modelling) and letting \emph{t} be the density of \bold{u} in this distribution.  

\bold{Selection of raw summary statistics:} See the code of the \code{Infusion:::..get_gof_stats} function for the method used. It requires that \code{ranger} has been used to produce the projectors, and that the latter include variable importance statistics (by default, \pkg{Infusion} calls \code{ranger} with argument \code{importance="permutation"}). \code{.get_gof_stats} then selects the raw summary statistics with \emph{least} importance over projections (this may not be optimal, and in particular appears redundant with the procedure described below to construct goodness-of-fit statistics from raw summary statistics; so this might change in a later version), and returns a vector of names of raw statistics, sorted by increasing least-importance. The number of summary statistics can be controlled by the global package option \code{gof_nstats_fn}, a function with arguments \code{nr} and \code{nstats} for, respectively, the number of simulations of the processus (as controlled by \code{goftest(.,nsim)}) and the total number of raw summary statistics used in the projections.

The \bold{diagnostic plot} will show a data frame of residuals \emph{u} of the summary statistics identified as the first elements of the vector returned by \code{Infusion:::..get_gof_stats}, i.e. again a set of raw statistics with least-importance over projectors. 

}
\value{
An object of class \code{goftest}, which is a\code{list} with element(s)
\item{pval }{The p-value of the test (NULL if the test is not feasible).}
\item{plotframe }{The data frame which is (by default) plotted by the function.
Its last line contains the residuals \emph{u} for the analyzed data, and other lines contain the bootstrap replicates.}
}
\examples{
### See end of example("example_reftable") for minimal example.

\dontrun{
### Performance of GoF test over replicate draws from data-generating process

# First, run 
example("example_reftable") 
# (at least up to the final 'slik_j' object), then

# as a shortcut, the same projections will be used in all replicates:
dprojectors <- slik_j$projectors 

set.seed(123)
gof_draws <- replicate(200, {
  cat(" ")
  dSobs <- blurred(mu=4,s2=1,sample.size=40) 
  ## ----Inference workflow-----------------------------------------------
  dprojSobs <- project(dSobs,projectors=dprojectors)
  dslik <- infer_SLik_joint(dprojSimuls,stat.obs=dprojSobs,verbose=FALSE)
  dslik <- MSL(dslik, verbose=FALSE, eval_RMSEs=FALSE)
  ## ----GoF test-----------------------------------------------
  gof <- goftest(dslik,nb_cores = 1L, plot.=FALSE,verbose=FALSE) 
  cat(unlist(gof))
  gof
})
# ~ uniform distribution under correctly-specified model: 
plot(ecdf(unlist(gof_draws)))
}

}
