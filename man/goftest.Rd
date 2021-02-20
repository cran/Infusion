\name{goftest}
\alias{goftest}
\title{
Assessing goodness of fit of inference using simulation
}
\description{
A goodness-of-fit test is performed in the case projected statistics have been used for inference. Otherwise some plots of limited interest are produced.
}
\usage{
goftest(object, nsim = 99L, method = "", stats=NULL, plot. = TRUE, nb_cores = NULL, 
        Simulate = attr(object$logLs, "Simulate"), 
        packages = attr(object$logLs, "packages"), 
        env = attr(object$logLs, "env"), verbose = interactive())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an \code{SLik} or \code{SLikp} object.
}
  \item{nsim}{
  Number of draws of summary statistics.
}
  \item{method}{
For development purposes, not documented.
}
\item{stats}{
Character vector, or NULL: the set of summary statistics to be used to construct the test. If NULL, the statistics used accross all projections are used.
}
\item{plot.}{
Control diagnostic plots. \code{plot.} can be of logical, character or numeric type. If \code{plot.} is \code{FALSE}, no plot is produced. If \code{plot.} is \code{TRUE} (the default), a data frame of up to 8 goodness-of-fit statistics (the statistics denoted \emph{u} in Details) is plotted. If more than eight raw summary statistics (denoted \emph{s} in Details) were used, then only the first eight \emph{u} are retained (the order or \emph{u} being deduced from the order of use of the \emph{s} in the different projections). If \code{plot.} is a \bold{numeric vector}, then \emph{u}\code{[plot.]} are retained (possibly more than 8 statistics, as in the next case). If \code{plot.} is a \bold{character vector}, then it is used to match the names of the \emph{u} statistics (not of \emph{s}) to be retained in the plot; the names of \emph{u} are built from names of \emph{s} by wrapping the latter within \code{"Res(".")"} (see axes labels of default plots for examples of valid names).    
}
  \item{nb_cores, Simulate, packages, env, verbose}{
See same-named \code{\link{add_simulation}} arguments.
}
}
\details{
The test is somewhat heuristic but appears to give reasonable results (the Example shows how this can be verified). It assumes that all summary statistics are reduced to projections predicting all model parameters. It is then conceived as if any projection \emph{p} predicting a parameter were a sufficient statistic for this parameter, given the information contained in the summary statistics \bold{s} (this is certainly the ideal objective of machine-learning regression methods). Then a statistic \emph{u} independent (under the fitted model) from all projections should be a suitable statistic for testing goodness of fit: if the model is correctly specified, the quantile of observed \emph{u}, in the distribution of \emph{u} under the fitted model, should be uniformly distributed over repeated sampling under the data-generating process. The procedure constructs statistics uncorrelated to all \bold{p} (over repeated sampling under the fitted model) and proceeds as if they were independent from \emph{p} (rather than simply uncorrelated). Statistics \emph{u} uncorrelated to \emph{p} are obtained as the residuals of the regression of each summary statistic to all projections, where the regression input is a simulation table of \code{nsim} replicates of \bold{s} under the fitted model, and of their projections \bold{p} (using the \dQuote{projectors} constructed from the full reference table). The latter regression involves one more, small-\code{nsim}, approximation (as it is the sample correlation that is zeroed) but using the residuals is crucially better than using the original summary statistics (as some ABC software may do). An additional feature of the procedure is to construct a single test statistic \emph{t} from joint residuals \bold{u}, by estimating their joint distribution (using Gaussian mixture modelling) and letting \emph{t} be the density of \bold{u} in this distribution.  
}
\value{
A list with currently a single element
\item{pval }{The p-value of the test (NULL if the test is not feasible).}
}
\examples{
\dontrun{
## Long example, despite minimal settings!
## Showing uniform distribution under correctly-specified model

# Normal(mu,sd) model, with inefficient summary statistics:
blurred <- function(mu,s2,sample.size) {
  s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
  s <- exp(s/4)
  return(c(mean=mean(s),var=var(s)))
}

# Construct reference table and projections once for all replicates:
set.seed(123)
parsp_j <- data.frame(mu=runif(5000L,min=2.8,max=5.2),
                      s2=runif(5000L,min=0.4,max=2.4),sample.size=40)
dsimuls <- add_reftable(,Simulate="blurred",par.grid=parsp_j,verbose=FALSE)
mufit <- project("mu",stats=c("mean","var"),data=dsimuls,verbose=FALSE,knotnbr=5000L)
s2fit <- project("s2",stats=c("mean","var"),data=dsimuls,verbose=FALSE,knotnbr=5000L)
dprojectors <- list(MEAN=mufit,VAR=s2fit)
dprojSimuls <- project(dsimuls,projectors=dprojectors,verbose=FALSE)

# Analysis of replicate draws from data-generating process
gof_draws <- replicate(50, {
  cat(".")
  dSobs <- blurred(mu=4,s2=1,sample.size=40) ## stands for the actual data to be analyzed
  ## ----dcflow--------------------------------------------------------------
  dprojSobs <- project(dSobs,projectors=dprojectors)
  dslik <- infer_SLik_joint(dprojSimuls,stat.obs=dprojSobs,verbose=FALSE)
  dslik <- MSL(dslik, verbose=FALSE, eval_RMSEs=FALSE)
  #dslik <- refine(dslik,maxit=2L,verbose=FALSE, update_projectors=TRUE)
  gof <- goftest(dslik,nsim = 99L,nb_cores = 1L, method="", plot.=FALSE,verbose=FALSE)
  cat(unlist(gof))
  gof
})
plot(ecdf(unlist(gof_draws)))
}
}
