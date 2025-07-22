\name{def_projectors}
\alias{def_projectors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper to generate projection functions for all parameters
}
\description{
Convenience function automating the initial \code{\link{project.character}} calls of a typical workflow. 
}
\usage{
def_projectors(reftable, pars, stats, npp = 1L, projNames = NULL, 
               npp_opt=Infusion.getOption("npp_opt"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reftable}{
A reference table of simulation as produced by \code{\link{add_reftable}}.
}
  \item{pars}{
Fitted parameter names for which projection functions should be returned.
}
  \item{stats}{
Names of summary statistics to be used for predicting parameter values.
}
  \item{npp, npp_opt}{
  For development purposes.
% A single positive integer, or a vector of such integers with one element for each parameter.
% Specifies the number of projections per parameter. The only non-default value tried so far is 2 (\code{2L}).  
}
  \item{projNames}{
Character vector of names of the projections. By default, names are given as \code{paste0("p", pars)} 
(or a more elaborate default when \code{npp > 1L}), and these names are checked for possible name conflicts 
with pre-existing variables in the \code{reftable}. Non-default values should be given in case of conflict. 
}
\item{\ldots}{Other arguments passed to \code{\link{project.character}}.}.
}
\details{
The \code{npp > 1L} feature is experimental and awaiting further documentation. % ____F I X M E_____
}
\value{
A list of \emph{projectors}. Each of them is either an object returned by \code{project.character}, 
or a \code{list} inheriting from class \code{"projector"} which includes such a return object as one of its elements. 
The more detailed structure of a \code{"projector"} object is not part of the API.
}
\examples{
## variation on example from help("example_reftable")

if (Infusion.getOption("example_maxtime")>2) {
  blurred <- function(mu,s2,sample.size) {
    s <- rnorm(n=sample.size,mean=mu,sd=sqrt(s2))
    s <- exp(s/4)
    return(c(mean=mean(s),var=var(s)))
  }
  
  set.seed(123)
  dSobs <- blurred(mu=4,s2=1,sample.size=40)
  
  parsp_j <- init_reftable(lower=c(mu=2.5, s2=0.25, sample.size=40), 
                           upper=c(mu=5.2, s2=2.4, sample.size=40))
  dsimuls <- add_reftable(,Simulate="blurred", parsTable=parsp_j,verbose=FALSE)

  ## Then instead of
  #
  # mufit <- project("mu",stats=c("mean","var"),data=dsimuls,verbose=FALSE)
  # s2fit <- project("s2",stats=c("mean","var"),data=dsimuls,verbose=FALSE)
  # dprojectors <- list(MEAN=mufit,VAR=s2fit)
  #
  ## one can use
  #
  dprojectors <- def_projectors(dsimuls, pars=c("mu","s2"), 
                   stats=c("mean","var"), npp = 1L, 
                   projNames = c("MEAN", "VAR"))
                   
  
  ## More speculatively, for non-default 'npp':                   
  dprojectors <- def_projectors(dsimuls, pars=c("mu","s2"), 
                   stats=c("mean","var"), npp = 2L)
}
}
