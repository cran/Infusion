\name{get_nbCluster_range}
\alias{get_nbCluster_range}
\alias{refine_nbCluster}
\alias{seq_nbCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Control of number of components in Gaussian mixture modelling
}
\description{
These functions implement the default values for the number of components tried in Gaussian mixture modelling (matching the \code{nbCluster} argument of \code{Rmixmod::mixmodCluster()}). \code{get_nbCluster_range} allows the user to reproduce the internal rules used by \pkg{Infusion} to determine this argument. \code{seq_nbCluster} is a wrapper to the function defined by the \code{nbClu_pow_rule_fn} global option of the package. Its default result is a sequence of integers determined by the number of rows of the data (see \code{\link{Infusion.options}}). \code{get_nbCluster_range()} uses additional criteria involving the number of columns of the data to determine the maximum number of clusters. This maximum is controlled by the function defined by the \code{maxnbCluster} global option of the package. 

\code{refine_nbCluster} controls the default number of clusters of \code{\link{refine}}: it gets the range from \code{seq_nbCluster} and keeps only the maximum value of this range if this maximum is higher than the \code{onlymax} argument.

Adventurous users can change the rules used by \pkg{Infusion} by changing the global options \code{nbClu_pow_rule_fn} and \code{maxnbCluster} (while conforming to the interfaces of these functions). Less ambitiously, they can for example use the maximum value of the result of \code{get_nbCluster_range()} as a single reasonable value for the \code{nbCluster} argument of \code{infer_SLik_joint}.
}
\details{
The default upper value of the \code{nbCluster} range is controlled by two rules:
  
\code{ * }The first rule sets the maximum number of clusters as function of the number of samples \eqn{n} in the reference table. The default rule, as implemented in \code{Infusion:::.nbClu_pow_rule}, reduces when more than three parameters are estimated to the rule \eqn{n^{0.31}} discussed in Supplementary Information of Rousset et al. (2026).   

\code{ * }This first rule is corrected by a second rule setting a maximum dependent also on the dimensions of the \code{projdata} (the one used internally for clustering, which typically differs from the dimensions of the user-level \code{data}, if projections have been applied, in particular). This second rule is controlled by the \code{maxnbCluster} option. 

For large number of points, experience shows that the maximum value derived from these two rules rules is practically always selected by AIC. So, in practice it is faster to only perform clustering with this maximum number of cluster, rather than to perform AIC-based selection among a range of number of clusters. This rule is implemented as the default for argument \code{nbCluster} of \code{\link{refine.default}}, by its default value specified by \code{refine_nbCluster}.
}
\usage{
seq_nbCluster(nr, nc=(nr/500+2)/3)
refine_nbCluster(nr, nc, onlymax=7)
get_nbCluster_range(projdata, nr = nrow(projdata), nc = ncol(projdata), 
                    nbCluster = seq_nbCluster(nr, nc), verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{projdata}{data frame: the data to be clustered, which typically include parameters and \bold{projected} summary statistics;}
  \item{nr}{integer: number of rows of the data to be clustered;}
  \item{onlymax}{integer: see Description;}
  \item{nc}{integer: number of columns of the data to be clustered, typically \bold{twice} the number of estimated parameters (except if latent variables are included);}
  \item{nbCluster}{integer or vector of integers: candidate values, which feasability is checked by the function.}
  \item{verbose}{boolean. Whether to print some information, or not.}
}
\value{ An integer vector}
\examples{
# Determination of number of clusters when attempting to estimate 
#   20 parameters from a reference table with 30000 rows:
seq_nbCluster(nr=30000L)
get_nbCluster_range(nr=30000L, nc=40L) # nc = *twice* the number of parameters
}
\references{
   Rousset F, Raphaël L, Estoup A, Marin J (2026). A new iterative framework for simulation-based population
  genetic inference with improved coverage properties of confidence intervals. \emph{Peer Community Journal}
  \bold{6}: e43. \doi{10.24072/pcjournal.721}\cr
  Supplementary Information: \doi{10.5281/zenodo.19615137}
}
