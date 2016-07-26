\name{options}
\alias{Infusion.options}
\alias{Infusion.getOption}
\title{Infusion options settings}
\description{
  Allow the user to set and examine a variety of \emph{options}
  which affect operations of the Infusion package. 
  However, typically these should not be modified, and if they are, not more than once in a data analysis.
}
\usage{
Infusion.options(\dots)

Infusion.getOption(x)

}
\arguments{
  \item{x}{a character string holding an option name.}
  \item{\dots}{A named value or a list of named values. The following values, with their defaults, 
    are used in \code{Infusion}:

  \describe{
\item{\code{projTrainingSize = 200}:}{default value of \code{trainingsize} argument of \code{\link{project.character}}.}
\item{\code{projKnotNbr = 300}:}{default value of \code{knotnbr} argument of \code{\link{project.character}}.}
\item{\code{logLname = "logL"}:}{default value of \code{logLname} argument of \code{\link{infer_logLs}}. The name given to the inferred log likelihoods in all analyses.}
\item{\code{LRthreshold= - qchisq(0.999,df=1)/2}:}{A value used internally by \code{\link{sample_volume}} to sample points 
      in the upper region of the likelihood surface, as defined by the given likelihood ratio threshold.
}
\item{\code{precision = 0.1}:}{default value of \code{precision} argument of \code{\link{refine}}. Targets RMSE of log L and log LR estimates.}
\item{\code{nRealizations=1000}:}{default value of \code{nRealizations} argument of \code{\link{add_simulation}}. Number of realizations for each empirical distribution.}
\item{\code{mixmodGaussianModel="Gaussian_pk_Lk_Dk_A_Dk"}:}{default models used in clustering by \code{Rmixmod}. Run \code{Rmixmod::mixmodGaussianModel()} for a list of possible models, and see the statistical documentation (Mixmod Team 2016) for explanations about them.}
\item{\code{nbCluster = quote(ceiling(nrow(data)^0.3))}:}{default value of \code{nbCluster} used in clustering by \code{Rmixmod}}
\item{\code{example_maxtime=2.5}:}{Used in the documentation to control whether the longer examples should be run. 
The approximate running time of given examples (or some very rough approximation for it) on one author's laptop is compared to this value.} 
} 
and possibly other undocumented values for development purposes.
  }

}
\details{
  Invoking \code{Infusion.options()} with no arguments returns a list with the
  current values of the options. Invoking \code{Infusion.getOption(<option name>)} returns the value of the option rather than a list.
  
  The default \code{nbCluster} value is the upper value of the range recommended in the \code{mixmod} statistical documentation (Mixmod Team, 2016). If clustering by the given number(s) of clusters fails, decreasing values are tryied until success.
}
\references{
Mixmod Team (2016). Mixmod Statistical Documentation. Université de Franche-Comté,
Besançon, France. Version: February 10, 2016 retrieved from \url{http://www.mixmod.org}.
}
\value{
  For \code{Infusion.getOption}, the current value set for option \code{x}, or
  \code{NULL} if the option is unset.

  For \code{Infusion.options()}, a list of all set options sorted by category.  For
  \code{Infusion.options(name)}, a list of length one containing the set value,
  or \code{NULL} if it is unset.  For uses setting one or more options,
  a list with the previous values of the options changed (returned
  invisibly).
}
\examples{
  Infusion.options()
  Infusion.getOption("LRthreshold")
  \dontrun{
  Infusion.options(LRthreshold=- qchisq(0.99,df=1)/2)
  }
}
