\name{densv}
\alias{densv}
\alias{densb}
\alias{saved_seed}
\docType{data}
\title{
Saved computations of inferred log-likelihoods
}
\description{
These are saved results from toy examples used in other documentation page for the package. They give estimates by simulation of log-likelihoods of the \code{(mu,s2)} parameters of a Gaussian distribution for a given sample of size 20 with mean 4.1416238 and (bias-corrected) variance 0.9460778. \code{densv} is based on the sample mean and sample variance as summary statistics, and \code{densb} on more contrived summary statistics.   
}
\usage{
data("densv")
data("densb")
}
\format{
  Data frames (with additional attributes) with observations on the following 5 variables.
  \describe{
    \item{\code{mu}}{a numeric vector; mean parameter of simulated Gaussian samples}
    \item{\code{s2}}{a numeric vector; variance parameter of simulated Gaussian samples}
    \item{\code{sample.size}}{a numeric vector; size of simulated Gaussian samples}
    \item{\code{logL}}{a numeric vector; log probability density of a given statistic vector inferred from simulated values for the given parameters}
    \item{\code{isValid}}{a boolean vector. See \code{\link{infer_logLs}} for its meaning.} 
  }
  Both data frames are return objects of a call to \code{\link{infer_logLs}}, and as such they includes attributes providing information about the parameter names and statistics names (not detailed here).
}
\seealso{
 See step (3) of the workflow in the Example on the main \code{\link{Infusion}} documentation page, showing how \code{densv} was produced, and the Example in \code{\link{project}} showing how \code{densb} was produced.
}
\keyword{datasets}
