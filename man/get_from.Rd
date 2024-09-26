\name{get_from}
\alias{get_from}
\alias{get_from.default}
\alias{get_from.SLik}
\alias{get_from.SLik_j}
\alias{get_from}
\title{
Backward-compatible extractor from summary-likelihood objects
}
\description{
A generic function, whose default method works for \code{list}, and with specific methods for objects inheriting from classes \code{SLik_j} and \code{SLik}. 
}
\usage{
get_from(object, which, ...)

## S3 methods with additional argument(s)
\method{get_from}{SLik}(object, which, raw=FALSE, force=FALSE, ...)
\method{get_from}{SLik_j}(object, which, raw=FALSE, force=FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Any object with a list structure.
}
  \item{which}{
  Character: identifier for the element to be extracted. See Examples for possible values.
}
  \item{raw}{
  Boolean: if TRUE, \code{object[[which]]} is returned, without any particular check of its value. By default, \code{raw} is FALSE and various operations may be performed on the extracted value (see \dQuote{example} below), including optional recomputation if \code{force} is TRUE.
}
  \item{force}{
  Boolean: if TRUE, the extracted element may be computed if it appears to be missing from the \code{object}. This is notably so for \code{which="RMSEs"} or \code{which="par_RMSEs"};  in these cases, the results of the computation are further saved in the original object.
}
  \item{\dots}{
    further arguments passed to or from other methods (currently not used).
}
}
\value{
Will depend on \code{which}, but aims to retain a convenient format backward-compatible with version 1.4.0.
}
\seealso{\code{\link{logLik}}.}
\examples{
  # ##### 0bserved summary statistics 
  # #  (raw data)
  #   get_from(slik, "raw_data") 
  # #  (projected data, with raw ones as attribute, if relevant)
  #   get_from(slik, which="obs") # or which="stat.obs" or "proj_data"
  # 
  # ##### Reference-table information 
  # #  (raw = unprojected):
  #   get_from(slik, "reftable_raw")
  # # : but this is NULL if no projections were performed.
  #
  # # Projected:
  #   get_from(slik, "reftable")
  # # : including all parameters, latent variables, statistics,
  # #   'cumul_iter' (the iteration in which each sample was added),
  # #   and attributes.
  #
  # ##### RMSEs
  # # On any summary-likelihood object 'slik':
  #   get_from(slik, which="par_RMSEs") # matrix
  # # despite <object>$par_RMSEs being an environment if 
  # #  'slik' was created by version > 1.4.0, as then shown by
  #   get_from(slik, which="par_RMSEs", raw=TRUE)
  #
  # # Further, if 
  #   get_from(slik, which="par_RMSEs")
  # # returns NULL because the element is absent from the object, 
  # # then one can force its computation by 
  #   get_from(slik, which="par_RMSEs", force=TRUE)
  # # The result are saved in the 'slik' object, so running again 
  #   get_from(slik, which="par_RMSEs")
  # # will no longer return NULL.
  #
  # ##### Other, less commonly needed elements
  #   get_from(slik, "Simulate"")
  #   get_from(slik, "control.Simulate")
  #   get_from(slik, "env")
  #   get_from(slik, "packages")
  # # => for these elements, the documentation and arguments of refine.default() 
  # #    provides meaning and context where they may be used.
  #
  # Number of elements of the multivariate gaussian mixture model:
  #   get_from(slik, "nbCluster") 
}