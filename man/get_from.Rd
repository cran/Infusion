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
A generic function, whose default method works for list, and with specific methods for objects inheriting from classes \code{SLik_j} and \code{SLik}. 
}
\usage{
get_from(object, which, ...)

## S3 methods with additional argument(s)
\method{get_from}{SLik}(object, which, raw, ...)
\method{get_from}{SLik_j}(object, which, raw, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  Any object with a list structure.
}
  \item{which}{
  Character: names of element to be extracted.
}
  \item{raw}{
  Boolean: if TRUE, \code{object[[which]]} is returned, ignoring any more specific processing.
}
  \item{\dots}{
    further arguments passed to or from other methods (currently not used).
}
}
\value{
Will depend on \code{which}, but aims to retain a convenient format backward compatible with version 1.4.0.
}
\seealso{\code{\link{logLik}}.}
\examples{
  # 0bserved summary statistics 
  #  (projected, with raw ones as attribute, if relevant)
  # get_from(slik, "obs") 
  # 
  # On any summary-likelihood object 'slik':
  # get_from(slik, which="par_RMSEs") # matrix
  # despite <object>$par_RMSEs being an environment if 
  #  'slik' was created by version > 1.4.0, as then shown by
  # get_from(slik, which="par_RMSEs", raw=TRUE)
}