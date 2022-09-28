\name{Infusion}
% ----- \Rdversion{1.1}
\alias{Infusion}
\alias{Infusion-package}

\docType{package}

\title{Inference using simulation}

\description{
  Implements a collection of methods to perform inferences based on simulation of realizations of the model considered. In particular it implements 
  \dQuote{summary likelihood}, an approach that effectively evaluates and uses the likelihood of simulated summary statistics. 
%It also implements a \dQuote{summary tail \emph{p}} method that corrects the small sample inaccuracies of likelihood methods.  
}
\details{
  The methods implemented in \code{Infusion} by default assume that the summary statistics have densities. Special values of some statistic, having discrete probability mass, can be handled using the \code{boundaries} attribute of the observed summary statistics (see \code{\link{handling_NAs}} for a further use of this attribute).              
}
\note{See examples \code{\link{example_reftable}} for the most complete example using up-to-date workflow, and \code{\link{example_raw_proj}} or \code{\link{example_raw}} for older workflows.}
\examples{
  ## see Note for links to examples.
}
\keyword{ package }
