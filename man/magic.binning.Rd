\name{multi_binning}
\alias{multi_binning}
\title{
Multivariate histogram
}
\description{
Constructs a multivariate histogram of the points. Optionally, first tests whether a given value is within the convex hull of input points and constructs the histogram only if this test is TRUE. This function is available for development purposes but is not required otherwise . It is sparsely documented and subject to changes without notice.  
}
%~ FR->FR doc to write here if the functio is more exposed
\usage{
multi_binning(m, subsize=trunc(nrow(m)^(Infusion.getOption("binningExponent"))),
              expand=5/100, focal=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
A matrix representing points in \emph{d}-dimensional space, where \emph{d} is the number of columns
}
  \item{subsize}{
 A control parameter for an undocumented algorithm
}
  \item{expand}{
 A control parameter for an undocumented algorithm
}
  \item{focal}{
  Value to be tested for inclusion within the convex hull. Its elements must have names.
}
}
\details{
The algorithm may be detailed later.
}
\value{
  Either NULL (if the optional test returned FALSE), or an histogram represented as a data frame each row of which represents an histogram cell by its barycenter (a point in \emph{d}-dimensional space), its \dQuote{binFactor} (the volume of the cell times the total number of observations) and its \dQuote{count} (the number of observations within the cell).
  The returned data frame has the following attributes: \code{attr(.,"stats")} are the column names of the \emph{d}-dimensional points; 
  \code{attr(.,"count")} is the column name of the count, and 
  \code{attr(.,"binFactor")} is the column name of the binFactor.
}
