\name{plot.SLik}
\alias{plot.SLik}
\alias{plot.SLik_j}
\alias{plot.SLikp}
\title{
Plot SLik or SLikp objects
}
\usage{
\method{plot}{SLik}(x, y, filled = FALSE, decorations = NULL,
                    color.palette = NULL, plot.axes = NULL, 
                    plot.title = NULL, ...)
\method{plot}{SLik_j}(x, y, filled = nrow(x$logLs)>5000L, decorations = NULL, 
                      color.palette = NULL, plot.axes = NULL, 
                      plot.title = NULL, from_refine=FALSE, ...)
%\method{plot}{SLikp}(x, y, filled=FALSE, log.=TRUE,...)
}
\description{
Mostly conceived for exposition purposes, for the two-parameters case.
The black-filled points are those for which the observed summary statistic was outside of the convex hull of the simulated empirical distribution. The crosses mark the estimated ML point and the confidence intervals points, that is, the outmost points on the contour defined by the profile likelihood threshold for the profile confidence intervals. There is a pair of CI points for each interval.
The smaller black dots mark points added in the latest iteration, if \code{refine} was used.
}
\arguments{
  \item{x}{
  An object of class \code{SLik} or \code{SLikp}
}
  \item{y}{
  Not used, but included for consistency with the \code{plot} generic. 
}
  \item{filled}{
  whether to plot a \code{\link[spaMM]{mapMM}} or a \code{\link[spaMM]{filled.mapMM}}.
}
  \item{decorations}{
  Graphic directives added to the default \code{decorations} value in calls of \code{\link[spaMM]{mapMM}} or a \code{\link[spaMM]{filled.mapMM}} (see the source code of \code{plot.SLik} for the latter default values). 
}
%  \item{log.}{
%  Whether to log transform the response variable (i.e, a frequency). Zero frequency values are %converted to half the minimum nonzero frequency value before log transformation. 
%}
  \item{color.palette}{
  Either NULL or a function that can replace the default color function used by \code{plot.SLik}. Thefunction must have a single argument, giving the number of color levels. 
}
\item{plot.title}{	
statements which replace the default titles to the main plot (see Details). 
}
\item{plot.axes}{	
statements which replace the default axes on the main plot (see Details). 
}
\item{from_refine}{	
For programming purposes, not documented. 
}
  \item{\dots}{
further arguments passed to or from other methods (currently can be used to pass a few arguments such as \code{map.asp} in all cases, 
or \code{variances} to \code{filled.mapMM}).
}
}
\value{Returns the plotted object invisibly.}
\details{Different graphic functions are called depending on the number of estimated parameters. For two parameters, \code{\link[spaMM]{mapMM}} or \code{\link[spaMM]{filled.mapMM}} are called. For more than two parameters, \code{\link[spaMM]{spaMM.filled.contour}} is called. See the documentation of these functions for the appropriate format of the \code{plot.title} and \code{plot.axes} arguments.  

}
\examples{
\dontrun{
## Using 'slik' object from the example in help("Infusion-package") 
  plot(slik, filled=TRUE,
       plot.title=quote(title("Summary-likelihood-ratio surface",
                              xlab=expression(mu),
                              ylab=expression(sigma^2))))
}
}
