\name{plot_proj}
\alias{plot_proj}
\alias{plot_importance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic plots for projections 
}
\description{
\code{plot_proj} and \code{plot_importance} are convenience functions providing diagnostic plots for projections, \code{plot_proj} plots predictions. This function is tailored for use on \pkg{ranger} projections, with out-of bag predictions for the training set shown in blue, and predictions for points outside the training set shown in black. Either \code{parm} or \code{proj} will then be needed to identify the projection to be plotted (perhaps both in some programming contexts, or if other projections methods are used). 

\code{plot_importance} plots the importance metric \code{variable.importance}  stored in a \code{ranger} object. Here, the projection may be identified by either \code{parm} or \code{proj}, or even directly provided as the \code{object}.
}
\details{
Projectors may be updated in an object after the projected statistics have been computed and included in the reference table, without any explicit action of the user on the object (see general information about \pkg{\link{Infusion}} for why and when this may occur). In particular the projectors stored in input and output fit objects of a \code{refine} call are stored in the same environment, and therefore the projectors stored in the input object are modified in light of new simulations. In that case, the plots produced by \code{plot_proj} may reflect properties either of the updated projections (if \code{new_rawdata} is used), or worse, may mix results of different projections (if \code{new_rawdata} is not used). A message or a warning may be issued when such events occur.   
}
\usage{
plot_proj(
  object, parm=NULL, proj, 
  new_rawdata, 
  use_oob=Infusion.getOption("use_oob"), is_trainset=FALSE, 
  xlab=NULL, ylab=NULL, ...) 

plot_importance(object, parm, proj, n.var = 30L, xlim=NULL, 
                xlab = "Variable Importance", ylab = "", main="", ...)
}
\arguments{
  \item{object}{An object of class \code{SLik_j}. For \code{plot_importance}, it may also be an object of class \code{ranger}.}
  \item{parm}{Character string: a parameter name. Either one of \code{parm} or \code{proj} is required to get a diagnostic plot for a specific projector. Otherwise, a multipanel plot will be produced for all projectors.}
  \item{proj}{Character string: name of projected statistic.}
  \item{new_rawdata}{Reference table on which projections should be computed. If NULL,
    projections are not recomputed. Instead,  the values stored in the \code{object} are used.}
  \item{use_oob, is_trainset}{ Passed to \code{\link{project.default}}. Ignored if \code{new_rawdata} is NULL.}
  \item{n.var}{Integer: (maximum) number of predictor variables to be included in the plot.}
  \item{xlim, xlab, ylab, main, \dots}{Passed to \code{plot} (for \code{plot_proj}) or to \code{\link[graphics]{dotchart}}  (for \code{plot_importance}).}
}
\value{
These functions are mainly used for their side effect (the plot). \code{plot_importance} returns the vector of importance values invisibly. \code{plot_proj} returns invisibly, for each parameter a list with the \code{x},\code{y} \code{xlab} and \code{ylab} elements of the \code{plot} call, or a structured list of such lists when plots are produced for several parameters.
}
\note{See workflow examples in \code{\link{example_reftable}}.}
\examples{
  ## see Note for links to examples.
}
