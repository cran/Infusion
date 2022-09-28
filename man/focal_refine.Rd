\name{focal_refine}
\alias{focal_refine}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Refine summary likelihood profile in focal parameter values
}
\description{
This function refines an \code{SLik_j} object in a focused way defined by focal parameter values. It is paticularly useful to check a suspect pattern in a likelihood profile. If there is a suspect dip or peak at value <somepar>=<somevalue>, \code{focal_refine(}<SLik_j object>\code{, focal}=c(<somepar>=<somevalue>)\code{, size}=<size>) will define <size> parameter points near \code{c(}<somepar>=<somevalue>) and (subject to these points being in the parameter bounds of the object) simulate new samples for these parameter points and refine the object using these new simulations. 
}
\usage{
focal_refine(object, focal, size, plotprof = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  An object of class \code{SLik_j}.
}
  \item{focal}{
Parameter value(s) (as a vector of named values) 
}
  \item{size}{
Target number of points to add to the reference table
}
  \item{plotprof}{
Whether to replot a likelihood profile (1D or 2D depending on the dimension of \code{focal}).
}
  \item{\dots}{
Further arguments passed to \code{profile.SLik_j} (but not including argument \code{return.optim}) or \code{refine}.
}
}
\value{
The updated \code{object}
}
\examples{
\dontrun{
# Using the slik_j object from the toy example in help("example_reftable"):

plot1Dprof(slik_j,"s2")
slik_fix <- focal_refine(slik_j,focal=c(s2=2), size=100)
plot1Dprof(slik_fix,"s2")

# In that case the effect is not spectacular because 
# there is no major problem in the starting profile. 
}
}

