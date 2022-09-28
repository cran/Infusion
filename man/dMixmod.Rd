\name{dMixmod}
\alias{dMixmod}
\alias{class:dMixmod}
\alias{dMixmod-class}
\alias{NULLorChar}
\alias{class:NULLorChar}
\alias{NULLorChar-class}
\alias{NULLorNum}
\alias{class:NULLorNum}
\alias{NULLorNum-class}
\alias{plot.dMixmod}
\docType{methods}
\title{
  Internal S4 classes. 
}
\description{
  The objects or methods referenced here are not to be called by the user, or are waiting for documentation to be written.  

\code{dMixmod} is an S4 class describing some distributions that extend the multivariate gaussian mixture models (MGMM) by possibly involving discrete probability masses for some variables and gaussian mixtures for other variables conditional on such discrete events. In terms of the represented probability models, and of its slots, is effectively extends the \code{MixmodResults} class from the \code{Rmixmod} package. But it does not formally extends this class in terms of OOP programming. It should not be considered as part of the programming interface, and may be subject to backward-incompatible modifications without notice. In the current implementation it cannot represent general mixtures of discrete probabilities and MGMMs, and may yield correct results only for the degenerate case of pure MGMMs or when inference can be based on the conditional density of continuous variables conditional on the (joint-, if relevant) discrete event observed in the data.
}
\usage{
# dMixmod: Don't try to use it! It's for programming only.
}
\value{
A \code{dMixmod} object has the same slots as a \code{MixmodResults} object, plus additional ones: \code{@freq} is the frequency of the conditioning event for the gaussian mixture model. In the \code{Infusion} code, this event is defined jointly by the \dQuote{observed} summary statistics and the reference simulation table: a probability mass for specific values \bold{v} is identified from the simulated distribution of summary statistics in the reference table, and \code{freq} is an estimate of the probability mass if the summary statistics match \bold{v}, or the converse probability if they do not match.
}
\note{
Use \code{str(attributes(.))} to see the slots of a \code{dMixmod} object if \code{str(.)} does not work.
}
\examples{
# The dMixmod object can be used internally to handle repeated and boundary values 
# of summary statistics. The user has to add an attribute to the observations,
# as explained in help("boundaries-attribute"):
Sobs <- c(mean=4.321, se=0.987) # hypothetical observation
attr(Sobs,"boundaries") <- c(someSummStat=-1)  
}