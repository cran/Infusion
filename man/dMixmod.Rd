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
\docType{methods}
\title{
  Internal S4 classes. 
}
\description{
  The objects or methods referenced here are not to be called by the user, or are waiting for documentation to be written.  

\code{dMixmod} is and S4 class describing distributions that involve discrete probability masses for some variables and gaussian mixtures for other variables conditional on such discrete events. In terms of the represented probability models, and of its slots, is effectively extends the \code{MixmodResults} class from the \code{Rmixmod} package. But it does not formally extends this class in termes of OOP programming. It should not be considered as part of the programming interface, and may be subject to backward-incompatible modifications without notice.
}
\usage{
# dMixmod: Don't try to use it!
}
\note{
Use \code{str(attributes(.))} to see the slots of a \code{dMixmod} object if \code{str(.)} does not work.
}
\value{
A \code{dMixmod} object has the same slots as a \code{MixmodResults} object, plus additional ones: \code{@freq} is the frequency of the conditioning event for the gaussian mixture model. In the \code{Infusion} code, this event is defined jointly by the \dQuote{observed} summary statistics and the reference simulation table: a probability mass for specific values \bold{v} is identified from the simulated distribution of summary statistics in the reference table, and \code{freq} is an estimate of the probability mass if the summary statistics match \bold{v}, or the converse probability if they do not match.
}