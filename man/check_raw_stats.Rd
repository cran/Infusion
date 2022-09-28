\name{check_raw_stats}
\alias{check_raw_stats}
\title{
Check linear dependencies among raw summary statistics
}
\description{
A convenient wrapper function for \code{caret::findLinearCombos}, allowing to detect linear dependencies among the statistics, 
and optionally to remove variables that induce them.
}
\usage{
check_raw_stats(x, statNames, remove = FALSE, verbose = interactive())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data frame (particularly inheriting from class \code{"reftable"}, i.e. a reference table of simulations); or possibly a matrix with column names}
  \item{statNames}{
Character vector: variables among which dependencies are sought. Must belong column names of \code{x}. For a \code{reftable}, this argument is optional and by default, all raw statistic are included. For other classes of input, this argument is required. 
}
  \item{remove}{Boolean: whether to return \code{x} with \dQuote{offending} columns removed, or other information.}
  \item{verbose}{Boolean: whether to display some messages.}
}
\value{
Return type depends on the availability of the \pkg{caret} package, and on the \code{remove} argument, as follows. if \code{remove=TRUE}, an object of the same class as \code{x} is returned (with redundant columns removed). If \code{remove=FALSE}, either the \pkg{caret} package is available, in which case a list is returned with the same structure as the return value of  \code{caret::findLinearCombos} but with column indices replaced by column names; or a message pointing that \pkg{caret} is not available is returned (and another is printed, only once per session).
}
