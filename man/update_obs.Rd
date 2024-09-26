\name{.update_obs}
\alias{.update_obs}
\title{
Updating an 'SLik_j' object for new data
}
\description{
\code{.update_obs} is an \emph{experimental} utility to recycle the information in an inference object to produce an inference for different data without repeating a simulation workflow. Beyond its experimental nature, its result is not expected to provide the same precision of inference as a standard iterative \pkg{Infusion} workflow, since the recycled simulations were adapted to the original data only (\code{.update_obs} indeed allows one to investigate the effect of using this non-adapted information). So the new results should not be used as a equivalent to a full \pkg{Infusion} iterative workflow.   
}
\usage{
.update_obs(object, new.obs, CIs = FALSE, eval_RMSEs = FALSE, ...)
}
\arguments{
  \item{object}{
an object of class \code{"SLik_j"} as createded by \code{\link{infer_SLik_joint}} (possible updated by \code{MSL} and \code{refine}).
}
  \item{new.obs}{
A named vector of summary statistics (ideally the projected ones, but raw ones are handled)
}
  \item{CIs, eval_RMSEs, \dots}{
further arguments passed to \code{\link{MSL}}.
}
}
\value{
An object of class \code{"SLik_j"} as the input \code{object}.
}
