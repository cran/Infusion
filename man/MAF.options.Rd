\name{MAF.options}
\alias{config_mafR}
\alias{MAF.options}
\title{
Control of MAF design and training
}
\description{
Masked Autoregressive Flows can be used by \pkg{Infusion} to infer various densities.
This functionality requires the \pkg{mafR} package, and is requested through the \code{using} argument of 
\code{infer.SLik.joint} or \code{refine} (see Details).

\code{config_mafR} is a convenient way to reset the Python session (erasing all results stored in its memory), and in particular to enforce GPU usage.

\code{MAF.options} is a wrapper for \code{Infusion.options} which facilitates the specification of 
options that control the design of Masked Autoregressive Flows and their training.
}
\usage{
config_mafR(torch_device, ...)
MAF.options(template = "zuko-like", ...)
}
\arguments{
  \item{torch_device}{
character: \code{"cpu"}, \code{"cuda"} (or \code{"cuda:0"}, etc., useful if several GPUs are available) or \code{"mps"}.
}
  \item{template}{
Defines a template list of options, subsequently modified by options specified through the \dots, if any.
Possible values of \code{template} are \code{"zuko-like"}, \code{"PSM19"}, \code{NULL}, or a \code{list} of options. The default value \code{"zuko-like"} reproduces \pkg{Infusion}'s defaults. 
}
  \item{\dots}{For \code{MAF.options}: a named value, or several of them, that override or complete the \code{template} values. For \code{config_mafR}: for development purposes; not documented.
}
}
\details{
\bold{Possible \code{using} values:}

With \code{using="c.mafR"} four different MAFs are computed in every iteration:
the density of statistics given parameters (providing the likelihood), the joint density, the instrumental density of parameters, and the \dQuote{instrumental posterior} density of parameters given the data.

With \code{using="MAFmix"}, MAFs are computed only for the joint density and the instrumental density. The likelihood is deduced from them and a multivariate Gaussian mixture model is used to infer the \dQuote{instrumental posterior} density.

\code{using="mafR"} can be used to let \pkg{Infusion} select one of the above options. \code{"MAFmix"} is currently called as it is faster, but this is liable to change in the future, so do not use this if you want a repeatable workflow.

\bold{Possible \code{template} values for \code{MAF.options}:}

\code{"PSM19"} defines a template list of control values recommended by Papamakarios et al. (2019, Section 5.1).
\code{"zuko-like"} defines a template list of values inspired from the tutorials of the \pkg{zuko} Python package, with some extra twists. Concretely, their definition is
\preformatted{
    if (template == "zuko-like") {
        optns <- list(design_hidden_layers = .design_hidden_layers_MGM_like, 
            MAF_patience = 30, MAF_auto_layers = 3L, Adam_learning_rate = 0.001)
    }
    else if (template == "PSM19") {
        optns <- list(design_hidden_layers = .design_hidden_layers_PSM19, 
            MAF_patience = 20, MAF_auto_layers = 5L, Adam_learning_rate = 1e-04)
    }
}
and any replacement value should match the types (function, numeric, integer...) of the shown values, and the formals of the internal functions, in order to avoid cryptic errors. 

The internal \code{.design_hidden_layers...} functions return a vector of numbers \eqn{H_i} of hidden values per layer \eqn{i} of the neural network. THevector has an attribute giving the resulting approximate number of parameters \eqn{P} of the deep-learning model according to Supplementary Table 1 of Papamakarios et al. 2017. \eqn{H_i=}\code{50L} for \code{"PSM19"}. For \code{"zuko-like"}, a typically higher value will be used. It is defined as a power of two such that \eqn{P} is of the order of 8 to 16 times the default number of parameters of the multivariate gaussian mixture model that could be used instead of MAFs.     


Other controls which can be modified through the \dots are \cr
\code{ * MAF_validasize}, a function which returns the size of the validation set,
whose default definition returns 5\% of its input value \code{nr} which is the number of samples in the reference table (consistently with Papamakarios et al., 2019);\cr
\code{ * MAF_batchsize}, a function that returns the batch size for the Adam optimizer. Its default simply returns \code{100L}, but non-default functions can be defined, with at least the \dots as formal arguments (more elaborate formals are possible but not part of the API).
}
\value{
\code{config_mafR} is used for its side effects. Returns NULL invisibly.

\code{MAF.options} returns the result of calling \code{Infusion.options} on the arguments defined by the \code{template} and the \dots . Hence, it is a list of \bold{previous} values of the affected options. 
}
\references{
Papamakarios, G., T. Pavlakou, and I. Murray (2017) Masked Autoregressive
Flow for Density Estimation. Pp. 2335–2344 in I. Guyon, U. V. Luxburg,
S. Bengio, H. Wallach, R. Fergus, S. Vishwanathan, and R. Garnett, eds.
Advances in Neural Information Processing Systems 30. Curran Associates,
Inc.\cr
http://papers.nips.cc/paper/6828-masked-autoregressive-flow-for-density-estimation

Rozet, F., Divo, F., Schnake, S (2023) Zuko: Normalizing flows in PyTorch. 
https://doi.org/10.5281/zenodo.7625672
}
\seealso{\code{\link{save_MAFs}} for saving and loading MAF objects.}
\examples{
  MAF.options(template = "zuko-like", 
              Adam_learning_rate=1e-4,
              MAF_batchsize = function(...) 100L)

\dontrun{
## These examples require the mafR package, 
## and a Python installation with cuda capabilities.

if (requireNamespace("mafR", quietly=TRUE)) {

config_mafR() # set default: "cpu", i.e. GPU not used
config_mafR("cuda") # sets cuda as GPU backend 
config_mafR() # reset python session, keeping latest backend


  config_mafR(torch_device="cuda")
  
  # function for sampling from N(.,sd=1)
  toyrnorm <- function(mu) {
    sam <- rnorm(n=40,mean=mu,sd=1)
    return(c(mean1=mean(sam)))
  } 

  # simulated data, standing for the actual data to be analyzed:  
  set.seed(123)
  Sobs <- toyrnorm(mu=4) 

  parsp <- init_reftable(lower=c(mu=2.8), 
                         upper=c(mu=5.2))
  simuls <- add_reftable(Simulate="toyrnorm", parsTable=parsp)

  MAF.options(template = "zuko-like")
  densv <- infer_SLik_joint(simuls,stat.obs=Sobs, using="mafR")

  # Usual workflow using inferred surface:
  slik_j <- MSL(densv, eval_RMSEs = FALSE) ## find the maximum of the log-likelihood surface
  # ETC.
  
  save_MAFs(slik_j, prefix = "toy_") # See its distinct documentation.
}
}
}
