\name{save_MAFs}
\alias{save_MAFs}
\alias{load_MAFs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Save or load MAF Python objects
}
\description{
\pkg{Infusion} fit results created using the \pkg{mafR} package contain (pointers to) Python objects, which are lost (their pointers being reduced to null pointers) when the fit object is saved on file and reloaded. The functions described here circumvent this issue. 

\code{save_MAFs} will save the Python objects in distinct files. \code{load_MAFs} will load them back into the fit object, \bold{not} overriding anyone pre-existing into the target fit object. 
The Python objects are saved under, and read from, files whose names are made of \code{"jointdens"}, \code{"pardens"}... and the given \code{prefix} and \code{ext}ension.  
}
\usage{
save_MAFs(object, ext="_MAF.pkl", prefix="")
load_MAFs(object, ext="_MAF.pkl", prefix="", set_path_only=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{"SLik_j"} as produced by the up-to-date \pkg{Infusion} workflow.
}
  \item{prefix, ext}{
Character: Prefix and extension for the filename of each saved MAF object.  
}
  \item{set_path_only}{
Boolean: if TRUE, checks the presence of the saved MAF files, but do not load them, and (re)set the path information in the \code{object} (see Details for a programming context of usage).  
}
}
\details{
Both functions can write file path information into the input \code{object}'s \code{load_MAFs_info} element (an environment), where it can be read afterwards. \code{save_MAFs} will do so in all cases, and \code{load_MAFs} will do so when called with argument \code{set_path_only=TRUE}. When a bootstrap is run using parallelisation, the child processes can thus load files using their location information stored in this way. When the files have been moved after the \code{save_MAFs} call, their location information must then be updated using \code{load_MAFs(., set_path_onlTRUE)}.
}
\value{
Both functions return the updated input \code{object}. 
}
\seealso{
\code{reticulate::\link[reticulate]{py_save_object}} is used to save the Python objects.
}
\examples{
## Given object 'slik_j' of class SLik_j
# save_MAFs(slik_j, prefix = "2024Feb30_")
## => '2024Feb30_jointdens_MAF.pkl', etc.
# save(slik_j, file="slik_j.rda") 
## => Objects from Python environments are not saved in the RData file. 

# and later:
# load(file="slik_j.rda") 
# slik_j <- load_MAFs(slik_j, prefix = "2024Feb30_")
## Objects from pkl files are put back at the right place in 'slik_j'. 

}
