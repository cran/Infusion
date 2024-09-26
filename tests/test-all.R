# tools::check_packages_in_dir(".", check_args = "--no-multiarch")
if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## not on CRAN
  if(require("testthat", quietly = TRUE)) {
    pkg   <- "Infusion"
    require(pkg, character.only=TRUE, quietly=TRUE)
    if (interactive())  {
      { ## block of standard tests
        # options(error=recover)
        while (dev.cur()>1L) dev.off()
        oldask <- devAskNewPage(ask=FALSE)
        if (FALSE) {
          testfiles <- dir(paste0(Infusion::projpath(),"/package/tests/testthat/"),pattern="*.R",full.names = TRUE)
        } else {
          testfiles <- dir(paste0(Infusion::projpath(),"/package/tests/testthat/"),full.names = TRUE)
          testfiles <- testfiles[grep("*.R$",testfiles)]
        }
        oldmaxt <- Infusion.options(example_maxtime=100, ## won't test much otherwise
                                    nb_cores=NULL) # enforce default nb_cores
        timings <- t(sapply(testfiles, function(fich){system.time(source(fich))}))
        Infusion.options(oldmaxt)
        print(colSums(timings))
        devAskNewPage(oldask) ## fixme: pb if dev.new has been called in-between 
      }
      if (FALSE) { ## tests not included in package (using unpublished data, etc.)
        if (TRUE) { # see above comment about Rstudio
          priv_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests_private/"),pattern="*.R$",full.names = TRUE)
        } else {
          priv_testfiles <- dir(paste0(spaMM::projpath(),"/package/tests_private/"),full.names = TRUE)
          priv_testfiles <- priv_testfiles[grep("*.R$",priv_testfiles)]
        }
        priv_timings <- t(sapply(priv_testfiles, function(fich){
          cat(cli::col_green(paste0("\n",fich)))
          gc()
          tps <- system.time(chk <- try(source(fich)), gcFirst=FALSE)
          if (inherits(chk,"try-error")) warning(paste0(fich," generated an error"))
          tps
        }))
        print(colSums(priv_timings)) # very roughly 1205.44 s elapsed for default maxtime (0.7) # much less after improving COMP
      }
      if (FALSE) { # slow and using several packages external to Infusion (but faster if one keeps the cache...)
        source(paste0(projpath(),"/../vignette/knitr.R"))
      }
    } else { ## for devtools or R CMD check , cf ?test_check for using library() here:
      if (TRUE) { # condition to characterize devtools::check?
        # this means that the test_all checks are run by R CMD Check (contrary to those in spaMM)
        # and acan detect issues from using global variables, that R CMD Check would not detect on spaMM.
        library("testthat")
        library(pkg, character.only = TRUE)
        report <- test_check(pkg) 
        print(warnings()) # TODO? catch most of these by expect_warning(..)
      }
    }
  } else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
  }
} else cat( "Environment variable _LOCAL_TESTS_ must be set to TRUE for the tests to run.\n" )
