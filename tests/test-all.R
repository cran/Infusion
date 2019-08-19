if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## not on CRAN
  if(require("testthat", quietly = TRUE)) {
    pkg   <- "Infusion"
    require(pkg, character.only=TRUE, quietly=TRUE)
    if (interactive())  {
      while (dev.cur()>1L) dev.off()
      oldask <- devAskNewPage(ask=FALSE)
      testfiles <- dir("C:/home/francois/travail/stats/InferentialSimulation/Infusion/package/tests/testthat/",pattern="*.R",full.names = TRUE)
      oldmaxt <- Infusion.options(example_maxtime=100) ## won't test much otherwise
      timings <- t(sapply(testfiles, function(fich){system.time(source(fich))}))
      Infusion.options(oldmaxt)
      print(colSums(timings))
      devAskNewPage(oldask) ## fixme: pb if dev.new has been called in-between 
    } else {
      test_check(pkg) ## for R CMD check
      print(warnings()) # TODO? catch most of these by expect_warning(..)
    }
  } else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
  }
} else cat( "Environment variable _LOCAL_TESTS_ must be set to TRUE for the tests to run.\n" )