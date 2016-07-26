if (Sys.getenv("_LOCAL_TESTS_")=="TRUE") { ## not on CRAN
  if(require("testthat", quietly = TRUE)) {
    pkg   <- "Infusion"
    require(pkg, character.only=TRUE, quietly=TRUE)
    test_check(pkg) ## for R CMD check
    print(warnings()) # TODO? catch most of these by expect_warning(..)
  } else {
    cat( "package 'testthat' not available, cannot run unit tests\n" )
  }
} else cat( "Environment variable _LOCAL_TESTS_ must be set to TRUE for the tests to run.\n" )