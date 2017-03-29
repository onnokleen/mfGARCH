library(devtools)
library(roxygen2)


#' @importFrom Rcpp sourceCpp
#' @useDynLib mfGARCH
NULL

devtools::use_rcpp()
Rcpp::compileAttributes()

devtools::document()
devtools::load_all()

# Clean up empty spaces
#formatR::tidy_dir("R")
# Check for inconsistencies
lintr::lint_package()

