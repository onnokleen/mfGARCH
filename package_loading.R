library("devtools")
library(roxygen2)

devtools::document()
devtools::load_all()

# Clean up empty spaces
formatR::tidy_dir("R")
# Check for inconsistencies
lintr::lint_package()

