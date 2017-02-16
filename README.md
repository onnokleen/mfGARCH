# mfGARCH - a tidyverse-garch package

An R package to estimate multiplicative mixed-frequency GARCH models (GARCH-MIDAS) as proposed in Engle et al. (2013) and is currently under development.

## Highlights
- A comprehensive method for estimating, plotting and forecasting using GARCH-MIDAS models
- Easy to use due to similar syntax as in the regression-lm-command

## Installation
Right now, you need the Rcpp-package to compile it on your own computer.
Development version (GitHub):
```r
# install.packages("devtools")
library(devtools)
install_github("onnokleen/mfGARCH")
```
## ToDo and Roadmap

ToDo
* Define class mfGARCH for generic functions, i.e. plot etc.
* Introduce non-full sample estimation by introducting sample-begin/end variable for forecasting
* Generate table function for list of mfGARCH models
* Hide internal functions
* Add examples in vignette and readme

Roadmap
* 2017-04 Feature-completeness, ToDo done
* 2017-05 Submission to CRAN

Future
* Introduce other weighting schemes
* Employ more than one variable
* Employ more than GJR-GARCH
* Broom-package plug-in
* Allow for different releases for improving real-time forecasting

## History
- January 2017: Initial submit
